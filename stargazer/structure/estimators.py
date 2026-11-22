# Adapted from SeaQuest (same author, commit 064d667): the plug-in entropy and mutual-information code in
# seaquest/structure/mi.py (`_mi`) and the H1 - H2 computation in seaquest/decipher/structure.py
# (`pareidolia_controls`). Changes: bias-corrected entropy estimators (Miller-Madow, Grassberger, Chao-Shen) replace
# the plug-in estimate; block and pair counting is vectorised with integer codes and np.bincount / np.unique instead
# of Python Counters; conditional entropies go to any order with an explicit undersampling rule; lagged mutual
# information takes an optional time axis with lag bins in time units; compression-based description lengths added.
"""Entropy, block-entropy, conditional-entropy, lagged mutual-information and compression estimators.

Every function takes a one-dimensional integer array `seq` with values in [0, k). Entropies and mutual
information are in bits; compressed sizes are in bits.
"""
from __future__ import annotations

import bz2
import lzma
import math
import zlib
from typing import Iterator, Sequence

import numpy as np
from scipy.special import digamma

__all__ = [
    "ENTROPY_METHODS", "block_counts", "block_entropy", "compression_bits", "compression_ratio", "conditional_entropy",
    "entropy", "entropy_chao_shen", "entropy_grassberger", "entropy_miller_madow", "entropy_plugin",
    "entropy_rate_curve", "lag_pair_index", "lagged_mutual_information", "max_trustworthy_order", "shuffled_mi_null",
]

_LN2 = math.log(2.0)
_DENSE_LIMIT = 1 << 22          # code spaces up to this size are counted with np.bincount
_CODE_LIMIT = 1 << 62           # beyond this, block codes are re-ranked so they never overflow int64
ENTROPY_METHODS = ("plugin", "miller_madow", "grassberger", "chao_shen")


# ----------------------------------------------------------------------------------------------------------------
# entropy estimators on counts
# ----------------------------------------------------------------------------------------------------------------

def _clean_counts(counts: np.ndarray | Sequence[int]) -> tuple[np.ndarray, int]:
    c = np.asarray(counts, dtype=np.int64).ravel()
    if c.size and c.min() < 0:
        raise ValueError("counts must be non-negative")
    c = c[c > 0]
    return c, int(c.sum())


def entropy_plugin(counts: np.ndarray | Sequence[int]) -> float:
    """Maximum-likelihood ("plug-in") entropy in bits, -sum p log2 p with p = n_i / N.

    Biased low by about (m - 1) / (2 N ln 2) bits, m the number of occupied bins. Returns 0 for empty counts.
    """
    c, n = _clean_counts(counts)
    if n == 0:
        return 0.0
    p = c / n
    return float(-(p * np.log2(p)).sum())


def entropy_miller_madow(counts: np.ndarray | Sequence[int]) -> float:
    """Miller-Madow entropy in bits: the plug-in estimate plus (m - 1) / (2 N) nats, m = occupied bins."""
    c, n = _clean_counts(counts)
    if n == 0:
        return 0.0
    return entropy_plugin(c) + (c.size - 1) / (2.0 * n * _LN2)


def _grassberger_g(n: np.ndarray) -> np.ndarray:
    """G(n) = psi(n) + (-1)^n * integral_0^1 t^(n-1) / (1 + t) dt, with the integral written as
    (psi((n + 1) / 2) - psi(n / 2)) / 2. G(1) = -gamma - ln 2, G(2) = 2 - gamma - ln 2, G(2m+1) = G(2m)."""
    n = n.astype(float)
    sign = np.where(n.astype(np.int64) % 2 == 0, 1.0, -1.0)
    return digamma(n) + sign * 0.5 * (digamma((n + 1.0) / 2.0) - digamma(n / 2.0))


def entropy_grassberger(counts: np.ndarray | Sequence[int]) -> float:
    """Grassberger entropy estimate in bits.

    Implements the estimator of P. Grassberger, "Entropy estimates from insufficient samplings",
    arXiv:physics/0307138 (2003, revised 2008), eq. (35) with the G_n of eq. (30):

        H = ln N - (1 / N) * sum_i n_i G(n_i),   G(n) = psi(n) + (-1)^n * int_0^1 t^(n-1) / (1 + t) dt,

    converted from nats to bits. It is nearly unbiased while most occupied bins hold a few counts; it is still
    biased low when the number of possible bins greatly exceeds N (see `entropy_rate_curve` for the rule used here).
    A histogram with a single occupied bin returns 0 (the formula's parity term would give a few hundredths of a
    bit either side of 0).
    """
    c, n = _clean_counts(counts)
    if n == 0 or c.size == 1:
        return 0.0
    values, mult = np.unique(c, return_counts=True)
    h = math.log(n) - float((mult * values * _grassberger_g(values)).sum()) / n
    return max(h, 0.0) / _LN2


def entropy_chao_shen(counts: np.ndarray | Sequence[int]) -> float:
    """Chao-Shen coverage-adjusted entropy in bits (Chao & Shen 2003, Environ. Ecol. Stat. 10, 429).

    Coverage C = 1 - f1 / N (f1 = singletons; f1 = N is replaced by N - 1 so C > 0), p_i = C n_i / N, and
    H = -sum p_i ln p_i / (1 - (1 - p_i)^N), the Horvitz-Thompson correction for unseen bins.
    """
    c, n = _clean_counts(counts)
    if n == 0:
        return 0.0
    f1 = int((c == 1).sum())
    if f1 == n:
        f1 = n - 1
    coverage = 1.0 - f1 / n
    pa = coverage * c / n
    denom = 1.0 - (1.0 - pa) ** n
    return float(-(pa * np.log(pa) / denom).sum()) / _LN2


_ESTIMATORS = {"plugin": entropy_plugin, "miller_madow": entropy_miller_madow, "grassberger": entropy_grassberger,
               "chao_shen": entropy_chao_shen}


def entropy(counts: np.ndarray | Sequence[int], method: str = "grassberger") -> float:
    """Entropy in bits of a histogram by one of "plugin", "miller_madow", "grassberger" (default) or "chao_shen"."""
    try:
        return _ESTIMATORS[method](counts)
    except KeyError:
        raise ValueError(f"unknown entropy method {method!r}; choose from {ENTROPY_METHODS}") from None


# ----------------------------------------------------------------------------------------------------------------
# block counting
# ----------------------------------------------------------------------------------------------------------------

def _as_symbols(seq: np.ndarray | Sequence[int], k: int) -> np.ndarray:
    s = np.asarray(seq).ravel()
    if s.size and not np.issubdtype(s.dtype, np.integer):
        if not np.all(np.equal(np.mod(s, 1), 0)):
            raise ValueError("seq must hold integer symbols")
    s = s.astype(np.int64, copy=False)
    if k < 1:
        raise ValueError("alphabet size k must be at least 1")
    if s.size and (s.min() < 0 or s.max() >= k):
        raise ValueError("symbols must lie in [0, k)")
    return s


def _count_codes(codes: np.ndarray, span: int) -> np.ndarray:
    """Counts of the distinct values in `codes` (all < span), zeros dropped, ordered by code."""
    if codes.size == 0:
        return np.zeros(0, dtype=np.int64)
    if span <= max(_DENSE_LIMIT, 4 * codes.size):
        c = np.bincount(codes, minlength=0)
        return c[c > 0]
    return np.unique(codes, return_counts=True)[1]


def _iter_block_codes(seq: np.ndarray, k: int, n_max: int) -> Iterator[tuple[int, np.ndarray, int]]:
    """Yield (n, codes, span) for n = 1..n_max: integer codes of all overlapping n-blocks, each below `span`.

    Codes are the base-k value of the block while that fits in int64; once k**n would overflow, the codes are
    replaced by their ranks (same equality structure), so counting stays exact for any alphabet and length.
    """
    if n_max < 1 or seq.size == 0:
        return
    codes, span = seq.copy(), k
    yield 1, codes, span
    for n in range(2, min(n_max, seq.size) + 1):
        if span * k > _CODE_LIMIT:
            uniq, codes = np.unique(codes, return_inverse=True)
            codes, span = codes.astype(np.int64).ravel(), int(uniq.size)
        codes = codes[:-1] * k + seq[n - 1:]
        span *= k
        yield n, codes, span


def block_counts(seq: np.ndarray | Sequence[int], k: int, n: int) -> np.ndarray:
    """Counts of every observed (overlapping) n-block of `seq`.

    Inputs: seq, integer symbols in [0, k); k, alphabet size; n, block length.
    Output: int64 array of the non-zero counts (ordered by block code), summing to len(seq) - n + 1.
    n = 0 returns [1] (one empty block); n > len(seq) returns an empty array.
    """
    s = _as_symbols(seq, k)
    if n < 0:
        raise ValueError("block length n must be non-negative")
    if n == 0:
        return np.ones(1, dtype=np.int64)
    out = np.zeros(0, dtype=np.int64)
    for m, codes, span in _iter_block_codes(s, k, n):
        if m == n:
            out = _count_codes(codes, span)
    return out


def block_entropy(seq: np.ndarray | Sequence[int], k: int, n: int, method: str = "grassberger") -> float:
    """Entropy in bits of the distribution of overlapping n-blocks (H(0) = 0)."""
    if n == 0:
        return 0.0
    return entropy(block_counts(seq, k, n), method)


def conditional_entropy(seq: np.ndarray | Sequence[int], k: int, order: int, method: str = "grassberger") -> float:
    """Conditional entropy in bits of a symbol given the `order` preceding symbols, H(order + 1) - H(order).

    The difference of block entropies (rather than an average of per-context entropies) is used because it is the
    standard finite-order entropy-rate estimate h_n = H_{n+1} - H_n and because the bias corrections then apply to
    two histograms only. Order 0 is the plain symbol entropy. No undersampling guard: see `entropy_rate_curve`.
    """
    if order < 0:
        raise ValueError("order must be non-negative")
    return block_entropy(seq, k, order + 1, method) - block_entropy(seq, k, order, method)


def max_trustworthy_order(n_symbols: int, k: int, max_order: int) -> int:
    """Largest order <= max_order whose number of possible contexts k**order does not exceed n_symbols / 5."""
    best = 0
    for order in range(max_order + 1):
        if k ** order * 5 <= n_symbols:
            best = order
        else:
            break
    return best


def entropy_rate_curve(seq: np.ndarray | Sequence[int], k: int, max_order: int,
                       method: str = "grassberger") -> np.ndarray:
    """Conditional entropies h_0..h_max_order in bits (h_m = H(m + 1) - H(m)).

    Undersampling rule: h_m is reported only while the number of possible contexts k**m is at most len(seq) / 5,
    i.e. while every possible context could be seen on average five times; beyond that the entries are NaN. Past
    that point the block histograms are so sparse that every estimator here (Grassberger included) is biased low
    by an amount comparable to the structure one is looking for, so a falling curve would mean nothing.

    Output: float array of length max_order + 1. Block codes are built incrementally so all orders cost one pass
    each (10**6 binary symbols to order 8 take well under a second).
    """
    s = _as_symbols(seq, k)
    out = np.full(max_order + 1, np.nan)
    top = max_trustworthy_order(s.size, k, max_order)
    if s.size == 0 or (k ** 0 * 5 > s.size):
        return out
    h_prev = 0.0
    for n, codes, span in _iter_block_codes(s, k, top + 1):
        h = entropy(_count_codes(codes, span), method)
        out[n - 1] = h - h_prev
        h_prev = h
    return out


# ----------------------------------------------------------------------------------------------------------------
# lagged mutual information
# ----------------------------------------------------------------------------------------------------------------

def _mi_from_pairs(a: np.ndarray, b: np.ndarray, k: int) -> tuple[float, int, int, int, int]:
    """Plug-in mutual information in bits of paired symbols, plus (m_a, m_b, m_ab, N) for bias corrections."""
    n = int(a.size)
    if n == 0:
        return float("nan"), 0, 0, 0, 0
    joint = a * k + b
    if k * k <= max(_DENSE_LIMIT, 4 * n):
        cab = np.bincount(joint, minlength=k * k).reshape(k, k)
        ca, cb = cab.sum(1), cab.sum(0)
        nz = cab > 0
        num = cab[nz]
        mi = float((num * np.log(num * n / np.outer(ca, cb)[nz])).sum()) / n / _LN2
        return mi, int((ca > 0).sum()), int((cb > 0).sum()), int(nz.sum()), n
    codes, cab = np.unique(joint, return_counts=True)
    ua, ia = np.unique(a, return_inverse=True)
    ca = np.bincount(ia)
    ub, ib = np.unique(b, return_inverse=True)
    cb = np.bincount(ib)
    pa = ca[np.searchsorted(ua, codes // k)]
    pb = cb[np.searchsorted(ub, codes % k)]
    mi = float((cab * np.log(cab * n / (pa * pb.astype(float)))).sum()) / n / _LN2
    return mi, int(ua.size), int(ub.size), int(codes.size), n


def _sorted_time(time: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(time, dtype=float).ravel()
    if t.size != n:
        raise ValueError("time must have the same length as seq")
    if not np.all(np.isfinite(t)):
        raise ValueError("time must be finite")
    order = np.argsort(t, kind="stable")
    return t[order], order


def _time_pairs(t: np.ndarray, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray]:
    """Index pairs (i, j), i < j, of a sorted time array with lo <= t[j] - t[i] < hi."""
    n = t.size
    idx = np.arange(n)
    start = np.maximum(np.searchsorted(t, t + lo, side="left"), idx + 1)
    stop = np.searchsorted(t, t + hi, side="left")
    cnt = np.clip(stop - start, 0, None)
    total = int(cnt.sum())
    if total == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    i = np.repeat(idx, cnt)
    offset = np.arange(total) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    j = np.repeat(start, cnt) + offset
    return i, j


def lag_pair_index(n: int, lags: np.ndarray | Sequence[int] | None = None, *, time: np.ndarray | None = None,
                   lag_edges: np.ndarray | Sequence[float] | None = None) -> tuple[list, np.ndarray]:
    """The pairing used by the lagged statistics.

    Index mode (time is None): one entry per lag, a pair of slices; time-separation mode: one entry per lag bin
    [lag_edges[b], lag_edges[b+1]), a pair of index arrays into the time-sorted sequence. Returns (pairings, order),
    where `order` sorts seq by time (identity in index mode).
    """
    if time is None and lag_edges is None:
        if lags is None:
            raise ValueError("give lags, or time together with lag_edges")
        lags = np.asarray(lags, dtype=np.int64).ravel()
        if lags.size and lags.min() < 1:
            raise ValueError("index lags must be positive integers")
        return [(slice(0, max(n - int(d), 0)), slice(int(d), n)) for d in lags], np.arange(n)
    if time is None or lag_edges is None:
        raise ValueError("time and lag_edges must be given together")
    t, order = _sorted_time(time, n)
    edges = np.asarray(lag_edges, dtype=float).ravel()
    if edges.size < 2 or np.any(np.diff(edges) <= 0) or edges[0] < 0:
        raise ValueError("lag_edges must be a non-negative, strictly increasing array of at least two edges")
    return [_time_pairs(t, lo, hi) for lo, hi in zip(edges[:-1], edges[1:])], order


def _mi_curve(s: np.ndarray, k: int, pairings: list) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Plug-in MI per pairing, the Miller-Madow correction term per pairing (bits), and the number of pairs."""
    mi = np.full(len(pairings), np.nan)
    mm = np.zeros(len(pairings))
    npairs = np.zeros(len(pairings), dtype=np.int64)
    for b, (ia, ib) in enumerate(pairings):
        value, ma, mb, mab, n = _mi_from_pairs(s[ia], s[ib], k)
        mi[b], npairs[b] = value, n
        if n:
            mm[b] = ((ma - 1) + (mb - 1) - (mab - 1)) / (2.0 * n * _LN2)
    return mi, mm, npairs


def shuffled_mi_null(s: np.ndarray, k: int, pairings: list, n_shuffle: int, rng: np.random.Generator) -> np.ndarray:
    """Plug-in MI for each of n_shuffle random permutations of the symbols: array (n_shuffle, n_lags)."""
    out = np.empty((n_shuffle, len(pairings)))
    for r in range(n_shuffle):
        out[r] = _mi_curve(rng.permutation(s), k, pairings)[0]
    return out


def lagged_mutual_information(seq: np.ndarray | Sequence[int], k: int, lags: np.ndarray | Sequence[int] | None, *,
                              time: np.ndarray | None = None, lag_edges: np.ndarray | Sequence[float] | None = None,
                              correction: str = "shuffle", rng: np.random.Generator | None = None,
                              n_shuffle: int = 200) -> np.ndarray:
    """Mutual information in bits between symbols separated by each lag.

    Index mode (default): MI between seq[i] and seq[i + lag] for each positive integer lag in `lags`.
    Time mode (`time` and `lag_edges` given; `lags` is ignored and may be None): the sequence is sorted by time and,
    for each bin [lag_edges[b], lag_edges[b + 1]) in the units of `time`, MI is computed over all pairs (i, j) with
    i before j and t_j - t_i in the bin. A symbol then appears in several pairs, so pairs are not independent; the
    shuffle correction accounts for that, the Miller-Madow correction (which counts pairs as samples) does not.

    correction: "shuffle" (default) subtracts the mean plug-in MI over `n_shuffle` random permutations of the
    symbols over the same positions or times (needs `rng`); "miller_madow" applies Miller-Madow to the three
    entropies, i.e. subtracts ((m_ab - 1) - (m_a - 1) - (m_b - 1)) / (2 N ln 2); "none" returns the plug-in value.
    Output: float array, one value per lag or bin; NaN where no pair exists. Corrected values can be negative.
    """
    s = _as_symbols(seq, k)
    pairings, order = lag_pair_index(s.size, lags, time=time, lag_edges=lag_edges)
    s = s[order]
    mi, mm, _ = _mi_curve(s, k, pairings)
    if correction == "none":
        return mi
    if correction == "miller_madow":
        return mi + mm
    if correction == "shuffle":
        if rng is None:
            raise ValueError("correction='shuffle' needs rng")
        if n_shuffle < 1:
            raise ValueError("n_shuffle must be positive")
        return mi - shuffled_mi_null(s, k, pairings, n_shuffle, rng).mean(axis=0)
    raise ValueError("correction must be 'shuffle', 'miller_madow' or 'none'")


# ----------------------------------------------------------------------------------------------------------------
# compression
# ----------------------------------------------------------------------------------------------------------------

def _pack(s: np.ndarray, k: int) -> bytes:
    if k <= 2:
        return np.packbits(s.astype(np.uint8)).tobytes()
    if k <= 256:
        return s.astype(np.uint8).tobytes()
    if k <= 1 << 16:
        return s.astype("<u2").tobytes()
    return s.astype("<u4").tobytes()


def compression_bits(seq: np.ndarray | Sequence[int], k: int, method: str = "lzma") -> float:
    """Compressed size in bits of the sequence at the compressor's highest level.

    Packing: eight symbols per byte when k <= 2, one byte per symbol when k <= 256, two bytes (little-endian) when
    k <= 65536, four bytes otherwise. method: "zlib" (level 9), "bz2" (level 9) or "lzma" (preset 9 | EXTREME, xz
    container). The container headers are included (a few hundred bits at most), so compare sequences of similar
    length or compare against surrogates rather than reading the absolute value.
    """
    s = _as_symbols(seq, k)
    data = _pack(s, k)
    if method == "zlib":
        out = zlib.compress(data, 9)
    elif method == "bz2":
        out = bz2.compress(data, 9)
    elif method == "lzma":
        out = lzma.compress(data, preset=9 | lzma.PRESET_EXTREME)
    else:
        raise ValueError("method must be 'zlib', 'bz2' or 'lzma'")
    return 8.0 * len(out)


def compression_ratio(seq: np.ndarray | Sequence[int], k: int, method: str = "lzma") -> float:
    """compression_bits divided by the raw size len(seq) * log2(k); NaN when that is zero (k = 1 or empty)."""
    s = _as_symbols(seq, k)
    raw = s.size * math.log2(k) if k > 1 else 0.0
    if raw == 0:
        return float("nan")
    return compression_bits(s, k, method) / raw
