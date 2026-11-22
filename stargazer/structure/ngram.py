# Adapted from SeaQuest (same author, commit 064d667), seaquest/structure/sequences.py. Changes: `order` now means
# the context length (order 0 = unigram; SeaQuest's NGramModel(order=n) used n - 1 context symbols); the alphabet
# size k is given, so the floor is uniform 1/k instead of 1/(V + 1); counting and scoring are vectorised with
# integer-encoded contexts held in sorted sparse tables (exposed as dictionaries) instead of Python Counters; the
# gain-in-bits bootstrap splits train/test first and resamples within each side (SeaQuest resampled before
# splitting, leaking test data into training); its null is >= 200 within-sequence shuffles with an empirical
# p-value instead of one shuffle; a single long sequence is split into contiguous blocks instead of being scored
# in-sample; transition-matrix z-scores are replaced by empirical two-sided p-values with Benjamini-Hochberg
# adjustment over the k^2 cells.
"""Interpolated Witten-Bell n-gram models, held-out bits per symbol, gain in bits, and transition matrices."""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .estimators import _as_symbols

__all__ = ["NGramModel", "benjamini_hochberg", "gain_in_bits", "held_out_bits_per_symbol", "transition_matrix"]

SequenceInput = np.ndarray | Sequence[np.ndarray] | Sequence[Sequence[int]]


def _as_sequences(sequences: SequenceInput, k: int) -> list[np.ndarray]:
    """A list of 1-D int64 arrays from one array (one sequence), a 2-D array (rows) or a list of sequences."""
    if isinstance(sequences, np.ndarray):
        if sequences.ndim == 1:
            return [_as_symbols(sequences, k)]
        return [_as_symbols(row, k) for row in sequences]
    seqs = list(sequences)
    if seqs and np.isscalar(seqs[0]):
        return [_as_symbols(np.asarray(seqs), k)]
    return [_as_symbols(s, k) for s in seqs]


def _block_codes(s: np.ndarray, k: int, length: int) -> np.ndarray:
    """Base-k codes of the overlapping blocks of the given length (length 0 gives zeros, one per position + 1)."""
    n = s.size
    if length == 0:
        return np.zeros(n + 1, dtype=np.int64)
    if n < length:
        return np.zeros(0, dtype=np.int64)
    codes = s[:n - length + 1].copy()
    for j in range(1, length):
        codes = codes * k + s[j:n - length + 1 + j]
    return codes


def _concat(seqs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate sequences; also return each position's offset within its own sequence."""
    if not seqs:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    if len(seqs) == 1:
        return seqs[0], np.arange(seqs[0].size)
    s = np.concatenate(seqs)
    lengths = np.array([x.size for x in seqs])
    starts = np.repeat(np.cumsum(lengths) - lengths, lengths)
    return s, np.arange(s.size) - starts


def _unique_counts(codes: np.ndarray, span: int) -> tuple[np.ndarray, np.ndarray]:
    """Sorted distinct codes and their counts (bincount when the code space is small, else a sort)."""
    if codes.size and span <= max(1 << 16, 4 * codes.size):
        c = np.bincount(codes, minlength=0)
        keys = np.nonzero(c)[0]
        return keys.astype(np.int64), c[keys]
    return np.unique(codes, return_counts=True)


def _lookup(keys: np.ndarray, values: np.ndarray, query: np.ndarray) -> np.ndarray:
    if keys.size == 0:
        return np.zeros(query.shape, dtype=values.dtype)
    idx = np.searchsorted(keys, query)
    idx = np.minimum(idx, keys.size - 1)
    return np.where(keys[idx] == query, values[idx], 0)


class NGramModel:
    """Interpolated Witten-Bell n-gram model over the alphabet {0, ..., k - 1}.

    order is the context length: P_order(x | c) uses the `order` preceding symbols. The recursion is
        P_-1(x) = 1 / k,
        P_m(x | c_m) = lam(c_m) * n(c_m, x) / n(c_m) + (1 - lam(c_m)) * P_{m-1}(x | c_{m-1}),
        lam(c) = n(c) / (n(c) + t(c)),
    with n(c) the number of times context c (the last m symbols) was followed by any symbol, t(c) the number of
    distinct symbols that followed it, and P_m = P_{m-1} for unseen contexts. Contexts never cross the boundaries of
    the training sequences, and there is no start padding: the first m symbols of a sequence are scored with the
    shorter contexts available.

    Storage is sparse: for each level m, the sorted integer codes of the observed (context, symbol) pairs and of the
    observed contexts with their counts. `joint_counts(m)` and `context_counts(m)` return them as dictionaries
    keyed by the integer-encoded context (base-k value of the context symbols, oldest first). k**(order + 1) must
    fit in 63 bits.
    """

    def __init__(self, order: int, k: int):
        if order < 0:
            raise ValueError("order must be non-negative")
        if k < 1:
            raise ValueError("k must be at least 1")
        if k ** (order + 1) >= 1 << 62:
            raise ValueError("k**(order + 1) must fit in a 64-bit integer")
        self.order, self.k = int(order), int(k)
        self._joint: list[tuple[np.ndarray, np.ndarray]] = []
        self._ctx: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        self.n_train = 0

    def fit(self, sequences: SequenceInput) -> "NGramModel":
        """Count (context, symbol) pairs at every level 0..order. Replaces any previous fit. Returns self."""
        return self._fit_concat(*_concat(_as_sequences(sequences, self.k)))

    def _fit_concat(self, s: np.ndarray, offset: np.ndarray) -> "NGramModel":
        """fit() on validated, concatenated sequences; offset[i] is position i's index within its own sequence."""
        self.n_train = int(s.size)
        self._joint, self._ctx = [], []
        for m in range(self.order + 1):
            codes = _block_codes(s, self.k, m + 1)
            if m and codes.size:
                codes = codes[offset[m:] >= m]          # drop blocks that cross a sequence boundary
            keys, counts = _unique_counts(codes, self.k ** (m + 1))
            ctx = keys // self.k
            ctx_keys, start, n_types = np.unique(ctx, return_index=True, return_counts=True)
            totals = np.add.reduceat(counts, start) if counts.size else np.zeros(0, dtype=np.int64)
            self._joint.append((keys, counts))
            self._ctx.append((ctx_keys, totals, n_types))
        return self

    def joint_counts(self, m: int) -> dict[int, int]:
        """{context_code * k + symbol: count} at level m (context length m)."""
        keys, counts = self._joint[m]
        return dict(zip(keys.tolist(), counts.tolist()))

    def context_counts(self, m: int) -> dict[int, tuple[int, int]]:
        """{context_code: (n(c), t(c))} at level m: tokens following the context and distinct followers."""
        keys, totals, types = self._ctx[m]
        return dict(zip(keys.tolist(), zip(totals.tolist(), types.tolist())))

    def level_log2_probs(self, sequence: SequenceInput) -> np.ndarray:
        """log2 P_m(x_i | context) for every level m = 0..order and position i: array (order + 1, n).

        Row m is what an order-m model fitted to the same data would give, so one fit scores all lower orders.
        `sequence` may also be a list of sequences; they are scored in one pass (concatenated, n = total length)
        with contexts that never cross from one sequence into the next.
        """
        if not self._joint:
            raise RuntimeError("fit the model first")
        return self._score_concat(*_concat(_as_sequences(sequence, self.k)))

    def _score_concat(self, s: np.ndarray, offset: np.ndarray) -> np.ndarray:
        """level_log2_probs() on validated, concatenated sequences (see _fit_concat for offset)."""
        n = s.size
        p = np.full(n, 1.0 / self.k)
        out = np.empty((self.order + 1, n))
        prev_block = np.zeros(n + 1, dtype=np.int64)
        for m in range(self.order + 1):
            if n > m:
                block = _block_codes(s, self.k, m + 1)
                ctx = prev_block[:n - m]
                ctx_keys, totals, types = self._ctx[m]
                n_c = _lookup(ctx_keys, totals, ctx).astype(float)
                t_c = _lookup(ctx_keys, types, ctx).astype(float)
                n_cx = _lookup(self._joint[m][0], self._joint[m][1], block).astype(float)
                seen = n_c > 0
                if m:
                    seen &= offset[m:] >= m             # the context would reach into the previous sequence
                tail = p[m:]
                safe_n = np.where(seen, n_c, 1.0)
                lam = n_c / (safe_n + t_c)
                tail[:] = np.where(seen, lam * n_cx / safe_n + (1.0 - lam) * tail, tail)
                prev_block = block
            out[m] = np.log2(p)
        return out

    def log2_prob(self, sequence: np.ndarray) -> float:
        """Total log2 probability (bits, <= 0) of `sequence` under the full order-`order` model."""
        return float(self.level_log2_probs(sequence)[-1].sum())


def held_out_bits_per_symbol(train: SequenceInput, test: SequenceInput, k: int, order: int) -> float:
    """Bits per symbol, -sum log2 P / (number of test symbols), of the test data under a model fitted to train."""
    model = NGramModel(order, k).fit(train)
    tests = _as_sequences(test, k)
    n = sum(t.size for t in tests)
    if n == 0:
        return float("nan")
    return float(-sum(model.log2_prob(t) for t in tests) / n)


def _unit_index(lengths: np.ndarray, units: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Gather positions (into the concatenation of all units) and within-unit offsets for the chosen units."""
    starts = np.cumsum(lengths) - lengths
    sel_len = lengths[units]
    total = int(sel_len.sum())
    first = np.repeat(np.cumsum(sel_len) - sel_len, sel_len)
    offset = np.arange(total) - first
    return np.repeat(starts[units], sel_len) + offset, offset


def gain_in_bits(sequences: SequenceInput, k: int, order: int, n_boot: int, rng: np.random.Generator, *,
                 test_frac: float = 0.3, n_shuffle: int = 200, n_blocks: int = 10) -> dict:
    """Held-out bits per symbol of an order-0 model minus an order-`order` model (Witten-Bell, see NGramModel).

    Units: with two or more sequences, each sequence is a unit; a single sequence is cut into `n_blocks`
    contiguous blocks (np.array_split), which become the units (contexts do not cross block joins). A random split
    puts round(test_frac * units) units (at least one, and at least one left for training) in the test set.
    Both models come from one Witten-Bell fit (the order-0 model is level 0 of the order-`order` model).

    Bootstrap (n_boot draws; 0 skips it): with the split fixed, training units are resampled with replacement from
    the training side and test units from the test side, so no test data can reach training (SeaQuest resampled
    first and split afterwards). ci_low / ci_high are the 2.5 and 97.5 percentiles (NaN if n_boot = 0); the
    interval describes resampling variability given this split.

    Null (n_shuffle draws, use at least 200): the symbols are permuted within each sequence (a single sequence is
    permuted as a whole and then re-cut into the same blocks), the same split is applied and the gain recomputed;
    p_value = (1 + #{null >= gain}) / (n_shuffle + 1). Under an i.i.d. null the data and its shuffles are
    exchangeable, so this p-value is exact up to its resolution 1 / (n_shuffle + 1).

    Output dict: gain_bits, bits_order0, bits_order (held-out bits per symbol), ci_low, ci_high, boot_gain_bits,
    null_gain_bits, null_mean, p_value, order, n_units, split ("sequences" or "blocks"), n_test_units.
    """
    seqs = _as_sequences(sequences, k)
    if order < 1:
        raise ValueError("order must be at least 1 (order 0 is the reference model)")
    blocked = len(seqs) < 2
    if blocked:
        s_all = seqs[0] if seqs else np.zeros(0, dtype=np.int64)
        if s_all.size < 2 * n_blocks:
            raise ValueError("a single sequence needs at least 2 * n_blocks symbols")
        lengths = np.array([b.size for b in np.array_split(s_all, n_blocks)])
        seq_id = np.zeros(s_all.size, dtype=np.int64)
    else:
        s_all = np.concatenate(seqs)
        lengths = np.array([x.size for x in seqs])
        seq_id = np.repeat(np.arange(len(seqs)), lengths)
    n_units = lengths.size
    n_test = int(min(max(1, round(test_frac * n_units)), n_units - 1))
    perm = rng.permutation(n_units)
    test_idx, train_idx = perm[:n_test], perm[n_test:]
    tr_pos, tr_off = _unit_index(lengths, train_idx)
    te_pos, te_off = _unit_index(lengths, test_idx)

    def bits_of(s: np.ndarray, trp: np.ndarray, tro: np.ndarray, tep: np.ndarray, teo: np.ndarray) -> np.ndarray:
        model = NGramModel(order, k)._fit_concat(s[trp], tro)
        lp = model._score_concat(s[tep], teo)
        return -lp.sum(axis=1) / max(lp.shape[1], 1)

    bits = bits_of(s_all, tr_pos, tr_off, te_pos, te_off)
    gain = float(bits[0] - bits[order])
    boots = np.empty(n_boot)
    for r in range(n_boot):
        trp, tro = _unit_index(lengths, rng.choice(train_idx, size=train_idx.size, replace=True))
        tep, teo = _unit_index(lengths, rng.choice(test_idx, size=test_idx.size, replace=True))
        b = bits_of(s_all, trp, tro, tep, teo)
        boots[r] = b[0] - b[order]
    null = np.empty(n_shuffle)
    for r in range(n_shuffle):
        if blocked:
            shuffled = rng.permutation(s_all)
        else:
            shuffled = s_all[np.lexsort((rng.random(s_all.size), seq_id))]
        b = bits_of(shuffled, tr_pos, tr_off, te_pos, te_off)
        null[r] = b[0] - b[order]
    p = (1.0 + float(np.count_nonzero(null >= gain))) / (n_shuffle + 1.0)
    ci = np.percentile(boots, [2.5, 97.5]) if n_boot else np.array([np.nan, np.nan])
    return {"gain_bits": gain, "bits_order0": float(bits[0]), "bits_order": float(bits[order]),
            "ci_low": float(ci[0]), "ci_high": float(ci[1]), "boot_gain_bits": boots, "null_gain_bits": null,
            "null_mean": float(null.mean()) if n_shuffle else float("nan"), "p_value": p, "order": int(order),
            "n_units": int(n_units), "split": "blocks" if blocked else "sequences", "n_test_units": n_test}


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values controlling the false discovery rate), same shape as p."""
    p = np.asarray(p, dtype=float)
    flat = p.ravel()
    m = flat.size
    if m == 0:
        return p.copy()
    order = np.argsort(flat)
    scaled = flat[order] * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adjusted, 1.0)
    return out.reshape(p.shape)


def _transition_counts(seqs: list[np.ndarray], k: int) -> np.ndarray:
    codes = [s[:-1] * k + s[1:] for s in seqs if s.size > 1]
    if not codes:
        return np.zeros((k, k), dtype=np.int64)
    return np.bincount(np.concatenate(codes), minlength=k * k).reshape(k, k)


def transition_matrix(sequences: SequenceInput, k: int, n_shuffle: int, rng: np.random.Generator) -> dict:
    """First-order transition counts and probabilities, with per-cell tests against within-sequence shuffles.

    For each cell (a, b) the observed count is compared with its distribution over n_shuffle shuffles (each
    sequence permuted independently, which keeps symbol frequencies). The two-sided empirical p-value is
    min(1, 2 * min(p_upper, p_lower)) with p_upper = (1 + #{null >= observed}) / (n + 1) and p_lower likewise;
    cells whose count cannot vary under shuffling get p = 1. p_adjusted applies Benjamini-Hochberg over all k^2
    cells. The table is dense (k x k), so this is meant for alphabets up to a few hundred symbols.

    Output dict: counts (k, k) int, probabilities (row-normalised; zero rows stay zero), expected (mean shuffled
    count), p_values, p_adjusted, n_shuffle.
    """
    seqs = _as_sequences(sequences, k)
    counts = _transition_counts(seqs, k)
    null = np.empty((n_shuffle, k, k), dtype=np.int64)
    for r in range(n_shuffle):
        null[r] = _transition_counts([rng.permutation(s) for s in seqs], k)
    upper = (1.0 + (null >= counts).sum(axis=0)) / (n_shuffle + 1.0)
    lower = (1.0 + (null <= counts).sum(axis=0)) / (n_shuffle + 1.0)
    p = np.minimum(1.0, 2.0 * np.minimum(upper, lower))
    rows = counts.sum(axis=1, keepdims=True)
    prob = np.divide(counts, rows, out=np.zeros((k, k)), where=rows > 0)
    return {"counts": counts, "probabilities": prob, "expected": null.mean(axis=0), "p_values": p,
            "p_adjusted": benjamini_hochberg(p), "n_shuffle": int(n_shuffle)}
