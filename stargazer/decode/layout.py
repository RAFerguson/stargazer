"""Geometric layout of an unknown bit stream: is it a picture sent row by row, and of what width; is it a volume;
where does a row start; is there a header in front.

The measure is code length. Every model here is a *context model*: each bit is coded with the Krichevsky-Trofimov
(KT) estimator of its context, a sequential, decodable adaptive code. Because the KT estimator is exchangeable,
the total code length of a context model depends only on the (context, bit) counts, so each model costs one
`np.bincount` and no Python loop over bits. A "gain" is therefore a real number of bits saved by a real code over
another real code, and the no-hypercompression inequality turns it into a probability bound with no simulation.

* The null (baseline) code: context = the previous m bits, m = 0..8, the best m chosen and paid for with log2(9)
  bits. It represents "a short-memory binary stream with no geometry".
* The 2-D code at width W: context = left (i-1), up-left (i-W-1), up (i-W) and up-right (i-W+1), optionally with
  i-2, or i-2 and i-2W, as well, or the small context (left, up, i-2); the best of these four is used and paid for
  with 2 bits.
* The 3-D code at width W and plane height H: each of those contexts plus the bit one plane back (i-W*H).

What this cannot see: layouts that are not row-major rasters with a fixed row length (variable-length lines,
interlacing, column-major order is seen as a raster of the other width, spiral or Hilbert scans), pictures whose
neighbouring pixels are not correlated (dithered or encrypted images), and any structure that a short-memory 1-D
code already captures as well (a raster that is mostly one colour is compressed by the baseline too, and only the
*extra* saving counts). A byte, character or frame structure also makes a stream fold well at the frame length
or a multiple: a significant gain says "the stream has period-W structure", not "the stream is a picture".
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.special import gammaln

from stargazer.core.rng import rng_for
from stargazer.core.types import BitStream

__all__ = [
    "FoldResult", "VolumeResult", "baseline_bits", "factorisations", "fold_scan", "folded_bits", "kt_bits",
    "render", "to_text", "volume_scan",
]

_LN2 = math.log(2.0)
_LOG2_PI = math.log2(math.pi)
_MAX_ORDER = 8


# ---------------------------------------------------------------------------------------------------------------
# Code lengths.
# ---------------------------------------------------------------------------------------------------------------


def kt_bits(n0: Any, n1: Any) -> np.ndarray | float:
    """Code length in bits of a binary sequence with n0 zeros and n1 ones under the Krichevsky-Trofimov estimator.

    -log2( G(n0 + 1/2) G(n1 + 1/2) / (pi G(n0 + n1 + 1)) ), with G the gamma function. This is the exact length of
    the sequential code that predicts each bit with P(1) = (ones so far + 1/2) / (bits so far + 1); it does not
    depend on the order of the bits. Vectorised: arrays broadcast, scalars give a float. Zero counts give 0 bits.
    """
    a = np.asarray(n0, dtype=float)
    b = np.asarray(n1, dtype=float)
    out = (gammaln(a + b + 1.0) - gammaln(a + 0.5) - gammaln(b + 0.5)) / _LN2 + _LOG2_PI
    if out.ndim == 0:
        return float(out)
    return out


def _context_cost(counts: np.ndarray) -> float:
    """Total KT code length of a context model from its counts, shape (n_contexts, 2)."""
    c = counts.reshape(-1, 2)
    used = (c[:, 0] + c[:, 1]) > 0
    return float(np.sum(kt_bits(c[used, 0], c[used, 1])))


def _prepare(bits: np.ndarray, known: np.ndarray | None) -> tuple[np.ndarray, np.ndarray | None]:
    """Bits as uint8 with erased bits set to 0 (so erased context bits read as 0), and the target mask."""
    b = np.asarray(bits).astype(np.uint8).ravel()
    if known is None:
        return b, None
    k = np.asarray(known, dtype=bool).ravel()
    if k.shape != b.shape:
        raise ValueError("known must match bits in length")
    if k.all():
        return b, None
    return np.where(k, b, 0).astype(np.uint8), k


def _padded(b: np.ndarray, pad: int) -> np.ndarray:
    """The stream with `pad` zeros in front, so that b[i - o] is padded[pad + i - o] ("missing = 0")."""
    return np.concatenate([np.zeros(pad, dtype=np.uint8), b])


def _shift(padded: np.ndarray, pad: int, n: int, offset: int) -> np.ndarray:
    """The bits at i - offset for i = 0..n-1 (zero before the start of the stream)."""
    return padded[pad - offset:pad - offset + n]


def _counts(ctx: np.ndarray, b: np.ndarray, mask: np.ndarray | None, n_ctx: int) -> np.ndarray:
    """(context, bit) counts over the coded (known) bits, shape (n_ctx, 2)."""
    code = ctx.astype(np.int64) * 2 + b
    if mask is not None:
        code = code[mask]
    return np.bincount(code, minlength=2 * n_ctx).reshape(n_ctx, 2)


def _baseline_counts(b: np.ndarray, mask: np.ndarray | None, max_order: int) -> list[np.ndarray]:
    n = b.size
    pad = max_order + 1
    p = _padded(b, pad)
    out = []
    ctx = np.zeros(n, dtype=np.int64)
    for m in range(max_order + 1):
        if m:
            ctx = ctx * 2 + _shift(p, pad, n, m)
        out.append(_counts(ctx, b, mask, 1 << m))
    return out


def baseline_bits(bits: np.ndarray, known: np.ndarray | None = None, max_order: int = _MAX_ORDER
                  ) -> tuple[float, int]:
    """Code length of the null model, and the order it chose: (bits, m).

    The null code is a KT context model whose context is the previous m bits (bits before the start read as 0),
    for the best m in 0..max_order, plus log2(max_order + 1) bits to say which m. It is the code of "a binary
    source with memory of at most max_order bits and no geometry". Erased bits (known False) are not coded and read
    as 0 when they are context. It cannot capture memory longer than max_order bits, so a stream with long but
    one-dimensional memory (text, for example) may let a 2-D model at a small width win for that reason alone.
    """
    b, mask = _prepare(bits, known)
    costs = [_context_cost(c) for c in _baseline_counts(b, mask, max_order)]
    m = int(np.argmin(costs))
    return costs[m] + math.log2(max_order + 1), m


# The 2-D contexts tried at every width, as the bits of the six-bit context they keep: 0 left (i-1),
# 1 up-left (i-W-1), 2 up (i-W), 3 up-right (i-W+1), 4 left-left (i-2), 5 up-up (i-2W).
CONTEXTS: dict[str, tuple[int, ...]] = {
    "left+up3": (0, 1, 2, 3),
    "left+up3+left2": (0, 1, 2, 3, 4),
    "left+up3+left2+up2": (0, 1, 2, 3, 4, 5),
    "left+up+left2": (0, 2, 4),
}
MODELS: tuple[str, ...] = tuple(CONTEXTS)
_CHOICE_BITS = math.log2(len(MODELS))


def _marginal(c: np.ndarray, keep: tuple[int, ...], batch: int = 0) -> np.ndarray:
    """Counts of a sub-context from full (context, bit) counts shaped (2,) * (k + 1), context bit j on axis k-1-j
    (after `batch` leading axes, which are kept)."""
    k = c.ndim - 1 - batch
    drop = tuple(batch + k - 1 - j for j in range(k) if j not in keep)
    return c.sum(axis=drop) if drop else c


def _batch_cost(c: np.ndarray) -> np.ndarray:
    """KT cost of each of a batch of context models, counts shaped (batch, ..., 2)."""
    return np.asarray(kt_bits(c[..., 0], c[..., 1])).reshape(c.shape[0], -1).sum(axis=1)


class _Folder:
    """A stream prepared for fast context counts at many widths (padded with zeros in front, erased bits 0)."""

    def __init__(self, bits: np.ndarray, known: np.ndarray | None, max_width: int) -> None:
        self.b, self.mask = _prepare(bits, known)
        self.n = self.b.size
        self.pad = 2 * max_width + 2
        self.p = _padded(self.b, self.pad)
        self.target = np.arange(self.n) if self.mask is None else np.flatnonzero(self.mask)

    def sh(self, offset: int) -> np.ndarray:
        return _shift(self.p, self.pad, self.n, offset)

    def ctx6(self, width: int) -> np.ndarray:
        """Six-bit context: bit 0 left, 1 up-left, 2 up, 3 up-right, 4 left-left (i-2), 5 up-up (i-2W)."""
        w = width
        return (self.sh(1) | (self.sh(w + 1) << 1) | (self.sh(w) << 2) | (self.sh(w - 1) << 3)
                | (self.sh(2) << 4) | (self.sh(2 * w) << 5))

    def costs(self, width: int) -> np.ndarray:
        """Code length at `width` under each model of MODELS (in that order), without the choice bits."""
        return self.costs_many(np.array([width]))[0]

    def folded(self, width: int) -> float:
        return float(self.costs(width).min()) + _CHOICE_BITS

    def costs_many(self, widths: np.ndarray, back: np.ndarray | None = None, block_elements: int = 1_000_000
                   ) -> np.ndarray:
        """Code lengths at many widths, shape (len(widths), len(MODELS)). With `back` (one lag per width), every
        model also conditions on the bit at i - back (the volume code). Short streams are counted many widths at
        a time (blocks of about `block_elements` entries); long ones one width at a time from contiguous slices."""
        widths = np.asarray(widths, dtype=np.int64)
        backs = None if back is None else np.asarray(back, dtype=np.int64)
        extra = 0 if backs is None else 1
        n_loc = 128 << extra
        counts = np.empty((widths.size, n_loc), dtype=np.int64)
        step = max(1, block_elements // max(self.n, 1))
        if self.n <= 4000:                           # gathering beats slicing only while Python overhead dominates
            idx = self.pad + self.target
            p = self.p
            near = (p[idx - 1].astype(np.int32) | (p[idx - 2].astype(np.int32) << 4)) * 2 + self.b[self.target]
            for a in range(0, widths.size, step):
                w = widths[a:a + step, None]
                far = (p[idx - w - 1] | (p[idx - w] << 1) | (p[idx - w + 1] << 2) | (p[idx - 2 * w] << 4))
                offset = (n_loc * np.arange(w.shape[0], dtype=np.int32))[:, None]
                code = near[None] + (far.astype(np.int32) << 2) + offset
                if backs is not None:
                    code += p[idx - backs[a:a + step, None]].astype(np.int32) << 7
                counts[a:a + step] = np.bincount(code.ravel(), minlength=n_loc * w.shape[0]).reshape(-1, n_loc)
        else:
            near = (self.sh(1) | (self.sh(2) << 4)).astype(np.int32) * 2 + self.b
            for j, w in enumerate(widths.tolist()):
                far = self.sh(w + 1) | (self.sh(w) << 1) | (self.sh(w - 1) << 2) | (self.sh(2 * w) << 4)
                code = near + (far.astype(np.int32) << 2)
                if backs is not None:
                    code += self.sh(int(backs[j])).astype(np.int32) << 7
                if self.mask is not None:
                    code = code[self.mask]
                counts[j] = np.bincount(code, minlength=n_loc)
        c = counts.reshape((widths.size,) + (2,) * (7 + extra))
        out = np.empty((widths.size, len(MODELS)))
        for k, keep in enumerate(CONTEXTS.values()):
            out[:, k] = _batch_cost(_marginal(c, keep + ((6,) if extra else ()), batch=1))
        return out


def folded_bits(bits: np.ndarray, width: int, known: np.ndarray | None = None) -> float:
    """Code length of the 2-D context model at row length `width`, in bits.

    The best of four KT context models (see CONTEXTS), plus log2(4) = 2 bits for the choice: (left, up-left, up,
    up-right), that is bits i-1, i-W-1, i-W, i-W+1; that context with i-2 added; with i-2 and i-2W added; and the
    small context (left, up, i-2), which wins on noisy or short pictures because it has fewer counts to learn.
    Bits before the start of the stream read as 0, and so do erased bits; erased bits are not coded.
    """
    if width < 2:
        raise ValueError("width must be at least 2")
    return _Folder(bits, known, width).folded(width)


# ---------------------------------------------------------------------------------------------------------------
# Factorisation.
# ---------------------------------------------------------------------------------------------------------------


def _prime_factors(n: int) -> list[int]:
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def factorisations(n: int) -> dict:
    """The prime factors of n and every way to write it as width x height with both at least 2.

    Returns {"n": n, "prime_factors": [p, ...] (with repetition), "pairs": [(w, h), ...] in increasing w,
    "semiprime": bool, "prime": bool}.
    """
    n = int(n)
    if n < 1:
        raise ValueError("n must be positive")
    pf = _prime_factors(n) if n > 1 else []
    pairs = [(w, n // w) for w in range(2, n // 2 + 1) if n % w == 0]
    return {"n": n, "prime_factors": pf, "pairs": pairs, "semiprime": len(pf) == 2, "prime": len(pf) == 1}


def _factor_note(n: int, width: int) -> str:
    f = factorisations(n)
    pf = f["prime_factors"]
    if n < 2:
        return f"{n}: no factorisation"
    if f["prime"]:
        note = f"{n} is prime: no exact rectangle"
    elif f["semiprime"]:
        a, b = pf
        note = f"{n} = {a} x {b}: both prime" if a != b else f"{n} = {a} x {a}: a prime squared"
    else:
        powers = []
        for p in sorted(set(pf)):
            k = pf.count(p)
            powers.append(f"{p}^{k}" if k > 1 else str(p))
        note = f"{n} = {' x '.join(powers)}: {len(f['pairs'])} exact rectangles"
    if n % width == 0:
        note += f"; width {width} divides it ({n // width} rows)"
    else:
        note += f"; width {width} leaves a remainder of {n % width} bits"
    return note


# ---------------------------------------------------------------------------------------------------------------
# Seam and header.
# ---------------------------------------------------------------------------------------------------------------


def _seam(fold: _Folder, width: int) -> tuple[int, float]:
    """Column where a row really starts, and the evidence for it (see FoldResult.seam)."""
    w, n = width, fold.n
    if n < 8 * w:
        return 0, 0.0
    x = fold.sh(1).astype(float)[1:]                 # the bit before ...
    y = fold.b.astype(float)[1:]                     # ... and the bit itself, for p = 1..n-1
    ok = np.ones(n - 1) if fold.mask is None else (fold.mask[1:] & fold.mask[:-1]).astype(float)
    col = np.arange(1, n) % w

    def colsum(v: np.ndarray) -> np.ndarray:
        return np.bincount(col, weights=v * ok, minlength=w)

    m = colsum(np.ones(n - 1))
    sx, sy, sd = colsum(x), colsum(y), colsum(np.abs(x - y))
    d = sd / np.maximum(m, 1.0)                      # fraction of rows in which columns c-1 and c differ
    # A column pair is informative only if both columns vary (at least 5 ones and 5 zeros in each): a vertical
    # line next to a blank column differs in every row without being a seam.
    varied = (np.minimum(sx, m - sx) >= 5) & (np.minimum(sy, m - sy) >= 5)
    if varied.sum() < 3:
        return 0, 0.0
    c = int(np.flatnonzero(varied)[np.argmax(d[varied])])
    med = float(np.median(d))
    p0 = min(max(med, 1.0 / m[c]), 0.5)
    z = float((d[c] - med) / math.sqrt(p0 * (1 - p0) / m[c]))
    return c, z


def _header(b: np.ndarray, mask: np.ndarray | None, width: int, min_saving: float = 10.0) -> tuple[int, float]:
    """Estimated length of an incompressible prefix, in bits, and the bits it saves (see FoldResult.header_bits).

    A two-part code: the first h bits sent raw (1 bit each, erased bits free), then the rest of the stream coded
    afresh by the 2-D code at `width`, as if the raster started at bit h. The h that minimises the total is the
    estimate; it is reported only if it saves at least `min_saving` bits over h = 0 after paying log2(number of h
    tried) to say which h. By the no-hypercompression argument the chance of a saving of s bits from a stream with
    no header (under the 2-D code as null) is at most 2^-s.
    """
    n = b.size
    h_max = min(n // 2, max(8 * width, 256), 4096)
    if h_max < 2:
        return 0, 0.0
    raw = np.concatenate([[0], np.cumsum(np.ones(n) if mask is None else mask.astype(float))])
    total = np.empty(h_max)
    for h in range(h_max):
        rest = _Folder(b[h:], None if mask is None else mask[h:], width)
        total[h] = raw[h] + rest.folded(width)
    h = int(np.argmin(total))
    saving = float(total[0] - total[h]) - math.log2(h_max)
    if h == 0 or saving < min_saving:
        return 0, max(saving, 0.0) if h else 0.0
    return h, saving


# ---------------------------------------------------------------------------------------------------------------
# The width scan.
# ---------------------------------------------------------------------------------------------------------------


@dataclass
class FoldResult:
    """The outcome of `fold_scan`.

    width: the chosen row length. Start from the width with the largest gain; move to W-1 or W+1 if the pure
        vertical code (the bit at i - W as the only context) prefers it, because the diagonal context bits let the
        neighbours of the true width rival it; then take the smallest divisor of that width whose gain is at least
        90% of the largest (a raster also folds at multiples of its width, a periodic stream at every multiple of
        its period).
    gain_bits: baseline_bits - folded_bits - log2(number of widths tried), at `width`. Positive means the 2-D code
        saves that many bits net of paying for the width search.
    p_bound: min(1, 2^-gain_bits): a conservative probability of a gain this large under the null (see fold_scan).
    n_bits, n_rows (= n_bits // width), complete (n_bits % width == 0).
    baseline_bits, baseline_order: the null code length and its chosen memory m.
    folded_bits: the 2-D code length at `width` (best of the four contexts of CONTEXTS, including 2 bits for the
        choice; meta["context"] names the winner).
    widths, gains: the whole scan (gain at every width tried, same definition as gain_bits).
    seam: the column (in the stream's own phase, offset 0) at which a picture row really starts: header_bits mod
        width when a header was found, otherwise 0 (no evidence). Detection does not depend on it: the 2-D code is
        the same whatever the phase.
    seam_contrast: the left-neighbour discontinuity statistic, kept for audit but not used for `seam`: for each
        column c, the fraction of rows in which the bit at c differs from the bit before it (column c-1, or the
        previous row's last bit for c = 0), over column pairs that both vary; seam_contrast is the excess of the
        largest such fraction over the median across columns in binomial standard errors (a z-score), and
        meta["seam_candidate"] is that column. On the benchmark pictures it points at picture content (the edge of
        a filled square, a border line next to a blank column) rather than at the row start, and a picture with a
        blank border has no seam signature at all, which is why it does not set `seam`.
    header_bits: the estimated length of an incompressible prefix before the raster, in bits (0 if none or if the
        stream was not detected): the h that minimises "h raw bits + the 2-D code of the rest of the stream coded
        from bit h", reported when it saves at least 10 bits net of log2(number of h tried); meta
        ["header_saving_bits"] holds the saving. Random header bits that happen to match the picture's background
        are absorbed into the picture, so the estimate can fall a few bits short.
    factor_note: the factorisation of n_bits, for example "1679 = 23 x 73: both prime; width 23 divides it".
    surrogate_p: (1 + #{surrogate max gain >= observed max gain}) / (n_surrogates + 1), or nan if not run.
    detected: gain_bits >= the threshold (default 20 bits, p_bound <= 1e-6).
    """

    width: int
    gain_bits: float
    p_bound: float
    n_bits: int
    n_rows: int
    complete: bool
    baseline_bits: float
    baseline_order: int
    folded_bits: float
    widths: np.ndarray
    gains: np.ndarray
    seam: int = 0
    seam_contrast: float = 0.0
    header_bits: int = 0
    factor_note: str = ""
    surrogate_p: float = float("nan")
    detected: bool = False
    meta: dict[str, Any] = field(default_factory=dict)

    def gain_at(self, width: int) -> float:
        """The gain at one width of the scan (nan if it was not tried)."""
        hit = np.flatnonzero(self.widths == width)
        return float(self.gains[hit[0]]) if hit.size else float("nan")


def _scan_gains(fold: _Folder, widths: np.ndarray, base: float) -> tuple[np.ndarray, np.ndarray]:
    """Gain at every width (net of the width search and the context choice) and the index of the context used."""
    penalty = math.log2(len(widths))
    c = fold.costs_many(widths)
    which = np.argmin(c, axis=1)
    costs = c[np.arange(len(widths)), which] + _CHOICE_BITS
    return base - costs - penalty, which


def _p_bound(gain: float) -> float:
    return 1.0 if gain <= 0 else float(2.0 ** -gain)


def _choose(fold: _Folder, lags: np.ndarray, gains: np.ndarray) -> int:
    """Index of the chosen lag (row length, or plane size for a volume) given the gain at every lag.

    Start from the largest gain. The neighbouring lags L-1 and L+1 see the true "up" (or "back") bit through
    their diagonal context bits, so they can rival L; the pure vertical code (context: the bit at i - L alone)
    tells them apart, so the best of L-1, L, L+1 under that code is taken. Then, because a raster also folds at
    multiples of its width, the smallest divisor of that lag (among those scanned) that keeps 90% of the largest
    gain is preferred.
    """
    j = int(np.argmax(gains))
    g_max = float(gains[j])
    if g_max <= 0:
        return j
    best, best_cost = j, np.inf
    for k in (j - 1, j, j + 1):
        if 0 <= k < len(lags):
            cost = _context_cost(_counts(fold.sh(int(lags[k])), fold.b, fold.mask, 2))
            if cost < best_cost:
                best, best_cost = k, cost
    l_star = int(lags[best])
    for k in np.flatnonzero((l_star % lags == 0) & (gains >= 0.9 * g_max)):
        return int(k)
    return best


def fold_scan(stream: BitStream, min_width: int = 2, max_width: int | None = None,
              rng: np.random.Generator | None = None, n_surrogates: int = 0, threshold_bits: float = 20.0
              ) -> FoldResult:
    """Scan row lengths and report the best raster width, its compression gain and a probability bound.

    For every width W in min_width..max_width (default min(n // 4, 4096), so at least four rows) the 2-D context
    code is compared with the null code (see `baseline_bits` and `folded_bits`), and
    gain(W) = baseline - folded(W) - log2(number of widths tried).

    What p_bound means. The null is the mixture over short-memory binary sources that the baseline code
    represents: the baseline code length L0 is -log2 of a probability distribution Q0 (KT mixtures over the
    order-m Markov sources, mixed uniformly over m = 0..8). The 2-D codes mixed uniformly over the widths are a
    second prefix code L1. By the no-hypercompression inequality, P_Q0(L0 - L1 >= g) <= 2^-g, so p_bound =
    min(1, 2^-gain) is valid for the whole width search with no further correction. For one fixed source of the
    null family (say a particular Markov chain of order m) the bound holds only up to the baseline's regret
    against that source, about 2^(m-1) log2 n + log2 9 bits (a few to a few tens of bits for small m), so treat
    gains below about 20 bits as unremarkable. The bound says nothing about sources with memory longer than 8
    bits but no geometry (text, compressed data with long matches, a strictly periodic stream): those also gain,
    at their period or a multiple, which is real structure but not necessarily a picture.

    The surrogate check (n_surrogates > 0) simulates Markov surrogates of order min(baseline order, 3) with
    `markov_surrogate` (erased positions keep the same mask) and recomputes the maximum gain; the empirical p-value
    is reported as surrogate_p. It is a check on the bound under a fitted, not a mixed, null.

    It cannot see the scan order or the orientation: column-major order folds at the height instead of the width.
    """
    b, mask = _prepare(stream.bits, stream.known)
    n = b.size
    if max_width is None:
        max_width = min(n // 4, 4096)
    max_width = int(min(max_width, max(n // 2, 2)))
    min_width = int(max(min_width, 2))
    if max_width < min_width:
        raise ValueError(f"stream of {n} bits is too short to scan widths {min_width}..{max_width}")
    widths = np.arange(min_width, max_width + 1)
    base, order = baseline_bits(b, mask)
    fold = _Folder(b, mask, max_width)
    gains, which = _scan_gains(fold, widths, base)
    g_max = float(gains.max())
    j = _choose(fold, widths, gains)
    width = int(widths[j])
    gain = float(gains[j])
    folded = base - gain - math.log2(len(widths))
    seam_candidate, contrast = _seam(fold, width)
    header, header_saving = (_header(b, mask, width) if gain >= threshold_bits else (0, 0.0))
    seam = header % width if header else 0
    surrogate_p = float("nan")
    if n_surrogates > 0:
        from stargazer.structure.surrogates import markov_surrogate
        rng = rng if rng is not None else rng_for(0, "decode.layout", "fold_scan", "surrogates")
        hits = 0
        for _ in range(n_surrogates):
            s = markov_surrogate(b, 2, min(order, 3), rng).astype(np.uint8)
            s_base, _ = baseline_bits(s, mask)
            s_gains, _ = _scan_gains(_Folder(s, mask, max_width), widths, s_base)
            hits += int(s_gains.max() >= g_max)
        surrogate_p = (1.0 + hits) / (n_surrogates + 1.0)
    return FoldResult(width=width, gain_bits=gain, p_bound=_p_bound(gain), n_bits=n,
                      n_rows=n // width, complete=(n % width == 0), baseline_bits=base, baseline_order=order,
                      folded_bits=folded, widths=widths, gains=gains, seam=seam, seam_contrast=contrast,
                      header_bits=header, factor_note=_factor_note(n, width), surrogate_p=surrogate_p,
                      detected=bool(gain >= threshold_bits),
                      meta={"context": MODELS[int(which[j])], "max_gain_bits": g_max,
                            "max_gain_width": int(widths[int(np.argmax(gains))]), "seam_candidate": seam_candidate,
                            "header_saving_bits": header_saving})


# ---------------------------------------------------------------------------------------------------------------
# The plane-height scan.
# ---------------------------------------------------------------------------------------------------------------


@dataclass
class VolumeResult:
    """The outcome of `volume_scan`.

    width: the fastest axis (given). height: the chosen plane height H, by the same rule as FoldResult.width
        applied to the plane size W*H (largest gain, then the neighbour preferred by the pure "back" code, then
        the smallest divisor keeping 90% of the largest gain).
    depth_layers: n_bits // (width * height). plane_bits: width * height.
    gain_bits: folded_bits(width) - volume code(width, H) - log2(number of heights tried): the bits saved by also
        conditioning on the bit one plane back. p_bound = min(1, 2^-gain_bits), under the null that the 2-D code at
        `width` represents (a picture with no correlation between planes), by the same argument as in fold_scan.
    detected: gain_bits >= threshold. heights, gains: the whole scan.
    """

    width: int
    height: int
    depth_layers: int
    plane_bits: int
    gain_bits: float
    p_bound: float
    detected: bool
    heights: np.ndarray
    gains: np.ndarray


def volume_scan(stream: BitStream, width: int, min_height: int = 2, max_height: int | None = None,
                threshold_bits: float = 20.0) -> VolumeResult:
    """Scan plane heights H for a 3-D raster whose fastest axis has length `width`.

    The volume code is the 2-D code at `width` (best of its four contexts, 2 bits for the choice) with one more
    context bit, back, the bit at i - width*H (one plane back); its gain is measured against the 2-D code at
    `width` and charged log2(number of heights tried), so p_bound covers the height search.
    max_height defaults to n // (2 * width), so at least two planes. A raster with a fixed row and plane length is
    assumed; it cannot see a volume whose fastest axis is not `width`, and a 2-D raster whose rows repeat every H
    rows (a periodic picture) also gains here.
    """
    b, mask = _prepare(stream.bits, stream.known)
    n = b.size
    if width < 2:
        raise ValueError("width must be at least 2")
    if max_height is None:
        max_height = n // (2 * width)
    min_height = max(int(min_height), 2)
    if max_height < min_height:
        raise ValueError("stream too short for a volume of this width")
    heights = np.arange(min_height, max_height + 1)
    fold = _Folder(b, mask, width * max_height)
    w = width
    two_d = fold.folded(w)
    penalty = math.log2(len(heights))
    costs = fold.costs_many(np.full(heights.size, w), back=heights * w).min(axis=1) + _CHOICE_BITS
    gains = two_d - costs - penalty
    j = _choose(fold, heights * w, gains)
    h = int(heights[j])
    gain = float(gains[j])
    return VolumeResult(width=int(width), height=h, depth_layers=n // (w * h), plane_bits=w * h, gain_bits=gain,
                        p_bound=_p_bound(gain), detected=bool(gain >= threshold_bits),
                        heights=heights, gains=gains)


# ---------------------------------------------------------------------------------------------------------------
# Rendering.
# ---------------------------------------------------------------------------------------------------------------


def render(stream: BitStream, width: int, offset: int = 0, height: int | None = None) -> np.ndarray:
    """The stream as a picture: rows of `width` bits starting at bit `offset`, last partial row dropped.

    Returns a uint8 array of shape (rows, width); if the stream has erasures, a numpy masked array whose mask is
    True at erased pixels. `height` limits the number of rows.
    """
    if width < 1:
        raise ValueError("width must be positive")
    bits = stream.bits[offset:]
    rows = bits.size // width
    if height is not None:
        rows = min(rows, int(height))
    img = bits[:rows * width].reshape(rows, width).astype(np.uint8)
    if stream.known is not None and not stream.known.all():
        erased = ~stream.known[offset:offset + rows * width].reshape(rows, width)
        return np.ma.MaskedArray(img, mask=erased)
    return img


def to_text(image: np.ndarray, on: str = "#", off: str = ".", erased: str = "?") -> str:
    """A picture as lines of text, one character per pixel (masked pixels as `erased`)."""
    data = np.asarray(np.ma.getdata(image))
    mask = np.ma.getmaskarray(image) if np.ma.isMaskedArray(image) else np.zeros(data.shape, dtype=bool)
    lines = []
    for r in range(data.shape[0]):
        lines.append("".join(erased if mask[r, c] else (on if data[r, c] else off) for c in range(data.shape[1])))
    return "\n".join(lines)
