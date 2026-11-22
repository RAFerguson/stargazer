"""Blind recognition of engineered redundancy in a bit stream: frames and sync words, linear block codes, LFSRs.

A sender can make a payload look like noise by compressing or encrypting it, but a receiver still has to find frame
boundaries and correct errors, so the framing and error-correction "envelope" stays visible. The four recognisers here
look for the three most common kinds of envelope that communication engineering adds on purpose:

* :func:`frame_scan` - a repeating frame length, from the lag-agreement spectrum (a sync word repeated every frame
  raises agreement at the frame length and all its multiples);
* :func:`sync_word` - the positions of a frame that stay constant (the synchronisation word) and their pattern;
* :func:`block_code_scan` - parity checks of a linear block code, by rank deficiency and a noise-tolerant
  random-subset null-space search (in the line of Sicot-Houcke-Barbier and Valembois; see the research note,
  section 2);
* :func:`lfsr_test` - a linear-feedback shift register (additive scrambler or PN sequence), by Berlekamp-Massey linear
  complexity, globally and on short windows.

Every result carries an estimate, a test statistic and a p-value that is corrected for the parameter search the
function performed; each docstring says how. "Surrogates" are always streams of independent bits with the same
fraction of ones as the input (and the same erasure mask), drawn from the generator passed in. When no generator is
passed, a fixed one from :func:`stargazer.core.rng.rng_for` is used so that results are reproducible.

All of this assumes GF(2)-linear, human-style coding. It cannot see non-linear codes or scramblers, variable-length
frames, bit slips, or codes whose parity checks are too heavy for the search (see :func:`block_code_scan`).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import fft as sfft
from scipy import stats

from stargazer.core.rng import rng_for
from stargazer.core.types import BitStream

from . import gf2

__all__ = [
    "FrameResult", "SyncResult", "BlockCodeResult", "LfsrResult",
    "frame_scan", "sync_word", "block_code_scan", "lfsr_test", "describe",
]

_LOG10 = np.log(10.0)


# --------------------------------------------------------------------------------------------------------- helpers

def _bits_known(stream: BitStream | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(stream, BitStream):
        stream = BitStream(np.asarray(stream))
    bits = stream.bits.astype(np.uint8)
    known = np.ones(bits.size, dtype=bool) if stream.known is None else stream.known.astype(bool)
    return bits, known


def _ones_fraction(bits: np.ndarray, known: np.ndarray) -> float:
    k = int(known.sum())
    return float(bits[known].sum()) / k if k else 0.5


def _rng(rng: np.random.Generator | None, name: str) -> np.random.Generator:
    return rng if rng is not None else rng_for(0, "decode.codes", name)


def _surrogate_bits(n: int, p: float, rng: np.random.Generator) -> np.ndarray:
    return (rng.random(n) < p).astype(np.uint8)


def _log_binom_tail(k: np.ndarray | int, n: np.ndarray | int, q: np.ndarray | float) -> np.ndarray:
    """Natural log of P(X >= k) for X ~ Binomial(n, q), accurate far into the tail.

    scipy's logsf is used where it is finite; where it underflows, the tail is bounded by the geometric series
    pmf(k) / (1 - r) with r = pmf(k+1)/pmf(k) < 1, which is tight in the far tail.
    """
    k_a, n_a, q_a = np.broadcast_arrays(np.asarray(k, float), np.asarray(n, float), np.asarray(q, float))
    out = np.zeros(k_a.shape, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        direct = stats.binom.logsf(k_a - 1, n_a, q_a)
        logpmf = stats.binom.logpmf(k_a, n_a, q_a)
        ratio = (n_a - k_a) / (k_a + 1) * q_a / np.maximum(1 - q_a, 1e-300)
        bound = logpmf - np.log1p(-np.minimum(ratio, 1 - 1e-12))
    good = np.isfinite(direct)
    out[good] = direct[good]
    out[~good] = bound[~good]
    out[k_a <= 0] = 0.0
    return np.minimum(out, 0.0)


def _p_from_log(log_p: float) -> float:
    return float(np.exp(min(log_p, 0.0)))


# ------------------------------------------------------------------------------------------------------ frame scan

@dataclass
class FrameResult:
    """Outcome of :func:`frame_scan`.

    frame_len is the fundamental frame length (the lag whose comb of multiples stands out most); z_score is the
    agreement z-score at that lag and comb_statistic the comb statistic that was maximised. p_value is corrected for
    the search over all lags (by surrogates when n_surrogates > 0, otherwise Bonferroni); p_bonferroni is always the
    analytic Bonferroni bound. lags/agreement/expected/z/comb are the full spectrum.
    """

    frame_len: int
    z_score: float
    comb_statistic: float
    p_value: float
    p_bonferroni: float
    method: str
    lags: np.ndarray
    agreement: np.ndarray
    expected: float
    z: np.ndarray
    comb: np.ndarray
    n_pairs: np.ndarray
    harmonic_lags: np.ndarray
    harmonic_z: np.ndarray
    n_surrogates: int
    surrogate_max: np.ndarray = field(default_factory=lambda: np.zeros(0))


def _agreement_counts(bits: np.ndarray, known: np.ndarray, max_lag: int) -> tuple[np.ndarray, np.ndarray]:
    """(agreeing pairs, known pairs) for lags 0..max_lag, all at once by FFT autocorrelation.

    bits may be one stream (n,) or a batch of streams (S, n) sharing the known mask; the result has the same
    leading shape as bits, with max_lag + 1 lags on the last axis.
    """
    n = bits.shape[-1]
    nfft = sfft.next_fast_len(2 * n)
    x = np.where(known, 2.0 * bits - 1.0, 0.0)
    xs = sfft.irfft(np.abs(sfft.rfft(x, nfft, axis=-1)) ** 2, nfft, axis=-1)[..., : max_lag + 1]
    ks = sfft.irfft(np.abs(sfft.rfft(known.astype(float), nfft)) ** 2, nfft)[: max_lag + 1]
    pairs = np.rint(ks)
    agree = np.rint((pairs + np.rint(xs)) / 2.0)
    return agree, np.broadcast_to(pairs, agree.shape)


def _comb_index(lags: np.ndarray, lag_max: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flattened indices (into the lag array) of the multiples of every lag, with group starts and group sizes."""
    first = int(lags[0])
    pieces, starts, sizes = [], [], []
    pos = 0
    for lag in lags:
        mult = np.arange(lag, lag_max + 1, lag) - first
        pieces.append(mult)
        starts.append(pos)
        sizes.append(mult.size)
        pos += mult.size
    return np.concatenate(pieces), np.asarray(starts), np.asarray(sizes, dtype=float)


def _frame_spectrum(bits: np.ndarray, known: np.ndarray, lags: np.ndarray, comb_idx: tuple[np.ndarray, ...]
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Agreement fraction, z-score, null mean, comb statistic and pair counts per lag (batched over leading axis)."""
    n = bits.shape[-1]
    n_known = max(int(known.sum()), 1)
    p = (np.where(known, bits, 0).sum(axis=-1, keepdims=True) / n_known).astype(float)
    q = 1.0 - p
    agree, pairs = _agreement_counts(bits, known, int(lags[-1]))
    agree, pairs = agree[..., lags], np.maximum(pairs[..., lags], 1.0)
    frac = agree / pairs
    mu = p * p + q * q
    # Indicators [b_i == b_{i+L}] and [b_{i+L} == b_{i+2L}] share a bit; their covariance is p^3 + q^3 - mu^2
    # (zero for p = 1/2). The number of such chained pairs is about pairs * (n - 2L) / (n - L).
    cov = p ** 3 + q ** 3 - mu * mu
    chained = pairs * np.clip(n - 2.0 * lags, 0, None) / np.maximum(n - lags, 1.0)
    var = (pairs * mu * (1 - mu) + 2.0 * chained * cov) / pairs ** 2
    z = (frac - mu) / np.sqrt(np.maximum(var, 1e-300))
    idx, starts, sizes = comb_idx
    comb = np.add.reduceat(z[..., idx], starts, axis=-1) / np.sqrt(sizes)
    return frac, z, mu, comb, pairs


def frame_scan(stream: BitStream, min_len: int = 8, max_len: int = 4096, rng: np.random.Generator | None = None,
               n_surrogates: int = 200) -> FrameResult:
    """Find a repeating frame length from the lag-agreement spectrum.

    Algorithm. For every lag L in [min_len, max_len] the agreement a(L) is the fraction of positions i, with both bits
    known, where bit[i] == bit[i + L]; all lags come from one FFT autocorrelation of the +-1 sequence and of the known
    mask. Under independence with ones-fraction p, a(L) has mean mu = p^2 + (1-p)^2 and the variance of a sum of
    pairwise-dependent indicators (exact for an unerased stream, approximate with erasures); z(L) is the standardised
    agreement. A sync word repeated every F bits raises a(F), a(2F), a(3F) ...; the comb statistic
    S(L) = sum_{j: jL <= max_len} z(jL) / sqrt(J) rewards a lag whose multiples all stand out, and the fundamental frame
    length is argmax S(L) (S of a sub-multiple F/m is diluted by the m-1 silent lags between multiples, and S of a
    multiple mF uses fewer terms).

    Null hypothesis and p-value. H0: independent bits with the stream's ones-fraction. The test statistic is
    max_L S(L); p_value is corrected for the search over every lag by rerunning the whole scan on n_surrogates
    surrogate streams (p = (1 + #{surrogate max >= observed}) / (1 + n_surrogates), so it cannot go below
    1/(1+n_surrogates)). p_bonferroni = n_lags * P(N(0,1) > max S) is also reported (approximate: z is Gaussian only
    for many pairs). With n_surrogates = 0 the Bonferroni value is used as p_value.

    Limits. Needs several frames (the excess agreement is (sync length / frame length), so a 32-bit word in a 512-bit
    frame needs roughly 10^4 bits). Variable-length frames, bit slips and a scrambler that covers the sync word
    defeat it. Lags beyond len(stream) - 64 are not scanned.
    """
    bits, known = _bits_known(stream)
    n = bits.size
    lag_max = int(min(max_len, n - 64))
    if lag_max < min_len:
        raise ValueError(f"stream of {n} bits is too short for lags {min_len}..{max_len}")
    lags = np.arange(min_len, lag_max + 1)
    comb_idx = _comb_index(lags, lag_max)
    frac, z, mu_arr, comb, pairs = _frame_spectrum(bits, known, lags, comb_idx)
    mu = float(mu_arr[0])
    best = int(np.argmax(comb))
    stat = float(comb[best])
    frame_len = int(lags[best])
    p_bonf = float(min(1.0, lags.size * stats.norm.sf(stat)))

    surrogate_max = np.zeros(0)
    if n_surrogates > 0:
        g = _rng(rng, "frame_scan")
        p1 = _ones_fraction(bits, known)
        chunks = []
        for start in range(0, n_surrogates, 25):
            sb = (g.random((min(25, n_surrogates - start), n)) < p1).astype(np.uint8)
            chunks.append(_frame_spectrum(sb, known, lags, comb_idx)[3].max(axis=-1))
        surrogate_max = np.concatenate(chunks)
        p_value = float((1 + np.sum(surrogate_max >= stat)) / (1 + n_surrogates))
        method = f"surrogates ({n_surrogates} independent-bit streams, whole lag search rerun)"
    else:
        p_value, method = p_bonf, f"Bonferroni over {lags.size} lags (Gaussian z)"

    harmonic_lags = np.arange(frame_len, lag_max + 1, frame_len)
    return FrameResult(frame_len=frame_len, z_score=float(z[best]), comb_statistic=stat, p_value=p_value,
                       p_bonferroni=p_bonf, method=method, lags=lags, agreement=frac, expected=float(mu), z=z,
                       comb=comb, n_pairs=pairs, harmonic_lags=harmonic_lags,
                       harmonic_z=z[harmonic_lags - lags[0]], n_surrogates=int(n_surrogates),
                       surrogate_max=surrogate_max)


# ------------------------------------------------------------------------------------------------------- sync word

@dataclass
class SyncResult:
    """Outcome of :func:`sync_word`.

    phase is the residue (stream index modulo frame_len) where the constant run starts; positions are the residues in
    the run that tested constant, pattern the majority bits over the whole run (offset 0 = phase), and
    pattern_positions the residues pattern refers to. agreement[r] is the majority fraction at residue r and
    position_p[r] its two-sided binomial p-value. run_length counts constant positions in the run, run_span its
    length including tolerated gaps. p_value is the corrected probability of so long a run under independence.
    """

    frame_len: int
    phase: int
    positions: np.ndarray
    pattern: np.ndarray
    pattern_positions: np.ndarray
    agreement: np.ndarray
    position_p: np.ndarray
    constant: np.ndarray
    run_length: int
    run_span: int
    p_value: float
    log10_p: float
    background_p: float
    n_frames: int
    tau: float

    @property
    def pattern_hex(self) -> str:
        if self.pattern.size == 0:
            return ""
        value = int("".join(str(int(b)) for b in self.pattern), 2)
        return f"0x{value:0{-(-self.pattern.size // 4)}X}"


def _majority_pvalues(ones: np.ndarray, total: np.ndarray, p: float) -> np.ndarray:
    """Two-sided P(max(X, t - X) >= observed majority) for X ~ Binomial(t, p), per position."""
    maj = np.maximum(ones, total - ones)
    minority = total - maj
    upper = np.exp(_log_binom_tail(maj, total, p))                    # P(X >= maj)
    lower = stats.binom.cdf(minority, total, p)                         # P(X <= t - maj)
    out = np.where(maj > minority, upper + lower, 1.0)
    out = np.where(total > 0, out, 1.0)
    return np.minimum(out, 1.0)


def _longest_cyclic_run(flags: np.ndarray, max_gap: int) -> tuple[int, int, int]:
    """(start, span, count) of the cyclic segment holding the most flagged positions, bridging gaps <= max_gap."""
    size = flags.size
    where = np.flatnonzero(flags)
    if where.size == 0:
        return 0, 0, 0
    if where.size == size:
        return 0, size, size
    step = np.diff(np.concatenate([where, [where[0] + size]]))
    breaks = np.flatnonzero(step > max_gap + 1)
    if breaks.size == 0:
        return int(where[0]), size, int(where.size)
    rot = (breaks[0] + 1) % where.size
    order = np.roll(where, -rot)
    gaps = np.roll(step, -rot)
    best = (0, 0, 0)
    i = 0
    while i < order.size:
        j = i
        while j < order.size - 1 and gaps[j] <= max_gap + 1:
            j += 1
        count = j - i + 1
        span = int((order[j] - order[i]) % size) + 1
        if count > best[2] or (count == best[2] and span < best[1]):
            best = (int(order[i]), span, count)
        i = j + 1
    return best


def sync_word(stream: BitStream, frame_len: int, tau: float = 1e-3, max_gap: int = 1) -> SyncResult:
    """Recover the synchronisation word of frames of known length.

    Algorithm. Cutting the stream into frames at phase 0..frame_len-1 and stacking them is the same as grouping the
    bits by residue r = index mod frame_len; at every residue we count ones and known bits. The background
    ones-fraction p_bg is the median over residues (most residues carry payload, so the median is what a residue's
    neighbours predict). A residue is *constant* when its majority count is improbable under Binomial(t, p_bg):
    two-sided p-value below tau. Bit errors only lower the majority fraction, so a 32-bit word at a few per cent BER
    still tests constant with a few dozen frames. The sync word is the cyclic run of residues holding the most
    constant positions, bridging gaps of up to max_gap non-constant residues (a position spoiled by errors); the phase
    is where the run starts, and the pattern is the majority bit at every residue of the run.

    Null hypothesis and p-value. H0: independent bits with ones-fraction p_bg. Each residue is flagged with
    probability at most tau, independently, so a run of span s holding c flagged residues has tail probability
    P(Binomial(s, tau) >= c). p_value multiplies this by frame_len (possible start phases) and by frame_len again
    (possible spans): Bonferroni over the run search.

    Limits. Needs the right frame length (see frame_scan); a sync word that changes between frames (counters,
    alternating words) shows only its constant bits; a header field that is constant also appears as "sync".
    """
    if frame_len < 2:
        raise ValueError("frame_len must be at least 2")
    bits, known = _bits_known(stream)
    residue = np.arange(bits.size) % frame_len
    total = np.bincount(residue[known], minlength=frame_len).astype(float)
    ones = np.bincount(residue[known], weights=bits[known].astype(float), minlength=frame_len)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac1 = np.where(total > 0, ones / np.maximum(total, 1), 0.5)
    p_bg = float(np.clip(np.median(frac1), 1e-6, 1 - 1e-6))
    pos_p = _majority_pvalues(ones, total, p_bg)
    agreement = np.where(total > 0, np.maximum(ones, total - ones) / np.maximum(total, 1), 0.0)
    constant = pos_p <= tau
    start, span, count = _longest_cyclic_run(constant, max_gap)
    run_res = (start + np.arange(span)) % frame_len
    pattern = (ones[run_res] * 2 > total[run_res]).astype(np.uint8)
    positions = run_res[constant[run_res]]
    if count > 0:
        log_p = float(_log_binom_tail(count, span, tau)) + 2.0 * np.log(frame_len)
    else:
        log_p = 0.0
    log_p = min(log_p, 0.0)
    return SyncResult(frame_len=frame_len, phase=int(start), positions=positions, pattern=pattern,
                      pattern_positions=run_res, agreement=agreement, position_p=pos_p, constant=constant,
                      run_length=int(count), run_span=int(span), p_value=_p_from_log(log_p), log10_p=log_p / _LOG10,
                      background_p=p_bg, n_frames=int(bits.size // frame_len), tau=tau)


# ---------------------------------------------------------------------------------------------- block code scan

@dataclass
class BlockCodeResult:
    """One detected linear block code at block length n and alignment offset.

    H holds the n - k independent parity checks (rows), lowest weight first; satisfaction and weights describe each
    row, check_log10_p its uncorrected binomial tail on all blocks. ber_estimate is the channel bit-error rate implied
    by the satisfaction fractions. p_value (and log10_p) is corrected for the candidate checks tried at this
    (n, offset) and for every (n, offset) pair scanned. evidence = sum over rows of H of -check_log10_p (the ranking
    key). exact_rank is the GF(2) rank of all blocks (rank < n means the code was seen without errors). collapsed
    lists the (n, offset) results that were folded into this one because their n is a multiple of this n, and
    aliases those folded in because they are partial echoes of this code at a non-multiple length (see _collapse).
    """

    n: int
    k: int
    offset: int
    H: np.ndarray
    satisfaction: np.ndarray
    weights: np.ndarray
    check_log10_p: np.ndarray
    ber_estimate: float
    p_value: float
    log10_p: float
    evidence: float
    n_blocks: int
    n_candidates: int
    exact_rank: int
    n_pairs_scanned: int
    collapsed: list[tuple[int, int]] = field(default_factory=list)
    aliases: list[tuple[int, int]] = field(default_factory=list)


def _null_sat_prob(weight: np.ndarray, p1: float) -> np.ndarray:
    """P(h . block = 0) for a weight-w check on independent bits with ones-fraction p1: (1 + (1 - 2 p1)^w) / 2."""
    return 0.5 * (1.0 + (1.0 - 2.0 * p1) ** np.asarray(weight, dtype=float))


def _ber_from_checks(sat: np.ndarray, weights: np.ndarray) -> float:
    """Pooled estimate of e from (1 - 2e)^w = 2 f - 1 over all checks."""
    excess = np.clip(2.0 * np.asarray(sat, float) - 1.0, 1e-12, 1.0)
    rho = float(np.exp(np.sum(np.log(excess)) / max(float(np.sum(weights)), 1.0)))
    return float(max(0.0, (1.0 - rho) / 2.0))


def _independent_subset(checks: list[int]) -> list[int]:
    """Greedy GF(2)-independent subset of the given checks, in the given order (Python-int Gaussian basis)."""
    basis: dict[int, int] = {}   # leading bit -> vector
    kept: list[int] = []
    for h in checks:
        v = h
        while v:
            top = v.bit_length() - 1
            if top in basis:
                v ^= basis[top]
            else:
                basis[top] = v
                kept.append(h)
                break
    return kept


def _span(vectors: list[int]) -> list[int]:
    out = [0]
    for v in vectors:
        out += [x ^ v for x in out]
    return out[1:]


def _blocks_for(bits: np.ndarray, known: np.ndarray, n: int, offset: int, max_blocks: int) -> np.ndarray:
    n_blocks = min((bits.size - offset) // n, max_blocks)
    if n_blocks <= 0:
        return np.zeros(0, dtype=np.uint64)
    seg = bits[offset: offset + n_blocks * n].reshape(n_blocks, n).astype(np.uint64)
    ok = known[offset: offset + n_blocks * n].reshape(n_blocks, n).all(axis=1)
    weights = np.uint64(1) << np.arange(n, dtype=np.uint64)
    return (seg[ok] * weights).sum(axis=1).astype(np.uint64)


def block_code_scan(stream: BitStream, n_min: int = 3, n_max: int = 64, rng: np.random.Generator | None = None,
                    max_blocks: int = 4000, n_trials: int = 16, extra_rows: int = 4, alpha: float = 0.01,
                    span_dim_max: int = 6) -> list[BlockCodeResult]:
    """Search for the parity checks of a linear block code at every block length and alignment.

    Algorithm. For each block length n in [n_min, n_max] (n <= 64, one machine word per block) and each offset
    0..n-1 the stream is cut into consecutive n-bit blocks (at most max_blocks; blocks with erased bits are dropped).
    A parity check h of a linear code satisfies h . block = 0 for every codeword, so the checks are the null space of
    the block matrix. Errors give the full matrix full rank, so the search is robust (random-subset method):
    n_trials times, draw a random subset of m = n + extra_rows blocks and take its null space (all subsets of all
    offsets of one n are eliminated together, see gf2.batched_null_space); an error-free subset that spans the code
    yields exactly the dual code, one with an error yields a subspace of it, and a subset of random blocks has an
    empty null space except with probability about 2^-extra_rows. When the null space has dimension <= span_dim_max
    every non-zero member is a candidate, else only its basis. Every distinct candidate is scored on the blocks
    *outside* the subset that produced it by the fraction s satisfying h . block = 0. A true check of weight w on a
    channel with bit-error rate e is satisfied with probability (1 + (1-2e)^w)/2; under H0 the probability is
    q0(w) = (1 + (1-2p)^w)/2 for ones-fraction p (1/2 for balanced bits). Checks that pass (below) are completed by
    the low-weight members of their span, reduced to an independent set lowest weight first, and k = n - #checks.
    When the data are error free this is exactly k = rank (exact_rank is reported for comparison).

    Null hypothesis and p-values. H0: independent bits with the stream's ones-fraction. A candidate is chosen from
    its subset only, so on the other blocks its satisfaction count is Binomial(N - m, q0(w)) under H0; the one-sided
    tail is Bonferroni-corrected for the K candidates tried at that (n, offset) and for all (n, offset) pairs scanned
    (P = sum over n of n). A check is accepted when this fully corrected value is below alpha, and a result's
    p_value is the corrected value of its best check (pass alpha=1 to see every pair with any candidate). Results are
    sorted by evidence (summed -log10 tail of the independent checks, which approximates -log10 of the joint
    probability of the whole rank deficit), and a result whose n is a multiple of a smaller detected n (a (7,4) code
    also shows at 14, 21 ...) is folded into the best result of that smaller n and listed in its `collapsed` field;
    weaker echoes at other lengths, whose checks hold only at some phases of the stronger code within the block (a
    block of 47 bits holds whole (3,1) codewords at a phase that cycles from block to block), go to its `aliases`
    field (chi-square test of equal satisfaction across those phases; see _collapse).
    Codes of length n also produce weaker results at misaligned offsets of the same n (a block that straddles two
    codewords can still contain a check of one of them); those are kept and rank lower.

    What this cannot see. Non-linear codes; convolutional codes (no fixed block length; their checks appear only
    approximately at long n); block lengths above 64; and long modern codes - LDPC, turbo, polar - whose block
    lengths are hundreds to thousands of bits and whose checks are too heavy or too many for a random subset of
    n + extra_rows blocks to be error free at realistic bit-error rates. The robust search needs subsets free of
    errors (or nearly so): with block-error rate b, a trial is clean with probability about (1-b)^(n+extra_rows),
    so it degrades with n and BER. Measured with 14,000-bit streams and default settings (10 streams each):
    Hamming(7,4) found at BER up to 0.05; Hamming(15,11) at 0.01 always, 0.02 in 8/10, 0.05 never; (8,7) parity
    at 0.01 always, 0.02 in 8/10, 0.05 in 2/10. Interleaved codes appear only at the interleaved length.

    Cost. About 1 s for 14,000 bits and 2-3 s for 100,000 bits (max_blocks caps the work per offset), dominated by
    the batched elimination of n_max * (n_max + 1) / 2 * n_trials small matrices.
    """
    if n_max > 64:
        raise ValueError("block_code_scan supports block lengths up to 64")
    bits, known = _bits_known(stream)
    g = _rng(rng, "block_code_scan")
    p1 = _ones_fraction(bits, known)
    n_values = list(range(max(2, n_min), n_max + 1))
    n_pairs = int(sum(n_values))
    log_pairs = np.log(max(n_pairs, 1))
    log_alpha = np.log(alpha)
    found: list[BlockCodeResult] = []

    for n in n_values:
        m = n + extra_rows
        per_offset = [_blocks_for(bits, known, n, o, max_blocks) for o in range(n)]
        usable = [o for o in range(n) if per_offset[o].size >= m + 16]
        if not usable:
            continue
        n_max_blocks = max(per_offset[o].size for o in usable)
        padded = np.zeros((len(usable), n_max_blocks), dtype=np.uint64)
        keys = np.full((len(usable), n_trials, n_max_blocks), 2.0)
        for row, o in enumerate(usable):
            nb = per_offset[o].size
            padded[row, :nb] = per_offset[o]
            keys[row, :, :nb] = g.random((n_trials, nb))
        subsets = np.argpartition(keys, m - 1, axis=2)[:, :, :m]                          # (U, T, m)
        rows = np.take_along_axis(padded[:, None, :], subsets, axis=2).reshape(-1, m)
        basis, valid = gf2.batched_null_space(rows, n)
        dims = valid.sum(axis=1).reshape(len(usable), n_trials)
        basis = basis.reshape(len(usable), n_trials, n)
        valid = valid.reshape(len(usable), n_trials, n)

        for row, o in enumerate(usable):
            if not dims[row].any():
                continue
            blocks = per_offset[o]
            n_blocks = blocks.size
            candidates: dict[int, int] = {}
            for t in np.flatnonzero(dims[row]):
                vecs = [int(v) for v in basis[row, t][valid[row, t]]]
                members = _span(vecs) if len(vecs) <= span_dim_max else vecs
                for h in members:
                    candidates.setdefault(h, int(t))
            cand = np.fromiter(candidates.keys(), dtype=np.uint64, count=len(candidates))
            n_cand = cand.size
            weights = np.bitwise_count(cand).astype(int)
            sat_all = n_blocks - np.count_nonzero(gf2.parity64(blocks[None, :] & cand[:, None]), axis=1)
            sat_rest = sat_all - m
            q0 = _null_sat_prob(weights, p1)
            log_tail = _log_binom_tail(sat_rest, n_blocks - m, q0)
            log_corr = log_tail + np.log(n_cand) + log_pairs
            accepted = log_corr < log_alpha
            if not accepted.any():
                continue
            best_log = float(log_corr.min())
            acc = [int(h) for h in cand[accepted]]
            # Complete with low-weight members of the accepted span, scored on all blocks (used for reporting the
            # basis only; the p-value above comes from the held-out scores).
            pool = set(acc)
            acc_basis = _independent_subset(acc)
            if len(acc_basis) <= 10:
                pool.update(_span(acc_basis))
            pool_arr = np.fromiter(pool, dtype=np.uint64, count=len(pool))
            pw = np.bitwise_count(pool_arr).astype(int)
            psat = n_blocks - np.count_nonzero(gf2.parity64(blocks[None, :] & pool_arr[:, None]), axis=1)
            plog = _log_binom_tail(psat, n_blocks, _null_sat_prob(pw, p1))
            keep = (plog + np.log(n_cand + pool_arr.size) + log_pairs) < log_alpha
            order = np.lexsort((plog, pw))
            ordered = [int(pool_arr[i]) for i in order if keep[i]]
            indep = _independent_subset(ordered)
            if not indep:
                continue
            idx = {int(h): i for i, h in enumerate(pool_arr)}
            sel = np.array([idx[h] for h in indep])
            H = gf2.unpack_rows(np.array(indep, dtype=np.uint64)[:, None], n)
            sat_frac = psat[sel] / n_blocks
            exact = gf2.rank(gf2.unpack_rows(blocks[:, None], n))
            found.append(BlockCodeResult(
                n=n, k=n - len(indep), offset=o, H=H, satisfaction=sat_frac, weights=pw[sel],
                check_log10_p=plog[sel] / _LOG10, ber_estimate=_ber_from_checks(sat_frac, pw[sel]),
                p_value=_p_from_log(best_log), log10_p=min(best_log, 0.0) / _LOG10,
                evidence=float(-np.sum(plog[sel]) / _LOG10), n_blocks=int(n_blocks), n_candidates=int(n_cand),
                exact_rank=int(exact), n_pairs_scanned=n_pairs))

    return _collapse(found, lambda n, o: _blocks_for(bits, known, n, o, max_blocks))


def _phase_inhomogeneity_log_p(r: BlockCodeResult, blocks: np.ndarray, period: int) -> float:
    """ln p-value that r's checks hold equally often in every class of block index modulo `period`.

    If r is an echo of a code of length n0 seen through blocks of another length n, the codeword phase inside a
    block cycles with period n0 / gcd(n0, n) and a check holds only in some phases; a genuine code at n gives the
    same satisfaction rate in every class. Chi-square homogeneity test per check, Bonferroni over the checks.
    """
    if period < 2 or blocks.size < 4 * period:
        return 0.0
    checks = (r.H.astype(np.uint64) << np.arange(r.n, dtype=np.uint64)).sum(axis=1).astype(np.uint64)
    cls = np.arange(blocks.size) % period
    totals = np.bincount(cls, minlength=period).astype(float)
    best = 0.0
    for h in checks:
        ok = gf2.parity64(blocks & h) == 0
        sat = np.bincount(cls[ok], minlength=period).astype(float)
        pooled = sat.sum() / totals.sum()
        if pooled <= 0.0 or pooled >= 1.0:
            continue
        expected = totals * pooled
        chi2 = float(np.sum((sat - expected) ** 2 / (expected * (1.0 - pooled))))
        best = min(best, float(stats.chi2.logsf(chi2, period - 1)))
    return best + np.log(max(len(checks), 1))


def _collapse(found: list[BlockCodeResult], blocks_of: Any) -> list[BlockCodeResult]:
    """Fold multiples and echoes of a detected block length into its strongest result.

    Pass 1, by increasing n: a result whose n is a multiple of an already kept smaller n is folded into the
    strongest kept result of that n (field `collapsed`). Pass 2, by decreasing evidence: a result at n is folded
    into a stronger kept result at n0 != n (field `aliases`) when its checks are satisfied at significantly
    different rates in the classes of block index modulo n0 / gcd(n0, n) (chi-square, ln p < ln 1e-6): that is the
    signature of blocks that catch whole codewords of the n0 code only at some phases. Same-n results at other
    offsets are kept.
    """
    stage: list[BlockCodeResult] = []
    for r in sorted(found, key=lambda r: (r.n, -r.evidence)):
        multiple = [b for b in stage if b.n < r.n and r.n % b.n == 0]
        if multiple:
            max(multiple, key=lambda b: b.evidence).collapsed.append((r.n, r.offset))
        else:
            stage.append(r)
    kept: list[BlockCodeResult] = []
    for r in sorted(stage, key=lambda r: (-r.evidence, r.log10_p)):
        host = None
        for base in kept:
            if base.n == r.n:
                continue
            period = base.n // int(np.gcd(base.n, r.n))
            if _phase_inhomogeneity_log_p(r, blocks_of(r.n, r.offset), period) < np.log(1e-6):
                host = base
                break
        if host is not None:
            host.aliases.append((r.n, r.offset))
        else:
            kept.append(r)
    return kept


# ------------------------------------------------------------------------------------------------------- LFSR test

@dataclass
class LfsrResult:
    """Outcome of :func:`lfsr_test`.

    Global part: linear complexity L of n bits, expected_L for random bits, deficit = n/2 - L, p_global (exact
    one-sided tail), and when significant the connection polynomial with regen_mismatches (bits that differ when the
    stream is regenerated from its first L bits). Windowed part: window_L (complexity of every window),
    window_L_min, window_agree (windows sharing the most common polynomial), window_polynomial, p_window
    (surrogates), p_window_analytic (exact bound for window_L_min), window_regen_mismatches. p_value combines both
    parts (Bonferroni over the two).
    """

    L: int
    n: int
    expected_L: float
    deficit: float
    p_global: float
    log10_p_global: float
    polynomial: np.ndarray | None
    regen_mismatches: int | None
    window_len: int
    n_windows: int
    window_L: np.ndarray
    window_L_min: int
    window_agree: int
    window_polynomial: np.ndarray | None
    p_window: float
    p_window_analytic: float
    window_regen_mismatches: int | None
    p_value: float
    method: str


def _log_p_linear_complexity(big_l: int, n: int) -> float:
    """ln P(L(s) <= big_l) for n independent fair bits: (2 * 4^l + 1) / (3 * 2^n) for l <= n/2 (Rueppel 1986)."""
    if 2 * big_l >= n:
        return 0.0
    return float(np.logaddexp(np.log(2.0) + 2 * big_l * np.log(2.0), 0.0) - np.log(3.0) - n * np.log(2.0))


def _regenerate_forward(seed: np.ndarray, poly: np.ndarray, length: int) -> np.ndarray:
    big_l = poly.size - 1
    taps = np.flatnonzero(poly[1:]) + 1
    out = np.zeros(length, dtype=np.uint8)
    out[:big_l] = seed[:big_l]
    for i in range(big_l, length):
        out[i] = np.bitwise_xor.reduce(out[i - taps]) if taps.size else 0
    return out


def _regen_mismatches(bits: np.ndarray, start: int, poly: np.ndarray) -> int:
    """Mismatches when the whole of bits is regenerated from the L bits at `start` (backwards too if c_L = 1)."""
    big_l = poly.size - 1
    if big_l == 0:
        return int(np.count_nonzero(bits))
    fwd = _regenerate_forward(bits[start:], poly, bits.size - start)
    mism = int(np.count_nonzero(fwd != bits[start:]))
    if start > 0:
        if poly[-1] == 1:
            # Run the recurrence backwards: s[i-L] = s[i] + sum_{j<L} c_j s[i-j].
            rev_poly = poly[::-1].copy()
            seq = bits[start: start + big_l][::-1]
            back = _regenerate_forward(seq, rev_poly, start + big_l)[::-1]
            mism += int(np.count_nonzero(back[:start] != bits[:start]))
        else:
            mism += start  # cannot run backwards; count the unreachable prefix as unverified mismatches
    return mism


def _longest_known_run(known: np.ndarray) -> tuple[int, int]:
    if known.all():
        return 0, known.size
    k = np.concatenate([[0], known.astype(np.int8), [0]])
    d = np.diff(k)
    starts, stops = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    if starts.size == 0:
        return 0, 0
    i = int(np.argmax(stops - starts))
    return int(starts[i]), int(stops[i] - starts[i])


def _window_stats(bits: np.ndarray, known: np.ndarray, window: int) -> tuple[np.ndarray, list[tuple[int, int, int]]]:
    """Linear complexity and polynomial (as int) of every fully known non-overlapping window."""
    out_l, polys = [], []
    for start in range(0, bits.size - window + 1, window):
        if not known[start: start + window].all():
            continue
        big_l, c_int = gf2.berlekamp_massey_int(bits[start: start + window])
        out_l.append(big_l)
        polys.append((big_l, c_int, start))
    return np.asarray(out_l, dtype=int), polys


def _window_summary(polys: list[tuple[int, int, int]]) -> tuple[int, int, tuple[int, int, int] | None]:
    """(smallest L, agreement count of the most common polynomial, a representative window of it)."""
    if not polys:
        return 0, 0, None
    counts = Counter((big_l, c) for big_l, c, _ in polys)
    (best_l, best_c), agree = min(counts.items(), key=lambda kv: (-kv[1], kv[0][0]))
    rep = next(p for p in polys if p[0] == best_l and p[1] == best_c)
    return min(p[0] for p in polys), int(agree), rep


def lfsr_test(stream: BitStream, rng: np.random.Generator | None = None, window: int = 256,
              n_surrogates: int = 500, alpha: float = 0.01) -> LfsrResult:
    """Test for a linear-feedback shift register by linear complexity, globally and on short windows.

    Global test. Berlekamp-Massey gives the linear complexity L of the (longest fully known run of the) stream. For
    n independent fair bits P(L <= l) = (2 * 4^l + 1) / (3 * 2^n) for l <= n/2 exactly (Rueppel 1986; mean about
    n/2 + 2/9, variance about 86/81), and that is p_global: H0 is independent fair bits, no search, no correction.
    (For a biased stream the fair-bit formula is only approximate.) When p_global < alpha the connection polynomial
    is returned and checked by regenerating the stream from its first L bits; regen_mismatches counts differences.

    Windowed test. One flipped bit lifts the global L to about n/2, but most short windows of a lightly corrupted
    LFSR stream are clean. Berlekamp-Massey runs on non-overlapping windows of `window` bits; the statistic is the
    number of windows that agree on the same connection polynomial (window_agree), and window_L_min the smallest L.
    p_window is from surrogates: the windowed analysis is rerun on n_surrogates independent-bit streams with the same
    ones-fraction and erasures, p = (1 + #{surrogate agree >= observed}) / (1 + n_surrogates). When window_agree <= 1
    no surrogate can do worse, so p_window = 1 without simulation. p_window_analytic = n_windows * P(L <= L_min) is
    the exact Bonferroni bound for the smallest window complexity. The modal polynomial is verified by regenerating
    the whole stream from one agreeing window (backwards as well when its top coefficient is 1); the mismatch count
    is then roughly the number of flipped bits.

    p_value = min(1, 2 * min(p_global, p_window)): Bonferroni over the two variants. With the default 500
    surrogates a detection by the windowed test alone therefore reports p_value = 2/501 ~ 0.004; raise n_surrogates
    for a smaller floor (the analytic bound p_window_analytic is not limited by it).

    Limits. Windows must exceed 2L, so with the default 256-bit windows LFSRs longer than about 120 cells are
    invisible to the windowed test (and the global one needs an error-free stream of more than 2L bits). Non-linear
    generators, multiplicative (self-synchronising) scramblers applied to non-trivial data, and additive scramblers
    over a payload that is itself random-looking are not detected: L of (LFSR xor payload) is that of the payload.
    """
    bits_all, known_all = _bits_known(stream)
    s0, length = _longest_known_run(known_all)
    bits = bits_all[s0: s0 + length]
    n = int(bits.size)
    big_l, c_int = gf2.berlekamp_massey_int(bits)
    log_pg = _log_p_linear_complexity(big_l, n)
    p_global = _p_from_log(log_pg)
    poly, regen = None, None
    if p_global < alpha:
        poly = gf2.int_to_bits(c_int, big_l + 1)
        regen = _regen_mismatches(bits, 0, poly)

    win_l, polys = _window_stats(bits_all, known_all, window)
    l_min, agree, rep = _window_summary(polys)
    n_windows = int(win_l.size)
    if n_windows:
        p_win_analytic = float(min(1.0, n_windows * np.exp(_log_p_linear_complexity(l_min, window))))
    else:
        p_win_analytic = 1.0
    if agree <= 1 or n_surrogates <= 0:
        p_window = 1.0
        w_method = "window agreement <= 1: minimum attainable, p = 1" if agree <= 1 else "no surrogates"
    else:
        g = _rng(rng, "lfsr_test")
        p1 = _ones_fraction(bits_all, known_all)
        exceed = 0
        for _ in range(n_surrogates):
            sb = _surrogate_bits(bits_all.size, p1, g)
            _, sp = _window_stats(sb, known_all, window)
            if _window_summary(sp)[1] >= agree:
                exceed += 1
        p_window = (1 + exceed) / (1 + n_surrogates)
        w_method = f"window agreement against {n_surrogates} surrogates"
    w_poly, w_regen = None, None
    if rep is not None and agree >= 2:
        w_poly = gf2.int_to_bits(rep[1], rep[0] + 1)
        w_regen = _regen_mismatches(bits_all, rep[2], w_poly)
    p_value = float(min(1.0, 2.0 * min(p_global, p_window)))
    return LfsrResult(L=int(big_l), n=n, expected_L=n / 2 + (4 + n % 2) / 18, deficit=n / 2 - big_l,
                      p_global=p_global, log10_p_global=log_pg / _LOG10, polynomial=poly, regen_mismatches=regen,
                      window_len=window, n_windows=n_windows, window_L=win_l, window_L_min=int(l_min),
                      window_agree=agree, window_polynomial=w_poly, p_window=float(p_window),
                      p_window_analytic=p_win_analytic, window_regen_mismatches=w_regen, p_value=p_value,
                      method=f"global: exact linear-complexity tail; windowed: {w_method}; combined: Bonferroni x2")


# -------------------------------------------------------------------------------------------------------- describe

def _poly_text(poly: np.ndarray | None) -> str | None:
    if poly is None:
        return None
    terms = [("1" if j == 0 else ("x" if j == 1 else f"x^{j}")) for j in np.flatnonzero(poly)]
    return " + ".join(reversed(terms))


def describe(stream: BitStream, rng: np.random.Generator | None = None, alpha: float = 0.01) -> dict[str, Any]:
    """Run all four recognisers and summarise them as a plain dictionary (JSON-serialisable values only).

    frame_scan always runs; sync_word runs on its frame length only when the frame p-value is below alpha;
    block_code_scan reports at most its five strongest results; lfsr_test reports both variants. Each section keeps
    its own search-corrected p-value; no further correction across the four recognisers is applied (multiply the
    smallest by 4 for a Bonferroni statement about the whole stream).
    """
    g = _rng(rng, "describe")
    out: dict[str, Any] = {"n_bits": len(stream), "ones_fraction": _ones_fraction(*_bits_known(stream))}
    try:
        fr = frame_scan(stream, rng=g)
        out["frame"] = {"frame_len": fr.frame_len, "z_score": fr.z_score, "comb_statistic": fr.comb_statistic,
                        "p_value": fr.p_value, "p_bonferroni": fr.p_bonferroni, "method": fr.method}
        if fr.p_value < alpha:
            sw = sync_word(stream, fr.frame_len)
            out["sync"] = {"phase": sw.phase, "length": sw.run_length, "span": sw.run_span,
                           "pattern": "".join(map(str, sw.pattern.tolist())), "pattern_hex": sw.pattern_hex,
                           "p_value": sw.p_value, "log10_p": sw.log10_p}
    except ValueError as err:
        out["frame"] = {"error": str(err)}
    codes = block_code_scan(stream, rng=g)
    out["block_codes"] = [{"n": r.n, "k": r.k, "offset": r.offset, "n_checks": int(r.H.shape[0]),
                           "check_weights": r.weights.tolist(), "ber_estimate": r.ber_estimate,
                           "p_value": r.p_value, "log10_p": r.log10_p, "evidence": r.evidence,
                           "collapsed": [list(c) for c in r.collapsed]} for r in codes[:5]]
    lf = lfsr_test(stream, rng=g)
    out["lfsr"] = {"L": lf.L, "n": lf.n, "deficit": lf.deficit, "p_global": lf.p_global,
                   "polynomial": _poly_text(lf.polynomial), "regen_mismatches": lf.regen_mismatches,
                   "window_L_min": lf.window_L_min, "window_agree": lf.window_agree, "n_windows": lf.n_windows,
                   "window_polynomial": _poly_text(lf.window_polynomial), "p_window": lf.p_window,
                   "p_value": lf.p_value}
    significant = []
    if out.get("frame", {}).get("p_value", 1.0) < alpha:
        significant.append("frame")
    if codes and codes[0].p_value < alpha:
        significant.append("block_code")
    if lf.p_value < alpha:
        significant.append("lfsr")
    out["significant"] = significant
    return out
