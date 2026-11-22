"""Anchor search: does an unknown bit stream contain something that sender and receiver could share without
meeting?

An "anchor" is a piece of knowledge that needs no prior contact: counting, the primes, arithmetic, mathematical
constants, and (for human-made signals) the known character codes. No unknown script or signal has ever been read
without one, and an unlimited search for anchors always finds one by chance. This module therefore runs a declared,
fixed, enumerated library of anchors, `LIBRARY`, and every probability it reports is corrected twice: inside each
anchor for every parameter setting that anchor tried (`AnchorHit.p_value`, with the count in `n_internal`), and by
`search` for the number of anchors run (`AnchorReport.p_corrected`, a Bonferroni factor of `n_anchors`).

The six anchors, in library order:

* ``run_length_sequence`` - runs of ones (or zeros) whose lengths follow a declared integer sequence;
* ``field_counting`` - fixed-width binary numbers that count upwards in steps of one;
* ``constants`` - a long exact match with the binary expansion of pi, e, the square root of two or the golden
  ratio;
* ``text_code`` - printable ASCII under 8-bit, 7-bit or asynchronous serial (8N1) readings;
* ``functional_relations`` - records of fixed-width fields in which one field is the sum, product or difference of
  two others, optionally depending on an "operator" field; this anchor also proposes a testable reading;
* ``repeated_block`` - the stream is one block repeated (a beacon; also what a pulsar or a rotating star produces).

Conventions shared by every anchor:

* Erased bits (``known`` False) are read as 0. No anchor models erasures; heavy erasure weakens every statistic.
* "Fraction of ones" is the stream's own fraction; "Markov fit" is an order-1 two-state chain fitted to the stream
  with add-one smoothing. Where an anchor uses several nulls (independent or exchangeable symbols, the Markov fit,
  and for some anchors a specificity or alignment control measured on the stream itself), it uses the largest of
  their chances: the pattern must be surprising under all of them.
* All p-values are upper bounds (union bounds, Bonferroni, binomial tails that dominate the Poisson-binomial), so a
  calibrated stream should give p-values that are, if anything, too large.
* Probabilities are computed in logarithms; ``p_value`` is floored at 1e-300 and the unfloored base-10 logarithm is
  kept in ``evidence["log10_p"]`` (used to order hits that both reach the floor).
* No anchor reads benchmark truth or any benchmark constant. The library is fixed: adding an anchor after looking
  at a stream breaks the correction.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Callable

import numpy as np
from scipy import fft as sfft
from scipy import special

from stargazer.core.rng import rng_for
from stargazer.core.types import BitStream

__all__ = ["AnchorHit", "AnchorReport", "LIBRARY", "search", "run_length_sequence", "field_counting", "constants",
           "text_code", "functional_relations", "repeated_block", "P_FLOOR", "READING_ALPHA"]

P_FLOOR = 1e-300
READING_ALPHA = 1e-3          # within-anchor corrected level at which functional_relations proposes a reading
_LN10 = math.log(10.0)
_NUMBER_WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}


# ============================================================ public containers

@dataclass
class AnchorHit:
    """What one anchor of the library found (or did not find) in one stream.

    anchor: library key, for example "run_length_sequence". family: "mathematics" or "human encoding".
    p_value: probability, under the null stated in the anchor's docstring, of a result at least this strong,
        already corrected for every parameter setting the anchor tried (but not for the number of anchors).
    statistic: the anchor's test statistic (its meaning is in the anchor's docstring).
    parameters: what was found (field width, offset, bit order, which sequence, where).
    sentence: one plain sentence stating only what was observed.
    evidence: the numbers behind the sentence (small arrays allowed).
    n_internal: how many parameter settings the anchor tried (already folded into p_value).
    reading: optional testable interpretation with held-out support (functional_relations only).
    """

    anchor: str
    family: str
    p_value: float
    statistic: float
    parameters: dict[str, Any]
    sentence: str
    evidence: dict[str, Any]
    n_internal: int
    reading: dict[str, Any] | None = None


@dataclass
class AnchorReport:
    """The result of `search`: one hit per anchor run, the declared budget and the corrected p-values.

    hits are ordered most significant first; p_corrected[i] = min(1, hits[i].p_value * n_anchors); best is hits[0]
    when its corrected p-value is at most alpha, else None.
    """

    hits: list[AnchorHit]
    n_anchors: int
    p_corrected: list[float]
    best: AnchorHit | None
    alpha: float


# ============================================================ shared helpers

def _stream_bits(stream: BitStream) -> np.ndarray:
    """The stream's bits as uint8, with erased bits read as 0."""
    bits = np.asarray(stream.bits, dtype=np.uint8)
    if stream.known is not None:
        bits = np.where(stream.known, bits, 0).astype(np.uint8)
    return bits


def _log_binom_sf(k: Any, n: Any, p: Any) -> np.ndarray:
    """Natural log of P(X >= k) for X ~ Binomial(n, p), vectorised; k and n may be non-integer (scaled counts).

    Uses the regularised incomplete beta function, and where that underflows the Chernoff bound
    -n KL(k/n || p), which is an upper bound on the tail and so keeps the result conservative.
    """
    k, n, p = np.broadcast_arrays(np.asarray(k, float), np.asarray(n, float), np.asarray(p, float))
    out = np.zeros(k.shape)
    act = (k > 0) & (n > 0)
    if not act.any():
        return out
    kk, nn, pp = np.minimum(k[act], n[act]), n[act], np.clip(p[act], 0.0, 1.0)
    res = np.zeros(kk.shape)
    res[pp <= 0] = -np.inf
    mid = (pp > 0) & (pp < 1)
    if mid.any():
        km, nm, pm = kk[mid], nn[mid], pp[mid]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
            sf = special.betainc(km, nm - km + 1.0, pm)
            x = km / nm
            t1 = np.where(x > 0, x * np.log(np.where(x > 0, x, 1.0) / pm), 0.0)
            t2 = np.where(x < 1, (1 - x) * np.log(np.where(x < 1, 1 - x, 1.0) / (1 - pm)), 0.0)
            chernoff = np.where(x > pm, -nm * (t1 + t2), 0.0)
            lsf = np.where(sf > 1e-290, np.log(np.maximum(sf, 1e-300)), np.minimum(chernoff, math.log(1e-290)))
        res[mid] = lsf
    out[act] = res
    return out


def _log_binom_sf_scalar(k: float, n: float, p: float) -> float:
    """Scalar version of _log_binom_sf (avoids array overhead in loops)."""
    if k <= 0 or n <= 0:
        return 0.0
    k = min(k, n)
    if p <= 0:
        return -math.inf
    if p >= 1:
        return 0.0
    sf = float(special.betainc(k, n - k + 1.0, p))
    if sf > 1e-290:
        return math.log(sf)
    x = k / n
    t1 = x * math.log(x / p) if x > 0 else 0.0
    t2 = (1 - x) * math.log((1 - x) / (1 - p)) if x < 1 else 0.0
    return min(-n * (t1 + t2), math.log(1e-290)) if x > p else math.log(1e-290)


def _finish(lp: float) -> tuple[float, float]:
    """(p_value floored at P_FLOOR, log10 p unfloored) from a natural-log p-value."""
    lp = min(0.0, float(lp))
    return (max(P_FLOOR, math.exp(lp)) if lp > -745 else P_FLOOR), lp / _LN10


@dataclass(frozen=True)
class _Chain:
    """A two-state order-1 Markov chain: t = (p00, p01, p10, p11), rows sum to one."""

    t: tuple[float, float, float, float]

    @property
    def T(self) -> np.ndarray:
        return np.array(self.t, dtype=float).reshape(2, 2)

    @property
    def pi(self) -> np.ndarray:
        p01, p10 = self.t[1], self.t[2]
        pi1 = p01 / (p01 + p10)
        return np.array([1.0 - pi1, pi1])

    @property
    def logT(self) -> np.ndarray:
        return np.log(self.T)

    @property
    def logpi(self) -> np.ndarray:
        return np.log(self.pi)

    @property
    def eigen(self) -> float:
        """The chain's second eigenvalue, 1 - p01 - p10: the lag-one correlation of its bits."""
        return 1.0 - self.t[1] - self.t[2]


def _fit_chain(bits: np.ndarray) -> _Chain:
    """Order-1 Markov fit with add-one smoothing on every transition count."""
    b = bits.astype(np.int64)
    c = np.bincount(2 * b[:-1] + b[1:], minlength=4).astype(float) + 1.0 if b.size > 1 else np.ones(4)
    t00, t01 = c[0] / (c[0] + c[1]), c[1] / (c[0] + c[1])
    t10, t11 = c[2] / (c[2] + c[3]), c[3] / (c[2] + c[3])
    return _Chain((round(t00, 12), round(t01, 12), round(t10, 12), round(t11, 12)))


@lru_cache(maxsize=64)
def _value_bits(w: int, lsb: bool) -> np.ndarray:
    """The bits of every w-bit value 0..2^w-1 in stream order, shape (2^w, w)."""
    v = np.arange(1 << w, dtype=np.int64)[:, None]
    shifts = np.arange(w) if lsb else np.arange(w - 1, -1, -1)
    return ((v >> shifts[None, :]) & 1).astype(np.int64)


class _FieldModel:
    """Exact probabilities of w-bit fields under a fitted chain, given the bits bordering each field.

    Under an order-1 Markov chain a field is independent of the rest of the stream given the bit just before it and
    the bit just after it, so P(field reads v | rest of stream) = P(v | before, after). Border codes are
    3 * before + after, with 2 standing for "no such bit" (the field touches an end of the stream).
    """

    def __init__(self, chain: _Chain) -> None:
        self.chain = chain
        self._internal: dict[bool, list[np.ndarray]] = {}
        self._edges: dict[int, np.ndarray] = {}

    def _log_internal(self, w: int, lsb: bool) -> np.ndarray:
        """log prod of the w - 1 transitions inside every w-bit value, built up one bit at a time."""
        tables = self._internal.setdefault(lsb, [np.zeros(1), np.zeros(2)])
        logT = self.chain.logT
        while len(tables) <= w:
            k = len(tables)
            v = np.arange(1 << k)
            if lsb:   # stream order is bit 0 first: the new (last) bit is bit k-1
                prev = tables[k - 1][v & ((1 << (k - 1)) - 1)]
                tables.append(prev + logT[(v >> (k - 2)) & 1, v >> (k - 1)])
            else:     # stream order is the top bit first: the new (last) bit is bit 0
                tables.append(tables[k - 1][v >> 1] + logT[(v >> 1) & 1, v & 1])
        return tables[w]

    def _edge_table(self, w: int) -> np.ndarray:
        """log P contributions of the borders, indexed by 4 * code + 2 * first bit + last bit."""
        if w not in self._edges:
            logT, logpi, T = self.chain.logT, self.chain.logpi, self.chain.T
            tw1 = np.linalg.matrix_power(T, w + 1)
            pitw = self.chain.pi @ np.linalg.matrix_power(T, w)
            out = np.zeros(36)
            for a in range(3):
                for b in range(3):
                    if a < 2 and b < 2:
                        norm = math.log(tw1[a, b])
                    elif a == 2 and b < 2:
                        norm = math.log(pitw[b])
                    else:
                        norm = 0.0
                    for first in range(2):
                        for last in range(2):
                            lin = logT[a, first] if a < 2 else logpi[first]
                            lout = logT[last, b] if b < 2 else 0.0
                            out[4 * (3 * a + b) + 2 * first + last] = lin + lout - norm
            self._edges[w] = out
        return self._edges[w]

    def prob(self, values: np.ndarray, codes: np.ndarray, w: int, lsb: bool) -> np.ndarray:
        """P(field reads `values` | its border code), elementwise; 0 for values outside 0..2^w - 1."""
        values = np.asarray(values)
        inr = (values >= 0) & (values < (1 << w))
        v = np.where(inr, values, 0)
        hi, lo = v >> (w - 1), v & 1
        first, last = (lo, hi) if lsb else (hi, lo)
        lp = self._log_internal(w, lsb)[v] + self._edge_table(w)[4 * codes + 2 * first + last]
        return np.where(inr, np.exp(lp), 0.0)


def _windows(bits: np.ndarray, w: int, lsb: bool) -> np.ndarray:
    """The w-bit value starting at every bit position (int64, length n - w + 1)."""
    n = bits.size - w + 1
    out = np.zeros(max(n, 0), dtype=np.int64)
    if n <= 0:
        return out
    for j in range(w):
        out += bits[j:j + n].astype(np.int64) << (j if lsb else w - 1 - j)
    return out


def _boundary_codes(bits: np.ndarray, w: int, off: int, n_f: int) -> np.ndarray:
    """For each field (from bit `off`, width w) the code 3 * (bit before) + (bit after); 2 means none."""
    starts = off + np.arange(n_f, dtype=np.int64) * w
    left = np.full(n_f, 2, dtype=np.int64)
    m = starts >= 1
    left[m] = bits[starts[m] - 1]
    right = np.full(n_f, 2, dtype=np.int64)
    ends = starts + w
    m = ends < bits.size
    right[m] = bits[ends[m]]
    return 3 * left + right


def _longest_true_run(mask: np.ndarray) -> tuple[int, int]:
    """(length, start) of the longest run of True values; (0, -1) if there is none."""
    if mask.size == 0 or not mask.any():
        return 0, -1
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    i = int(np.argmax(e - s))
    return int(e[i] - s[i]), int(s[i])


def _longest_true_run_rows(mask: np.ndarray) -> np.ndarray:
    """The longest run of True in each row of a 2-D boolean array."""
    if mask.shape[1] == 0:
        return np.zeros(mask.shape[0], dtype=np.int64)
    idx = np.arange(1, mask.shape[1] + 1)
    last_false = np.maximum.accumulate(np.where(mask, 0, idx[None, :]), axis=1)
    return (idx[None, :] - last_false).max(axis=1)


def _list_phrase(values: list[int], head: int = 5) -> str:
    if len(values) <= head + 2:
        return ", ".join(str(v) for v in values)
    return ", ".join(str(v) for v in values[:head]) + f", ..., {values[-1]}"


def _order_word(lsb: bool) -> str:
    return "least significant bit first" if lsb else "most significant bit first"


def _hit(anchor: str, family: str, lp: float, statistic: float, parameters: dict, sentence: str, evidence: dict,
         n_internal: int, reading: dict | None = None) -> AnchorHit:
    p, l10 = _finish(lp)
    evidence = dict(evidence)
    evidence["log10_p"] = l10
    return AnchorHit(anchor, family, p, float(statistic), parameters, sentence, evidence, int(n_internal), reading)


# ============================================================ 1. run lengths

SEQUENCES: tuple[str, ...] = ("naturals", "primes", "squares", "powers_of_two", "fibonacci", "triangular")
_SEQUENCE_PHRASE = {"naturals": "consecutive natural numbers", "primes": "consecutive prime numbers",
                    "squares": "consecutive square numbers", "powers_of_two": "consecutive powers of two",
                    "fibonacci": "consecutive Fibonacci numbers", "triangular": "consecutive triangular numbers"}
_MIN_STRETCH = 4
_N_PERMUTATIONS = 1_000


def _sequence_terms(name: str, limit: int) -> np.ndarray:
    """The terms (distinct, increasing, at least 1) of a declared sequence up to `limit`."""
    if limit < 1:
        return np.zeros(0, dtype=np.int64)
    if name == "naturals":
        return np.arange(1, limit + 1, dtype=np.int64)
    if name == "primes":
        sieve = np.ones(limit + 1, dtype=bool)
        sieve[:2] = False
        for p in range(2, int(limit ** 0.5) + 1):
            if sieve[p]:
                sieve[p * p::p] = False
        return np.flatnonzero(sieve).astype(np.int64)
    terms: list[int] = []
    if name == "squares":
        k = 1
        while k * k <= limit:
            terms.append(k * k)
            k += 1
    elif name == "powers_of_two":
        k = 1
        while k <= limit:
            terms.append(k)
            k *= 2
    elif name == "fibonacci":
        a, b = 1, 2
        while a <= limit:
            terms.append(a)
            a, b = b, a + b
    elif name == "triangular":
        k = 1
        while k * (k + 1) // 2 <= limit:
            terms.append(k * (k + 1) // 2)
            k += 1
    else:
        raise ValueError(f"unknown sequence {name!r}")
    return np.array(terms, dtype=np.int64)


def _runs(bits: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(value, length, start) of every maximal run of equal bits."""
    if bits.size == 0:
        z = np.zeros(0, dtype=np.int64)
        return z, z, z
    change = np.flatnonzero(np.diff(bits)) + 1
    starts = np.concatenate([[0], change]).astype(np.int64)
    ends = np.concatenate([change, [bits.size]]).astype(np.int64)
    return bits[starts].astype(np.int64), ends - starts, starts


def run_length_sequence(stream: BitStream, rng: np.random.Generator) -> AnchorHit:
    """Runs of ones (or of zeros) whose lengths follow a declared integer sequence: counting in unary.

    Statistic: L, the longest stretch of consecutive runs (of one polarity) whose lengths are consecutive terms of
    one of the declared sequences (naturals n, n+1, ...; primes; squares; powers of two; Fibonacci numbers 1, 2, 3,
    5, ...; triangular numbers), from any starting term. Twelve settings: two polarities times six sequences.

    Null: run lengths are independent draws from the stream's own empirical run-length distribution q. By the union
    bound the chance of a stretch of length L is at most 12 * (number of start positions) * sum over start terms t
    of prod_j q(a_{t+j}). Stretches shorter than four terms are never reported (p = 1). As secondary evidence only,
    the fraction of 1,000 random permutations of the run lengths that reach a stretch as long is given.

    Cannot see: unary numbers with separators longer than one symbol mixed with the counting symbol, numbers in any
    other notation, sequences outside the declared six, or a stretch broken by a single corrupted run (a bit flip
    splits or merges runs). The first and last runs may be cut short by the ends of the stream.
    """
    bits = _stream_bits(stream)
    rng = rng if rng is not None else rng_for(0, "anchors", "run_length_sequence")
    values, lengths, starts = _runs(bits)
    n_settings = 2 * len(SEQUENCES)
    best: dict[str, Any] | None = None
    for pol in (1, 0):
        sel = values == pol
        lens, st = lengths[sel], starts[sel]
        if lens.size < _MIN_STRETCH:
            continue
        top = int(lens.max())
        q = np.bincount(lens, minlength=top + 1) / lens.size
        with np.errstate(divide="ignore"):
            logq = np.log(q)
        for name in SEQUENCES:
            terms = _sequence_terms(name, top)
            index = np.full(top + 1, -1, dtype=np.int64)
            index[terms] = np.arange(terms.size)
            idx = index[lens]
            step = (idx[1:] == idx[:-1] + 1) & (idx[:-1] >= 0)
            n_steps, s = _longest_true_run(step)
            length = n_steps + 1 if n_steps else int((idx >= 0).any())
            lp = 0.0
            if length >= _MIN_STRETCH:
                lq = logq[terms]
                bad = ~np.isfinite(lq)
                cs = np.concatenate([[0.0], np.cumsum(np.where(bad, 0.0, lq))])
                cb = np.concatenate([[0], np.cumsum(bad)])
                win = cs[length:] - cs[:-length]
                ok = (cb[length:] - cb[:-length]) == 0
                lse = special.logsumexp(win[ok]) if ok.any() else -np.inf
                lp = math.log(n_settings) + math.log(lens.size - length + 1) + float(lse)
            cand = {"lp": min(0.0, lp), "length": length, "polarity": pol, "sequence": name, "idx": idx,
                    "start_run": s if n_steps else -1, "lens": lens, "starts": st, "q": q}
            if best is None or (cand["lp"], -cand["length"]) < (best["lp"], -best["length"]):
                best = cand
    if best is None or best["length"] < _MIN_STRETCH:
        longest = 0 if best is None else best["length"]
        return _hit("run_length_sequence", "mathematics", 0.0, float(longest), {},
                    "No stretch of four or more consecutive runs of ones or zeros has lengths that follow the "
                    "natural numbers, primes, squares, powers of two, Fibonacci or triangular numbers.",
                    {"longest_stretch": longest, "n_runs": int(lengths.size)}, n_settings)
    L, s = best["length"], best["start_run"]
    seq_vals = [int(v) for v in best["lens"][s:s + L]]
    # Secondary evidence: permutations of the run lengths of the winning polarity, scored on the winning sequence.
    idx = best["idx"]
    reach = 0
    for lo in range(0, _N_PERMUTATIONS, 100):
        m = min(100, _N_PERMUTATIONS - lo)
        perm = rng.permuted(np.tile(idx, (m, 1)), axis=1)
        step = (perm[:, 1:] == perm[:, :-1] + 1) & (perm[:, :-1] >= 0)
        reach += int(np.count_nonzero(_longest_true_run_rows(step) + 1 >= L))
    pol_word = "ones" if best["polarity"] == 1 else "zeros"
    start_bit = int(best["starts"][s])
    sentence = (f"{L} consecutive runs of {pol_word}, starting at bit {start_bit:,}, have lengths "
                f"{_list_phrase(seq_vals)}, which are {_SEQUENCE_PHRASE[best['sequence']]}.")
    params = {"polarity": pol_word, "sequence": best["sequence"], "first_term": seq_vals[0], "length": L,
              "start_run": s, "start_bit": start_bit}
    q = best["q"]
    evidence = {"run_lengths": np.array(seq_vals), "n_runs": int(best["lens"].size),
                "q_of_terms": np.array([q[v] for v in seq_vals]), "permutation_fraction": reach / _N_PERMUTATIONS,
                "n_permutations": _N_PERMUTATIONS}
    return _hit("run_length_sequence", "mathematics", best["lp"], L, params, sentence, evidence, n_settings)


# ============================================================ 2. field counting

COUNT_WIDTHS = tuple(range(2, 17))
COUNT_STRIDES = (1, 2, 3)
_COUNT_CONTROLS = (-3, -2, -1, 2, 3)
_MIN_COUNT_STEPS = 3                      # four consecutive numbers, as for run_length_sequence


def _longest_strided_run(hit: np.ndarray, s: int) -> tuple[int, int]:
    """(length, start index) of the longest chain hit[i], hit[i + s], hit[i + 2s], ... that is all True."""
    if hit.size == 0 or not hit.any():
        return 0, -1
    rows = -(-hit.size // s)
    pad = np.zeros(rows * s, dtype=bool)
    pad[:hit.size] = hit
    grid = pad.reshape(rows, s).T                     # row c holds hit[c], hit[c + s], ...
    idx = np.arange(1, rows + 1)
    last_false = np.maximum.accumulate(np.where(grid, 0, idx[None, :]), axis=1)
    run = idx[None, :] - last_false
    c, end = np.unravel_index(int(np.argmax(run)), run.shape)
    length = int(run[c, end])
    return length, int(c + (end - length + 1) * s)


def field_counting(stream: BitStream, rng: np.random.Generator) -> AnchorHit:
    """Fixed-width binary numbers that count upwards: v[i + s] == v[i] + 1, in a long unbroken stretch.

    Grid (810 settings for a long stream): field width w = 2..16 bits, alignment offset 0..w-1, both bit orders,
    and stride s = 1, 2 or 3 fields (a stride of two or three finds numbers separated by constant delimiter
    fields). Statistics: K, the number of positions i with v[i + s] == v[i] + 1 out of n - s; and S, the longest
    unbroken chain of such steps (S steps in a row cover S + 1 numbers counting by one).

    Null chance r of a +1 step, the largest of three: (a) the fields are exchangeable, r = sum_v q(v) q(v+1), q the
    empirical field distribution; (b) the Markov fit, the mean over positions of P(field i + s reads v[i] + 1 | the
    two bits bordering field i + s); (c) specificity: the largest observed rate of a step of -3, -2, -1, +2 or +3
    (structure that merely keeps neighbouring numbers close raises those as much as +1). The p-value of a setting
    is the larger of the binomial tail of K at rate r (conservative for a sum of unequal chances) and the union
    bound n r^S on a chain of S steps; both must be small, so counting must be both frequent and serial. Chains
    shorter than three steps (four numbers) are never reported, as for run_length_sequence. Bonferroni over the
    grid.

    Cannot see: counting by other steps, counting downwards, numbers wider than 16 bits or of variable width, or a
    count so riddled with errors that no long chain survives.
    """
    bits = _stream_bits(stream)
    model = _FieldModel(_fit_chain(bits))
    rows: list[tuple[int, bool, int, int, int, int, float, float, float, int, int]] = []
    for w in COUNT_WIDTHS:
        size = 1 << w
        for lsb in (False, True):
            win = _windows(bits, w, lsb)
            if win.size == 0:
                continue
            for off in range(w):
                v = win[off::w]
                n_f = v.size
                if n_f < 8:
                    continue
                codes = _boundary_codes(bits, w, off, n_f)
                q = np.bincount(v, minlength=size) / n_f
                r_perm = float((q[:-1] * q[1:]).sum())
                for s in COUNT_STRIDES:
                    n_t = n_f - s
                    if n_t < 4:
                        continue
                    dif = v[s:] - v[:-s]
                    hit = dif == 1
                    k = int(np.count_nonzero(hit))
                    if not k:
                        rows.append((w, lsb, off, s, 0, n_t, r_perm, 1.0, 1.0, 0, -1))
                        continue
                    r_m = float(model.prob(v[:-s] + 1, codes[s:], w, lsb).mean())
                    small = np.abs(dif) <= 3
                    near = np.bincount(dif[small] + 3, minlength=7)
                    ctrl = float(max(near[c + 3] for c in _COUNT_CONTROLS)) / n_t
                    run, start = _longest_strided_run(hit, s)
                    rows.append((w, lsb, off, s, k, n_t, r_perm, r_m, ctrl, run, start))
    n_settings = len(rows)
    if not n_settings:
        return _hit("field_counting", "mathematics", 0.0, 0.0, {}, "The stream is too short to cut into fields.",
                    {}, 0)
    arr = np.array([r[4:10] for r in rows], dtype=float)
    k_, n_, r_ = arr[:, 0], arr[:, 1], np.maximum.reduce([arr[:, 2], arr[:, 3], arr[:, 4]])
    lp_k = _log_binom_sf(k_, n_, r_)
    with np.errstate(divide="ignore"):
        lp_s = np.where(arr[:, 5] >= _MIN_COUNT_STEPS, np.log(n_) + arr[:, 5] * np.log(r_), 0.0)
    lp = np.minimum(0.0, np.maximum(lp_k, lp_s))
    # Most significant first; among settings that all reach the reporting floor, the narrower field, then the
    # longer stretch, then the shorter stride (the simplest reading; the p-value is the same floor either way).
    widths = np.array([r[0] for r in rows])
    strides = np.array([r[3] for r in rows])
    floored = np.maximum(lp + math.log(n_settings), math.log(P_FLOOR))
    i = int(np.lexsort((strides, -arr[:, 5], widths, floored))[0])
    w, lsb, off, s, k, n_t, r_perm, r_m, ctrl, run, start = rows[i]
    r_null = max(r_perm, r_m, ctrl)
    lp_corr = float(lp[i]) + math.log(n_settings)
    v = _windows(bits, w, lsb)[off::w]
    stretch = [int(x) for x in v[start:start + (run + 1) * s:s]] if run else []
    params = {"field_bits": w, "bit_order": "lsb_first" if lsb else "msb_first", "offset": off, "stride": s,
              "start_field": start, "start_bit": off + start * w if run else -1}
    evidence = {"k": k, "n": n_t, "r_exchangeable": r_perm, "r_markov": r_m, "r_control_steps": ctrl,
                "expected": n_t * r_null, "longest_stretch_steps": run, "log10_p_count": float(lp_k[i]) / _LN10,
                "log10_p_stretch": float(lp_s[i]) / _LN10,
                "stretch_values": np.array(stretch[:16] + (stretch[-4:] if len(stretch) > 16 else []))}
    where = "the next number" if s == 1 else f"the number {_NUMBER_WORDS.get(s, s)} places on"
    if lp_corr > math.log(0.05):
        sentence = (f"No counting found: the strongest reading ({w}-bit numbers, {_order_word(lsb)}, from bit "
                    f"{off}) has {where} one more than the current one at {k:,} of {n_t:,} places "
                    f"({n_t * r_null:.1f} expected), in unbroken stretches of at most {run} steps.")
    else:
        sentence = (f"Read as {w}-bit numbers ({_order_word(lsb)}, from bit {off}), {where} is one more than the "
                    f"current one at {k:,} of {n_t:,} places, where about {n_t * r_null:.1f} would be expected; "
                    f"the longest unbroken stretch counts from {stretch[0]} to {stretch[-1]} in {run} steps.")
    return _hit("field_counting", "mathematics", lp_corr, k, params, sentence, evidence, n_settings)


# ============================================================ 3. constants

CONSTANT_BITS = 16_384
_GUARD = 64


@lru_cache(maxsize=1)
def _constant_expansions() -> dict[str, tuple[str, np.ndarray]]:
    """For pi, e, sqrt(2) and the golden ratio: (integer part in binary, bits of integer part followed by the first
    CONSTANT_BITS fractional bits). Computed with Python integers and CONSTANT_BITS + 64 guard bits."""
    prec = CONSTANT_BITS + _GUARD
    one = 1 << prec

    def arctan_inv(x: int) -> int:
        total = term = one // x
        x2, k, sign = x * x, 1, -1
        while term:
            term //= x2
            total += sign * (term // (2 * k + 1))
            sign, k = -sign, k + 1
        return total

    pi = 16 * arctan_inv(5) - 4 * arctan_inv(239)          # Machin's formula
    e, term, k = 0, one, 0
    while term:
        e += term
        k += 1
        term //= k
    sqrt2 = math.isqrt(2 << (2 * prec))
    phi = (one + math.isqrt(5 << (2 * prec))) // 2
    out = {}
    for name, x in (("pi", pi), ("e", e), ("sqrt2", sqrt2), ("golden_ratio", phi)):
        integer = x >> prec
        frac = (x >> _GUARD) & ((1 << CONSTANT_BITS) - 1)
        int_bits = format(integer, "b")
        text = int_bits + format(frac, f"0{CONSTANT_BITS}b")
        out[name] = (int_bits, np.frombuffer(text.encode(), dtype=np.uint8) - ord("0"))
    return out


def _window64(bits: np.ndarray) -> np.ndarray:
    """The 64 bits starting at every position as one uint64 (first bit most significant; zero-padded at the
    end)."""
    n = bits.size
    pad = np.concatenate([bits.astype(np.uint64), np.zeros(64, dtype=np.uint64)])
    out = np.zeros(n, dtype=np.uint64)
    for j in range(64):
        out |= pad[j:j + n] << np.uint64(63 - j)
    return out


@lru_cache(maxsize=1)
def _constant_windows() -> dict[str, np.ndarray]:
    return {name: _window64(ref) for name, (_, ref) in _constant_expansions().items()}


@lru_cache(maxsize=512)
def _sorted_kgrams(name: str, k: int) -> np.ndarray:
    """The distinct k-bit windows of one constant's expansion, sorted (cached: the expansions never change)."""
    ref = _constant_expansions()[name][1]
    return np.unique(_constant_windows()[name][:ref.size - k + 1] >> np.uint64(64 - k))


def _longest_match_with(name: str, wa: np.ndarray, na: int) -> tuple[int, int, int]:
    """(length, position in the stream, position in the expansion) of the longest exact common run of bits between
    the stream (given its 64-bit windows) and one constant's expansion. Lengths up to 64 by binary search on exact
    k-bit window values; longer ones by chaining consecutive exact 64-bit matches along one diagonal."""
    ref = _constant_expansions()[name][1]
    wb, nb = _constant_windows()[name], ref.size

    def check(k: int) -> tuple[bool, int, int]:
        if k > na or k > nb:
            return False, -1, -1
        sb = _sorted_kgrams(name, k)
        a = wa[:na - k + 1] >> np.uint64(64 - k)
        pos = np.minimum(np.searchsorted(sb, a), sb.size - 1)
        hit = sb[pos] == a
        if not hit.any():
            return False, -1, -1
        i = int(np.argmax(hit))
        return True, i, int(np.argmax((wb[:nb - k + 1] >> np.uint64(64 - k)) == a[i]))

    ok64, _, _ = check(64)
    if not ok64:
        lo, hi, loc = 0, 63, (0, 0)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            ok, i, j = check(mid)
            if ok:
                lo, loc = mid, (i, j)
            else:
                hi = mid - 1
        return lo, loc[0], loc[1]
    a, b = wa[:na - 63], wb[:nb - 63]
    order = np.argsort(b, kind="stable")
    sb = b[order]
    left, right = np.searchsorted(sb, a, "left"), np.searchsorted(sb, a, "right")
    pairs_i, pairs_j = [], []
    for i in np.flatnonzero(right > left).tolist():
        for t in range(left[i], right[i]):
            pairs_i.append(i)
            pairs_j.append(int(order[t]))
    pi_, pj_ = np.array(pairs_i), np.array(pairs_j)
    d = pi_ - pj_
    o = np.lexsort((pi_, d))
    pi_, d = pi_[o], d[o]
    brk = np.flatnonzero((np.diff(d) != 0) | (np.diff(pi_) != 1)) + 1
    seg_s = np.concatenate([[0], brk])
    seg_e = np.concatenate([brk, [pi_.size]])
    g = int(np.argmax(seg_e - seg_s))
    i0 = int(pi_[seg_s[g]])
    return int(seg_e[g] - seg_s[g]) + 63, i0, i0 - int(d[seg_s[g]])


_CONSTANT_PHRASE = {"pi": "pi", "e": "e", "sqrt2": "the square root of two", "golden_ratio": "the golden ratio"}


def constants(stream: BitStream, rng: np.random.Generator) -> AnchorHit:
    """A long exact match with the binary expansion of pi, e, sqrt(2) or the golden ratio.

    Statistic: L, the longest run of bits the stream shares exactly with the integer part followed by the first
    16,384 fractional bits of any of the four constants, at any alignment in both (so a match may begin in the
    natural reading "11.001001..." for pi, or anywhere in the fraction).

    Null: the stream's bits are independent with its own fraction of ones p1; the chance that a given alignment
    matches L bits is at most max(p1, 1 - p1)^L, so by the union bound p <= n * m * 4 * max(p1, 1 - p1)^L, with n
    and m the numbers of start positions in the stream and in each expansion. n_internal is 4 * n * m alignments.

    Cannot see: the constants in another base, with bits inverted or reversed, beyond the first 16,384 bits, in
    fixed-width digit groups (BCD, ASCII digits), or with any bit error inside the match (a single flip halves the
    longest run).
    """
    bits = _stream_bits(stream)
    n = bits.size
    p1 = float(bits.mean()) if n else 0.5
    pmax = max(p1, 1.0 - p1)
    wa = _window64(bits)
    best = (0, 0, 0, "pi")
    m_total = 0
    for name, (int_bits, ref) in _constant_expansions().items():
        m_total = ref.size
        L, i, j = _longest_match_with(name, wa, n)
        if L > best[0]:
            best = (L, i, j, name)
    L, i, j, name = best
    int_len = len(_constant_expansions()[name][0])
    n_internal = 4 * n * m_total
    if L == 0 or pmax >= 1.0:
        lp = 0.0
    else:
        lp = math.log(4) + math.log(max(n - L + 1, 1)) + math.log(max(m_total - L + 1, 1)) + L * math.log(pmax)
    frac_pos = j - int_len
    if frac_pos < 0:
        where = (f"from the start of its natural reading (integer part first)" if j == 0 else
                 f"from {-frac_pos} bits before the binary point")
    else:
        where = f"from fractional bit {frac_pos + 1:,}"
    params = {"constant": name, "length": L, "stream_position": i, "constant_position": j,
              "fraction_position": frac_pos, "natural_reading_start": j == 0}
    evidence = {"p_one": p1, "match_bits": bits[i:i + min(L, 512)].copy(), "n": n, "m": m_total}
    if lp > math.log(0.05):
        sentence = (f"No constant found: the longest exact run of bits shared with pi, e, the square root of two "
                    f"or the golden ratio is {L} bits (with {_CONSTANT_PHRASE[name]}), no longer than chance "
                    "gives.")
    else:
        sentence = (f"Bits {i:,} to {i + L - 1:,} of the stream ({L:,} bits) match the binary expansion of "
                    f"{_CONSTANT_PHRASE[name]} exactly, {where}.")
    return _hit("constants", "mathematics", lp, L, params, sentence, evidence, n_internal)


# ============================================================ 4. text codes

PRINTABLE = np.zeros(256, dtype=bool)
PRINTABLE[32:127] = True
PRINTABLE[[9, 10, 13]] = True
TEXT_READINGS: tuple[tuple[str, int, bool], ...] = tuple(
    [("8bit", o, lsb) for lsb in (False, True) for o in range(8)]
    + [("7bit", o, lsb) for lsb in (False, True) for o in range(7)]
    + [("8n1", o, True) for o in range(10)])


def _text_null_probabilities(kind: str, lsb: bool, p1: float, chain: _Chain) -> tuple[float, float]:
    """Chance that one character of a reading counts (printable, and well framed for 8N1) under independent bits
    with fraction of ones p1, and under the Markov fit (exact, by summing over every character code)."""
    w = 7 if kind == "7bit" else 8
    pat = _value_bits(w, lsb if kind != "8n1" else True)
    printable = PRINTABLE[:1 << w]
    ones = pat.sum(axis=1)
    p_iid = float(np.sum(np.where(printable, p1 ** ones * (1 - p1) ** (w - ones), 0.0)))
    if kind == "8n1":
        p_iid *= (1 - p1) * p1
        pat = np.concatenate([np.zeros((pat.shape[0], 1), np.int64), pat, np.ones((pat.shape[0], 1), np.int64)], 1)
    logT, logpi = chain.logT, chain.logpi
    lp = logpi[pat[:, 0]] + logT[pat[:, :-1], pat[:, 1:]].sum(axis=1)
    return p_iid, float(np.sum(np.exp(lp[printable])))


def _text_reading(bits: np.ndarray, kind: str, off: int, lsb: bool) -> tuple[np.ndarray, np.ndarray]:
    """(character codes, counted mask) of one reading."""
    if kind == "8n1":
        n = (bits.size - off) // 10
        idx = off + 10 * np.arange(n)
        data = _windows(bits, 8, True)[idx + 1] if n else np.zeros(0, np.int64)
        ok = (bits[idx] == 0) & (bits[idx + 9] == 1) & PRINTABLE[data] if n else np.zeros(0, bool)
        return data, ok
    w = 7 if kind == "7bit" else 8
    chars = _windows(bits, w, lsb)[off::w]
    return chars, PRINTABLE[chars]


def _reading_phrase(kind: str, off: int, lsb: bool) -> str:
    if kind == "8n1":
        return (f"asynchronous serial frames (a 0 start bit, eight data bits least significant first, a 1 stop "
                f"bit) from bit {off}")
    return f"{7 if kind == '7bit' else 8}-bit characters ({_order_word(lsb)}, from bit {off})"


def _text_sliding(bits: np.ndarray, kind: str, lsb: bool) -> float:
    """Fraction of all bit positions (every alignment) at which a character of this kind would count."""
    if kind == "8n1":
        n = bits.size - 9
        if n <= 0:
            return 0.0
        data = _windows(bits, 8, True)[1:n + 1]
        return float(((bits[:n] == 0) & (bits[9:9 + n] == 1) & PRINTABLE[data]).mean())
    win = _windows(bits, 7 if kind == "7bit" else 8, lsb)
    return float(PRINTABLE[win].mean()) if win.size else 0.0


def text_code(stream: BitStream, rng: np.random.Generator) -> AnchorHit:
    """Printable ASCII under a small declared set of character readings.

    Readings (40): 8-bit groups at offsets 0..7 in both bit orders; 7-bit groups at offsets 0..6 in both orders;
    asynchronous serial frames (start bit 0, 8 data bits least significant first, stop bit 1) at offsets 0..9.
    Statistic: K, the number of characters that are printable ASCII (32..126, tab, line feed, carriage return); for
    serial frames a character counts only if its start and stop bits are also right.

    Null chance that one character counts, the largest of three: independent bits with the stream's fraction of
    ones and the Markov fit (each computed exactly by summing over all character codes), and the alignment null:
    the fraction of all bit positions, at every alignment, at which a character of that kind would count. Text has
    an alignment, so its characters must be printable at that alignment more often than at an arbitrary one; a line
    code such as Manchester, whose chips make almost any window printable, does not pass. K is binomial; Bonferroni
    over the 40 readings.

    Cannot see: text under a line code, an error-correcting code or a scrambler (it is hidden by design, though
    some codes leave so much of each character in place that a reading still passes, which the preview shows),
    other character sets (EBCDIC, Baudot, Unicode beyond ASCII), framing with parity or two stop bits, or text
    interleaved with other data. Printability is necessary but not sufficient for text; the evidence shows the
    decoded characters so a reader can judge.
    """
    bits = _stream_bits(stream)
    chain = _fit_chain(bits)
    p1 = float(bits.mean()) if bits.size else 0.5
    nulls: dict[tuple[str, bool], tuple[float, float, float]] = {}
    rows = []
    for kind, off, lsb in TEXT_READINGS:
        key = (kind, lsb)
        if key not in nulls:
            nulls[key] = _text_null_probabilities(kind, lsb, p1, chain) + (_text_sliding(bits, kind, lsb),)
        chars, ok = _text_reading(bits, kind, off, lsb)
        rows.append((kind, off, lsb, int(ok.sum()), int(ok.size)) + nulls[key])
    arr = np.array([r[3:] for r in rows], dtype=float)
    lp = _log_binom_sf(arr[:, 0], arr[:, 1], arr[:, 2:5].max(axis=1))
    frac = np.where(arr[:, 1] > 0, arr[:, 0] / np.maximum(arr[:, 1], 1), 0.0)
    i = int(np.lexsort((-frac, lp))[0])
    kind, off, lsb, k, n, p_iid, p_m, p_slide = rows[i]
    lp_corr = float(lp[i]) + math.log(len(TEXT_READINGS))
    chars, ok = _text_reading(bits, kind, off, lsb)
    shown = np.where(ok & (chars >= 32), chars, ord("."))
    shown = np.where(ok & (chars < 32), ord(" "), shown)
    preview = bytes(shown[:160].astype(np.uint8).tolist()).decode("ascii")
    counted = chars[ok]
    common = np.bincount(counted, minlength=256) if counted.size else np.zeros(256, np.int64)
    space_top = bool(counted.size) and int(np.argmax(common)) == 32
    expected = max(p_iid, p_m, p_slide)
    params = {"reading": kind, "offset": off, "bit_order": "lsb_first" if lsb else "msb_first"}
    evidence = {"preview": preview, "printable_fraction": k / n if n else 0.0, "n_chars": n, "k": k,
                "expected_fraction_iid": p_iid, "expected_fraction_markov": p_m,
                "expected_fraction_any_alignment": p_slide, "space_most_common": space_top}
    space = "; the commonest character is the space" if space_top else ""
    what = "frames are well formed and hold printable ASCII" if kind == "8n1" else "are printable ASCII"
    if lp_corr > math.log(0.05):
        sentence = (f"No text found: the best reading, {_reading_phrase(kind, off, lsb)}, gives {k:,} of {n:,} "
                    f"characters that {'are well framed and printable' if kind == '8n1' else 'are printable'}, "
                    f"about the {100 * expected:.0f}% that this stream gives at any alignment or by chance.")
    else:
        sentence = (f"Read as {_reading_phrase(kind, off, lsb)}, {k:,} of {n:,} "
                    f"{'' if kind == '8n1' else 'characters '}{what}, against about {100 * expected:.0f}% at an "
                    f"arbitrary alignment or by chance{space}.")
    return _hit("text_code", "human encoding", lp_corr, k, params, sentence, evidence, len(TEXT_READINGS))


# ============================================================ 5. functional relations

REL_WIDTHS = (4, 8, 12, 16)
REL_RECORDS = tuple(range(3, 9))
_OPS = ("+", "*", "-")
_OP_WORD = {"+": "plus", "*": "times", "-": "minus"}
_SPAN = max(REL_RECORDS)                  # fields 0..7 relative to a record position
_MIN_RECORDS = 8
_NO_OP = -99
_EXACT_BUDGET = 200
_REL_CONTROLS = (-2, -1, 1, 2)
_RANK_MIN_FRACTION = 0.3                  # stage one ranks only tests holding in at least this share of records


def _build_combos() -> list[tuple[str, int, int, int, int]]:
    """Every relative field triple (dx, dy, dz): distinct, smallest offset 0, span d = 2..7, for each operation;
    for + and * only dx < dy (the two orders are the same test)."""
    combos = []
    for d in range(2, _SPAN):
        for op in _OPS:
            for a in range(d + 1):
                for b in range(d + 1):
                    for c in range(d + 1):
                        distinct = len({a, b, c}) == 3 and min(a, b, c) == 0 and max(a, b, c) == d
                        if distinct and (op == "-" or a < b):
                            combos.append((op, a, b, c, d))
    return combos


_COMBOS = _build_combos()
_COMBO_INDEX = {c[:4]: i for i, c in enumerate(_COMBOS)}
_COMBO_SPAN = np.array([c[4] for c in _COMBOS])


def _build_pairs() -> dict[str, tuple[np.ndarray, ...]]:
    """Per operation: operand pairs (dx, dy) within 0..7, and for each combo its pair index, dz and combo id."""
    out = {}
    for op in _OPS:
        pairs = [(a, b) for a in range(_SPAN) for b in range(_SPAN) if a != b and (op == "-" or a < b)]
        index = {pq: i for i, pq in enumerate(pairs)}
        cp, cz, cid = [], [], []
        for a, b in pairs:
            for c in range(_SPAN):
                key = (op, a, b, c)
                if key in _COMBO_INDEX:
                    cp.append(index[(a, b)])
                    cz.append(c)
                    cid.append(_COMBO_INDEX[key])
        out[op] = (np.array([pq[0] for pq in pairs]), np.array([pq[1] for pq in pairs]), np.array(cp),
                   np.array(cz), np.array(cid))
    return out


_PAIRS = _build_pairs()


def _apply(op: str, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return x + y if op == "+" else (x * y if op == "*" else x - y)


@dataclass
class _RelScan:
    """Everything the relation search counted. A test is identified by (setting, R, residue, combo, e, u): the
    setting (width, bit order, offset), the record length R, the residue modulo R of the record position whose
    field is the triple's offset 0, the combo (operation and relative positions), and, for tests within an operator
    value, the operator field's position e relative to that field and its value u (e = _NO_OP for all records)."""

    bits: np.ndarray
    model: _FieldModel
    settings: list[tuple[int, bool, int]] = field(default_factory=list)
    fields: list[np.ndarray] = field(default_factory=list)
    classes: list[dict[int, list[tuple[str, Any]]]] = field(default_factory=list)
    tests: dict[str, np.ndarray] = field(default_factory=dict)
    n_tests: int = 0
    cache: dict[tuple[int, ...], tuple[float, dict[str, Any]]] = field(default_factory=dict)
    codes: dict[int, np.ndarray] = field(default_factory=dict)

    def border_codes(self, sid: int) -> np.ndarray:
        if sid not in self.codes:
            w, _, off = self.settings[sid]
            self.codes[sid] = _boundary_codes(self.bits, w, off, self.fields[sid].size)
        return self.codes[sid]


def _classify_all(f: np.ndarray, n7: int, w: int) -> dict[int, list[tuple[str, Any]]]:
    """Column roles for every record length R with at least _MIN_RECORDS records, from one count of all columns.

    ("constant", value) if one value fills at least 90% of the column (a delimiter); ("operator", values) if 2..4
    values with at least 5% share each together fill at least 90%; else ("value", None)."""
    Rs = [R for R in REL_RECORDS if n7 // R >= _MIN_RECORDS]
    if not Rs:
        return {}
    keys, col_rows, col_of = [], [], []
    col_id = 0
    for R in Rs:
        rows = n7 // R
        table = f[:rows * R].reshape(rows, R).astype(np.int64)
        keys.append(((col_id + np.arange(R))[None, :] << w) + table)
        col_rows += [rows] * R
        col_of += [(R, c) for c in range(R)]
        col_id += R
    allkeys = np.concatenate([k.ravel() for k in keys])
    if w <= 12:
        counts_all = np.bincount(allkeys, minlength=col_id << w)
        vals = np.flatnonzero(counts_all)
        counts = counts_all[vals]
    else:
        vals, counts = np.unique(allkeys, return_counts=True)
    col = vals >> w
    share = counts / np.array(col_rows)[col]
    value = vals & ((1 << w) - 1)
    starts = np.flatnonzero(np.concatenate([[True], col[1:] != col[:-1]]))
    top = np.maximum.reduceat(share, starts)
    freq = share >= 0.05
    n_freq = np.add.reduceat(freq.astype(np.int64), starts)
    s_freq = np.add.reduceat(np.where(freq, share, 0.0), starts)
    out: dict[int, list[tuple[str, Any]]] = {R: [("value", None)] * R for R in Rs}
    for g, st in enumerate(starts.tolist()):
        end = starts[g + 1] if g + 1 < starts.size else vals.size
        R, c = col_of[int(col[st])]
        if top[g] >= 0.9:
            role: tuple[str, Any] = ("constant", int(value[st:end][np.argmax(share[st:end])]))
        elif 2 <= n_freq[g] <= 4 and s_freq[g] >= 0.9:
            role = ("operator", tuple(int(v) for v in value[st:end][freq[st:end]]))
        else:
            role = ("value", None)
        out[R] = out[R][:c] + [role] + out[R][c + 1:]
    return out


@lru_cache(maxsize=4096)
def _valid_grid(op: str, R: int, vm_bytes: bytes) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For one operation, record length and pattern of value columns: the residue and combo id of every
    (combo, residue) test, flattened combo-major, and whether it is valid (span fits, x, y and z are value
    columns)."""
    vm = np.frombuffer(vm_bytes, dtype=bool)
    PX, PY, CP, CZ, CID = _PAIRS[op]
    rho = np.arange(R)
    span_ok = _COMBO_SPAN[CID] <= R - 1
    valid = (span_ok[:, None] & vm[(rho[None, :] + PX[CP][:, None]) % R] & vm[(rho[None, :] + PY[CP][:, None]) % R]
             & vm[(rho[None, :] + CZ[:, None]) % R])
    return np.tile(rho, CID.size), np.repeat(CID, R), valid.ravel()


def _relations_scan(bits: np.ndarray) -> _RelScan:
    """Stage one of the relation search: for every distinct test, count the informative records (n) and those that
    satisfy the relation (K). Only tests with K >= 1 are kept, but every valid test is counted in n_tests. A cheap
    ranking score (a binomial tail at the chance of the relation for fields drawn from the pooled field
    distribution, estimated on a fixed sample of records) orders the tests for the exact stage two."""
    scan = _RelScan(bits=bits, model=_FieldModel(_fit_chain(bits)))
    blocks: list[tuple[int, int, int, int, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    sample_rng = np.random.default_rng(0)          # fixed: only used for the ranking score, never for p-values

    def keep(sid: int, R: int, rho: np.ndarray, combo: np.ndarray, e: int, u: int, K: np.ndarray, N: np.ndarray,
             r0: float, valid: np.ndarray) -> None:
        sel = valid & (K > 0)
        scan.n_tests += int(np.count_nonzero(valid))
        if sel.any():
            blocks.append((sid, R, e, u, r0, rho[sel], combo[sel], K[sel], N[sel]))

    for w in REL_WIDTHS:
        size = 1 << w
        dtype = np.int16 if w <= 4 else (np.int32 if w <= 12 else np.int64)
        for lsb in (False, True):
            win = _windows(bits, w, lsb).astype(dtype)
            for off in range(w):
                f = win[off::w]
                n7 = f.size - (_SPAN - 1)
                if n7 < _MIN_RECORDS * min(REL_RECORDS):
                    continue
                classes = _classify_all(f, n7, w)
                vmask = {R: np.array([c[0] == "value" for c in cl]) for R, cl in classes.items()}
                vmask = {R: vm for R, vm in vmask.items() if vm.sum() >= 3}
                if not vmask:
                    continue
                sid = len(scan.settings)
                scan.settings.append((w, lsb, off))
                scan.fields.append(f)
                scan.classes.append(classes)
                view = np.lib.stride_tricks.sliding_window_view(f, _SPAN)[:n7]
                q_pool = np.bincount(f, minlength=size) / f.size
                sample = sample_rng.integers(0, n7, size=min(n7, 256))
                for op in _OPS:
                    PX, PY, CP, CZ, CID = _PAIRS[op]
                    X, Y = view[:, PX], view[:, PY]
                    F = _apply(op, X, Y)
                    noninf = (F == X) | (F == Y)
                    Fs, inf_s = F[sample].ravel(), ~noninf[sample].ravel()
                    ok_s = inf_s & (Fs >= 0) & (Fs < size)
                    r0 = float(q_pool[np.where(ok_s, Fs, 0)][ok_s].sum() / max(int(inf_s.sum()), 1))
                    nc = CID.size
                    hits = np.flatnonzero(view[:, CZ] == F[:, CP])
                    pp, cc = np.divmod(hits, nc)
                    good = ~noninf[pp, CP[cc]]
                    pp, cc = pp[good], cc[good]
                    ip, ipair = np.divmod(np.flatnonzero(noninf), PX.size)
                    DXc, DYc = PX[CP], PY[CP]
                    nc = CID.size
                    for R, vm in vmask.items():
                        m = (n7 // R) * R
                        rows = m // R
                        span_ok = _COMBO_SPAN[CID] <= R - 1
                        s1 = pp < m
                        P1, C1 = pp[s1], cc[s1]
                        K = np.bincount(C1 * R + P1 % R, minlength=nc * R)
                        s2 = ip < m
                        I2, J2 = ip[s2], ipair[s2]
                        ninf = np.bincount(J2 * R + I2 % R, minlength=PX.size * R).reshape(PX.size, R)
                        rho_flat, cid_flat, valid = _valid_grid(op, R, vm.tobytes())
                        keep(sid, R, rho_flat, cid_flat, _NO_OP, -1, K, (rows - ninf[CP]).ravel(), r0, valid)
                        ops_cols = [(o, vals) for o, (role, vals) in enumerate(classes[R]) if role == "operator"]
                        if not ops_cols:
                            continue
                        res1, res2 = P1 % R, I2 % R
                        for o, values in ops_cols:
                            for e in range(-(R - 3), R):
                                r_ = (o - e) % R
                                e_ok = span_ok & (e >= -(R - 1 - _COMBO_SPAN[CID]))
                                vx = e_ok & vm[(r_ + DXc) % R] & vm[(r_ + DYc) % R] & vm[(r_ + CZ) % R]
                                if not vx.any():
                                    continue
                                pos = np.arange(r_, m, R)
                                pos = pos[pos + e >= 0]
                                a1 = (res1 == r_) & (P1 + e >= 0)
                                a2 = (res2 == r_) & (I2 + e >= 0)
                                op1, op2 = f[P1[a1] + e], f[I2[a2] + e]
                                c1, j2 = C1[a1], J2[a2]
                                for u in values:
                                    n_u = int(np.count_nonzero(f[pos + e] == u))
                                    k_u = np.bincount(c1[op1 == u], minlength=nc)
                                    ni = np.bincount(j2[op2 == u], minlength=PX.size)
                                    keep(sid, R, np.full(nc, r_), CID, e, u, k_u, n_u - ni[CP], r0, vx)
    sizes = np.array([b[5].size for b in blocks], dtype=np.int64)
    t: dict[str, np.ndarray] = {}
    for j, key in enumerate(("sid", "R", "e", "u", "r0")):
        t[key] = np.repeat(np.array([b[j] for b in blocks]), sizes) if blocks else np.zeros(0)
    for j, key in enumerate(("rho", "combo", "K", "n")):
        t[key] = np.concatenate([b[5 + j] for b in blocks]) if blocks else np.zeros(0, dtype=np.int64)
    t["lp_hat"] = np.zeros(t["K"].size)
    if t["K"].size:
        n = np.maximum(t["n"].astype(float), 1.0)
        rank = t["K"] >= _RANK_MIN_FRACTION * n
        if rank.any():
            t["lp_hat"][rank] = _log_binom_sf(t["K"][rank], n[rank], t["r0"][rank])
    scan.tests = t
    return scan


def _exact(scan: _RelScan, spec: tuple[int, int, int, int, int, int]) -> tuple[float, dict[str, Any]]:
    """Stage two: the exact log p-value of one test, the largest of three nulls, with its numbers.

    The unit is the distinct statement: informative records with identical (x, y, z) are counted once, so a record
    that recurs (a sync word, a repeated word of text, a fixed header) cannot pass for many independent truths.
    Chances per distinct statement: (a) permutation: z shuffled against (x, y), chance q_z(x op y) with q_z the
    distribution of z over the distinct statements; (b) the Markov fit, chance P(field z reads x op y | the two bits
    bordering field z); (c) specificity: the largest observed rate of a near miss z == x op y + c, c in
    {-2, -1, 1, 2}. K is binomial at the mean chance, and p = 1 unless K >= n / 2."""
    if spec in scan.cache:
        return scan.cache[spec]
    sid, R, rho, cid, e, u = spec
    w, lsb, _ = scan.settings[sid]
    f = scan.fields[sid]
    op, dx, dy, dz, _ = _COMBOS[cid]
    m = ((f.size - (_SPAN - 1)) // R) * R
    p = np.arange(rho, m, R)
    if e != _NO_OP:
        p = p[p + e >= 0]
        p = p[f[p + e] == u]
    x, y, z = f[p + dx], f[p + dy], f[p + dz]
    F = _apply(op, x, y)
    inf = (F != x) & (F != y)
    info: dict[str, Any] = {"k_records": int(np.count_nonzero((F == z) & inf)), "n_informative": int(inf.sum()),
                            "n_records": int(p.size)}
    # Distinct statements: identical (x, y, z) triples are counted once.
    trip = np.stack([x[inf], y[inf], z[inf]], axis=1)
    if trip.shape[0]:
        _, first = np.unique(trip, axis=0, return_index=True)
    else:
        first = np.zeros(0, dtype=np.int64)
    Fi, zi, pi_ = F[inf][first], z[inf][first], p[inf][first]
    n = int(first.size)
    k = int(np.count_nonzero(Fi == zi))
    info.update({"k": k, "n_distinct": n, "holds_in_majority": 2 * k >= n})
    if k == 0 or 2 * k < n:
        scan.cache[spec] = (0.0, info)
        return scan.cache[spec]
    q = np.bincount(zi) / zi.size
    okq = (Fi >= 0) & (Fi < q.size)
    r_perm = float(np.where(okq, q[np.clip(Fi, 0, q.size - 1)], 0.0).mean())
    r_m = float(scan.model.prob(Fi, scan.border_codes(sid)[pi_ + dz], w, lsb).mean())
    D = zi - Fi
    r_c = max(int(np.count_nonzero(D == c)) for c in _REL_CONTROLS) / n
    r = max(r_perm, r_m, r_c)
    lp = _log_binom_sf_scalar(k, n, r)
    info.update({"r_permutation": r_perm, "r_markov": r_m, "r_control": r_c, "expected": n * r})
    scan.cache[spec] = (lp, info)
    return scan.cache[spec]


def _spec(scan: _RelScan, t: int) -> tuple[int, int, int, int, int, int]:
    tt = scan.tests
    return tuple(int(tt[k][t]) for k in ("sid", "R", "rho", "combo", "e", "u"))  # type: ignore[return-value]


def _best_test(scan: _RelScan, mask: np.ndarray | None = None) -> tuple[tuple[int, ...] | None, float]:
    """The test with the smallest exact log p among the _EXACT_BUDGET best by the ranking score. Evaluating only a
    subset can miss the best test (a loss of power), never create a false one: the Bonferroni count is over all."""
    lp_hat = scan.tests["lp_hat"]
    if mask is not None:
        lp_hat = np.where(mask, lp_hat, 0.0)
    cand = np.flatnonzero(lp_hat < 0)
    if cand.size > _EXACT_BUDGET:
        cand = cand[np.argpartition(lp_hat[cand], _EXACT_BUDGET)[:_EXACT_BUDGET]]
    best, best_lp = None, 0.0
    for t in cand[np.argsort(lp_hat[cand], kind="stable")].tolist():
        spec = _spec(scan, t)
        lp, _ = _exact(scan, spec)
        if lp < best_lp:
            best, best_lp = spec, lp
    return best, best_lp


def _build_reading(scan: _RelScan) -> dict[str, Any] | None:
    """A testable reading from the most significant relation within an operator value, or None if none reaches the
    READING_ALPHA level (within-anchor corrected). Each frequent value of the operator column is assigned the
    operation (on the same operand and result fields) that is significant at that level, if any."""
    if not scan.tests["K"].size:
        return None
    spec, lp = _best_test(scan, scan.tests["e"] != _NO_OP)
    log_n = math.log(max(scan.n_tests, 1))
    if spec is None or lp + log_n > math.log(READING_ALPHA):
        return None
    sid, R, rho, cid, e, _ = spec
    w, lsb, off = scan.settings[sid]
    op, dx, dy, dz, _ = _COMBOS[cid]
    if op != "-":
        dx, dy = min(dx, dy), max(dx, dy)
    values = scan.classes[sid][R][(rho + e) % R][1]
    operators: dict[int, str] = {}
    for u in values:
        options = []
        for op2 in _OPS:
            key = (op2, min(dx, dy), max(dx, dy), dz) if op2 != "-" else (op2, dx, dy, dz)
            if key in _COMBO_INDEX:
                options.append((_exact(scan, (sid, R, rho, _COMBO_INDEX[key], e, u))[0], op2))
        lp_u, op_u = min(options)
        if lp_u + log_n <= math.log(READING_ALPHA):
            operators[u] = op_u
    if not operators:
        return None
    s_rel = min(0, e)
    s0 = (rho + s_rel) % R
    delimiters = {}
    for c in range(R):
        role, val = scan.classes[sid][R][(s0 + c) % R]
        if role == "constant":
            delimiters[c] = val
    return {"field_bits": w, "bit_order": "lsb_first" if lsb else "msb_first", "bit_offset": off + s0 * w,
            "record_fields": R, "operator_column": e - s_rel, "operators": operators, "delimiters": delimiters,
            "operands": (dx - s_rel, dy - s_rel), "result": dz - s_rel}


def _reading_table(bits: np.ndarray, reading: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    """(records x fields table, start bit of each record) of the complete records a reading defines."""
    w, R, off = reading["field_bits"], reading["record_fields"], reading["bit_offset"]
    f = _windows(bits, w, reading["bit_order"] == "lsb_first")[off::w]
    n_rec = f.size // R
    return f[:n_rec * R].reshape(n_rec, R), off + np.arange(n_rec) * R * w


def _reading_correct(table: np.ndarray, reading: dict[str, Any]) -> np.ndarray:
    """For each record: do its delimiters match, is its operator value known, and is the result right?"""
    ok = np.ones(table.shape[0], dtype=bool)
    for c, v in reading["delimiters"].items():
        ok &= table[:, c] == v
    x, y = table[:, reading["operands"][0]], table[:, reading["operands"][1]]
    z, opv = table[:, reading["result"]], table[:, reading["operator_column"]]
    right = np.zeros(table.shape[0], dtype=bool)
    for u, op in reading["operators"].items():
        right |= (opv == u) & (z == _apply(op, x, y))
    return ok & right


def functional_relations(stream: BitStream, rng: np.random.Generator) -> AnchorHit:
    """Arithmetic stated as records: one field is the sum, product or difference of two others.

    Grid: field width w in {4, 8, 12, 16}, alignment offset 0..w-1, both bit orders, record length R = 3..8 fields,
    and every record phase. Columns of the records x R table are classified from their own value distributions:
    constant (delimiters: one value in at least 90% of records), operator (2..4 values with at least 5% each
    covering at least 90%) and value columns (all others). For every ordered triple (x, y, z) of distinct value
    columns in one record, z == x + y, z == x * y and z == x - y are tested over all records, and separately within
    each frequent value of each operator column of the same record.

    A record is counted only if it is informative: x op y differs from both x and y (so 0 + 0 = 0, x * 1 = x and
    similar coincidences, which repeated or correlated fields satisfy for reasons other than arithmetic, are left
    out), and identical (x, y, z) statements are counted once (a recurring record, such as a sync word, is one
    statement, not many). Statistic: K, the number of distinct informative statements satisfying the relation, out
    of n. A relation is reported only if it is the rule, holding in at least half of the statements tested (K >=
    n/2; otherwise p = 1): records that state arithmetic state true arithmetic, while a deterministic but
    non-arithmetic dependence between fields (a scrambler's linear map, a code's parity) agrees with x op y only
    now and then.

    Null chance per statement, the largest of three: (a) permutation: z shuffled against (x, y), chance q_z(x op y)
    with q_z the empirical distribution of z in the records tested; (b) the Markov fit, chance P(field z reads
    x op y | the two bits bordering field z); (c) specificity: the largest observed rate of a near miss
    z == x op y + c, c in {-2, -1, 1, 2} (dependence that merely makes z close to x op y raises those as well).
    K is binomial at that chance (conservative for unequal chances); Bonferroni over every distinct test, whose
    count is n_internal. A test is identified by width, offset, order, R, the residue modulo R of its first field,
    the relative positions of x, y and z, the operation, and the operator field's relative position and value;
    record phases that place the same fields in one record give the same test and are counted once, as are x + y
    and y + x. The exact p-value is computed for the 200 tests ranked best by a cheap score (record counts against
    the pooled field distribution, among tests holding in at least 30% of records); missing a better test outside
    them loses power but cannot create a false detection.

    Reading: when relations hold within the values of an operator column, the anchor proposes a reading (field
    width, bit order, bit offset of the first record, record length, operator column and the operation each of its
    values stands for, delimiter columns and values, operand and result columns; columns counted from 0). It is
    tested out of sample: the whole search is re-run on the first half of the stream only, and the reading found
    there is scored on the records of the second half (``held_out``: how many it predicts exactly, meaning matching
    delimiters, a known operator value and the right result).

    Cannot see: operations other than +, *, - (division, powers, comparisons), relations across records, fields of
    other widths or of variable width, numbers in other notations (unary, BCD, ASCII digits), or modular arithmetic
    that wraps at the field width.
    """
    bits = _stream_bits(stream)
    scan = _relations_scan(bits)
    n_tests = max(scan.n_tests, 1)
    spec, lp = _best_test(scan) if scan.tests["K"].size else (None, 0.0)
    if spec is None:
        return _hit("functional_relations", "mathematics", 0.0, 0.0, {},
                    "No arithmetic found: no field of any record layout is the sum, product or difference of two "
                    "others more often than chance gives.", {"n_tests": scan.n_tests}, scan.n_tests)
    lp_corr = lp + math.log(n_tests)
    sid, R, rho, cid, e, u = spec
    w, lsb, off = scan.settings[sid]
    op, dx, dy, dz, _ = _COMBOS[cid]
    info = _exact(scan, spec)[1]
    params = {"field_bits": w, "bit_order": "lsb_first" if lsb else "msb_first", "offset": off, "record_fields": R,
              "residue": rho, "x": dx, "y": dy, "z": dz, "operation": op,
              "operator_offset": None if e == _NO_OP else e, "operator_value": None if e == _NO_OP else u}
    evidence: dict[str, Any] = dict(info, n_tests=scan.n_tests)
    reading = _build_reading(scan) if lp_corr <= math.log(READING_ALPHA) else None
    if reading is not None:
        table, starts = _reading_table(bits, reading)
        (xc, yc), zc, oc = reading["operands"], reading["result"], reading["operator_column"]
        parts, per_value = [], {}
        for j, (uv, opu) in enumerate(reading["operators"].items()):
            sub = table[table[:, oc] == uv]
            k_u = int(np.count_nonzero(sub[:, zc] == _apply(opu, sub[:, xc], sub[:, yc])))
            per_value[uv] = {"operation": opu, "satisfied": k_u, "records": int(sub.shape[0])}
            lead = f"field {zc + 1} equals " if j == 0 else ""
            parts.append(f"{lead}field {xc + 1} {_OP_WORD[opu]} field {yc + 1} whenever field {oc + 1} is {uv} "
                         f"({k_u} of {sub.shape[0]} records)")
        half = bits.size // 2
        first = _build_reading(_relations_scan(bits[:half]))
        if first is None:
            held = {"n": int(np.count_nonzero(starts >= half)), "correct": 0, "reading_found_in_first_half": False}
        else:
            t2, s2 = _reading_table(bits, first)
            m2 = s2 >= half
            held = {"n": int(m2.sum()), "correct": int(_reading_correct(t2[m2], first).sum()),
                    "reading_found_in_first_half": True, "same_as_full_reading": first == reading}
        reading = dict(reading, held_out=held)
        tail = "".join(f"; field {c + 1} is always {v}" for c, v in reading["delimiters"].items())
        sentence = (f"In records of {_NUMBER_WORDS.get(R, R)} {w}-bit fields, " + " and ".join(parts) + tail
                    + f" (from bit {reading['bit_offset']}, {_order_word(lsb)}).")
        evidence.update({"per_operator_value": per_value, "first_records": table[:5].copy(),
                         "n_records": int(table.shape[0])})
    elif lp_corr <= math.log(0.05):
        cond = "" if e == _NO_OP else f", whenever the field {e:+d} places from field {dx + 1} is {u},"
        sentence = (f"In records of {_NUMBER_WORDS.get(R, R)} {w}-bit fields ({_order_word(lsb)}, from bit "
                    f"{off}), field {dz + 1} equals field {dx + 1} {_OP_WORD[op]} field {dy + 1}{cond} in "
                    f"{info['k']} of {info['n_distinct']} distinct informative records, where about "
                    f"{info['expected']:.1f} would "
                    f"be expected (fields counted from the first of the three).")
    else:
        sentence = (f"No arithmetic found: the strongest relation tested ({w}-bit fields, records of {R}) holds "
                    f"in {info['k']} of {info['n_distinct']} distinct informative records, about as often as "
                    f"chance gives.")
    return _hit("functional_relations", "mathematics", lp_corr, float(info["k"]), params, sentence, evidence,
                scan.n_tests, reading)


# ============================================================ 6. repeated block

def repeated_block(stream: BitStream, rng: np.random.Generator) -> AnchorHit:
    """The stream is one block repeated: the beacon.

    Statistic: for each lag L in 2..n/2, A(L), the number of positions i with b[i] == b[i + L], out of n - L; the
    most significant lag is reported.

    Null: the Markov fit. A short-memory source agrees with itself at small lags, so the chance of agreement at lag
    L is taken from the fitted chain, pi0^2 + pi1^2 + 2 pi0 pi1 lambda^L (lambda its second eigenvalue). Agreements
    at neighbouring positions are correlated in such a source, so the count is treated as binomial with an
    effective number of trials (n - L) (1 - lambda^2) / (1 + lambda^2). Bonferroni over the lags. The block length
    reported is the smallest divisor of the best lag that agrees nearly as well (within three standard errors).

    A repeating block is what a pulsar, a rotating star or a beacon all produce: repetition alone says nothing
    about origin, and a raster image or framed telemetry also agrees with itself at its line or frame length. This
    anchor cannot see repetition with drifting period, inserted or deleted bits, or a block longer than half the
    stream.
    """
    bits = _stream_bits(stream)
    n = bits.size
    if n < 8:
        return _hit("repeated_block", "mathematics", 0.0, 0.0, {},
                    "The stream is too short to test for repetition.", {}, 0)
    chain = _fit_chain(bits)
    s = bits.astype(float) * 2 - 1
    size = sfft.next_fast_len(2 * n)
    spec = sfft.rfft(s, size)
    corr = sfft.irfft(spec * np.conj(spec), size)[:n]
    lags = np.arange(2, n // 2 + 1)
    trials = (n - lags).astype(float)
    agree = np.rint((trials + corr[lags]) / 2)
    pi0, pi1 = chain.pi
    lam = chain.eigen
    p_null = pi0 ** 2 + pi1 ** 2 + 2 * pi0 * pi1 * lam ** lags.astype(float)
    vif = (1 + lam ** 2) / (1 - lam ** 2) if abs(lam) < 1 else 1e12
    frac = agree / trials
    above = frac > p_null
    lp = np.zeros(lags.size)
    if above.any():
        lp[above] = _log_binom_sf(agree[above] / vif, trials[above] / vif, p_null[above])
    i = int(np.lexsort((-frac, lp))[0])
    lp_corr = float(lp[i]) + math.log(lags.size)
    L = int(lags[i])
    block = L
    se = math.sqrt(max(frac[i] * (1 - frac[i]), 1e-12) / trials[i])
    for dvs in range(2, L // 2 + 1):
        if L % dvs == 0:
            j = dvs - 2
            if frac[j] >= frac[i] - 3 * se:
                block = dvs
                break
    reps = n // block
    view = bits[:reps * block].reshape(reps, block)
    consensus = (view.sum(axis=0) * 2 > reps).astype(np.uint8)
    agree_cons = float((view == consensus).mean())
    j = block - 2
    params = {"block_length": block, "best_lag": L}
    evidence = {"agreement": float(frac[j]), "agreement_null": float(p_null[j]), "repeats": n / block,
                "full_repeats": reps, "agreement_with_consensus": agree_cons, "effective_trial_factor": 1 / vif,
                "consensus_block": consensus.copy() if block <= 4096 else None}
    if lp_corr > math.log(0.05):
        sentence = (f"No repetition found: the stream agrees with itself shifted by {L:,} bits at "
                    f"{100 * frac[i]:.1f}% of positions, as a short-memory model of the stream predicts "
                    f"({100 * p_null[i]:.1f}%).")
    else:
        sentence = (f"The stream agrees with itself shifted by {block:,} bits at {100 * frac[j]:.1f}% of "
                    f"positions, where a short-memory model of the stream predicts {100 * p_null[j]:.1f}%; it "
                    f"reads as a {block:,}-bit block repeated {n / block:.1f} times.")
    return _hit("repeated_block", "mathematics", lp_corr, float(frac[j]), params, sentence, evidence, lags.size)


# ============================================================ the library

LIBRARY: dict[str, Callable[[BitStream, np.random.Generator], AnchorHit]] = {
    "run_length_sequence": run_length_sequence,
    "field_counting": field_counting,
    "constants": constants,
    "text_code": text_code,
    "functional_relations": functional_relations,
    "repeated_block": repeated_block,
}


def search(stream: BitStream, rng: np.random.Generator | None = None, anchors: list[str] | None = None,
           alpha: float = 1e-3) -> AnchorReport:
    """Run the declared anchor library (or the named subset) on a bit stream.

    Every anchor always returns a hit, so the budget is visible. p_corrected multiplies each anchor's p-value by
    the number of anchors run (Bonferroni); `best` is the most significant hit if its corrected p-value is at most
    alpha. Choosing a subset after looking at results defeats the correction: the budget is only honest when
    declared in advance. Each anchor draws from its own child generator, fixed by its place in LIBRARY, so a subset
    run gives the same per-anchor results as the full run.
    """
    names = list(LIBRARY) if anchors is None else list(anchors)
    unknown = [a for a in names if a not in LIBRARY]
    if unknown:
        raise ValueError(f"unknown anchors {unknown}; the library is {list(LIBRARY)}")
    if rng is None:
        rng = rng_for(0, "anchors")
    seeds = rng.integers(0, 2 ** 62, size=len(LIBRARY))
    keys = list(LIBRARY)
    hits = [LIBRARY[a](stream, np.random.default_rng(int(seeds[keys.index(a)]))) for a in names]
    hits.sort(key=lambda h: (h.p_value, h.evidence.get("log10_p", 0.0)))
    n = len(hits)
    p_corr = [min(1.0, h.p_value * n) for h in hits]
    best = hits[0] if hits and p_corr[0] <= alpha else None
    return AnchorReport(hits, n, p_corr, best, alpha)
