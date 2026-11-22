# Adapted from SeaQuest (same author, commit 064d667), seaquest/structure/segmentation.py. Changes: the O(N^2)
# batch "Sequitur" loop is replaced by Re-Pair (Larsson & Moffat 2000) on a linked list with a lazy max-heap, near
# linear in the input; the uniform-code MDL gain is replaced by a two-part Krichevsky-Trofimov code length and judged
# against surrogates with an empirical p-value; the transitional-probability-dip chunker is vectorised and gets a
# Markov-1 null (SeaQuest had none and placed about as many boundaries on noise as on signal).
"""Grammar-based structure: Re-Pair compression, its description-length gain against surrogates, and the
transitional-probability-dip chunker with a null distribution."""
from __future__ import annotations

import heapq
import math

import numpy as np
from scipy.special import gammaln

from .estimators import _as_symbols
from .surrogates import make_surrogate

__all__ = ["grammar_gain", "kt_code_bits", "repair", "tp_segmentation", "tp_segmentation_null"]

_LN2 = math.log(2.0)


def kt_code_bits(symbols: np.ndarray) -> float:
    """Krichevsky-Trofimov (adaptive add-half) code length in bits of a symbol string over the m symbols it uses.

    L = -log2[ Gamma(m/2) * prod_s Gamma(n_s + 1/2) / (Gamma(1/2)^m * Gamma(N + m/2)) ], with n_s the count of
    symbol s and N the length. This is an order-0 entropy code that includes the cost of learning the symbol
    frequencies (about (m - 1)/2 * log2 N bits beyond N times the empirical entropy), so codes over large and small
    alphabets can be compared fairly. Which m symbols are used is not charged for.
    """
    s = np.asarray(symbols).ravel()
    if s.size == 0:
        return 0.0
    counts = np.unique(s, return_counts=True)[1].astype(float)
    m = counts.size
    log_p = (gammaln(m / 2.0) + gammaln(counts + 0.5).sum() - m * gammaln(0.5) - gammaln(s.size + m / 2.0))
    return float(-log_p / _LN2)


def _universal_bits(n: int) -> float:
    """Elias-gamma length in bits of the positive integer n + 1 (used for the rule count and sequence length)."""
    return 2.0 * math.floor(math.log2(n + 1)) + 1.0


def repair(seq: np.ndarray, k: int) -> dict:
    """Re-Pair grammar compression.

    Repeatedly replaces the most frequent adjacent pair (counting non-overlapping occurrences, left to right) by a
    new symbol k, k + 1, ... until no pair occurs twice. Ties go to the pair with the smaller (left, right) code,
    so the result is deterministic. Uses a doubly linked list over positions, a set of start positions per pair and
    a lazily updated max-heap, so each replacement costs O(log N); 10**5 symbols take about a second.

    Inputs: seq, integer symbols in [0, k); k, alphabet size.
    Output dict:
      rules: int array (R, 2); row j is the body (left, right) of symbol k + j.
      sequence: the final (top-level) sequence, int array.
      rule_uses: how often each rule symbol appears in the final sequence and in rule bodies.
      description_bits: two-part description length
          DL = KT(sequence) + KT(rule bodies flattened) + gamma(R) + gamma(len(sequence)),
        KT the adaptive code of `kt_code_bits` (each part coded over the symbols it uses) and gamma the Elias-gamma
        length of a count. raw_bits: KT code length of the input; n_rules; input_length.
    """
    s = _as_symbols(seq, k)
    n = int(s.size)
    if n >= 1 << 31 or k + n >= 1 << 31:
        raise ValueError("sequence too long for 32-bit symbol codes")
    sym = s.tolist()
    nxt = list(range(1, n + 1))
    prv = list(range(-1, n - 1))
    if n:
        nxt[-1] = -1
    occ: dict[int, set[int]] = {}
    for i in range(n - 1):
        key = (sym[i] << 32) | sym[i + 1]
        bucket = occ.get(key)
        if bucket is None:
            occ[key] = {i}
        else:
            bucket.add(i)
    heap = [(-len(v), key) for key, v in occ.items() if len(v) >= 2]
    heapq.heapify(heap)
    rules: list[tuple[int, int]] = []
    next_symbol = k
    mask = (1 << 32) - 1

    while heap:
        neg, key = heapq.heappop(heap)
        bucket = occ.get(key)
        actual = len(bucket) if bucket is not None else 0
        if actual != -neg:      # stale entry (counts only fall between pushes): re-queue at the true count
            if actual >= 2:
                heapq.heappush(heap, (-actual, key))
            continue
        a, b = key >> 32, key & mask
        positions = sorted(bucket)
        if a == b:
            used, last_end = 0, -2
            for i in positions:
                if i != last_end:
                    used += 1
                    last_end = nxt[i]
            if used < 2:
                continue        # only overlapping occurrences (a run of three): re-queued if the pair gains more
        del occ[key]
        x = next_symbol
        replaced = 0
        touched: set[int] = set()
        for i in positions:
            if sym[i] != a:
                continue
            j = nxt[i]
            if j < 0 or sym[j] != b:
                continue
            p, q = prv[i], nxt[j]
            if p >= 0:
                old = occ.get((sym[p] << 32) | a)
                if old is not None:
                    old.discard(p)
                    if not old:
                        del occ[(sym[p] << 32) | a]
                new = (sym[p] << 32) | x
                bucket = occ.get(new)
                if bucket is None:
                    occ[new] = {p}
                else:
                    bucket.add(p)
                touched.add(new)
            if q >= 0:
                old = occ.get((b << 32) | sym[q])
                if old is not None:
                    old.discard(j)
                    if not old:
                        del occ[(b << 32) | sym[q]]
                new = (x << 32) | sym[q]
                bucket = occ.get(new)
                if bucket is None:
                    occ[new] = {i}
                else:
                    bucket.add(i)
                touched.add(new)
                prv[q] = i
            sym[i] = x
            sym[j] = -1
            nxt[i] = q
            replaced += 1
        for new in touched:     # only pairs containing the new symbol gained occurrences
            bucket = occ.get(new)
            if bucket is not None and len(bucket) >= 2:
                heapq.heappush(heap, (-len(bucket), new))
        if replaced:
            rules.append((a, b))
            next_symbol += 1
    final = []
    i = 0 if n else -1
    while i >= 0:
        final.append(sym[i])
        i = nxt[i]
    final_arr = np.asarray(final, dtype=np.int64)
    rules_arr = np.asarray(rules, dtype=np.int64).reshape(-1, 2)
    n_rules = rules_arr.shape[0]
    uses = np.zeros(n_rules, dtype=np.int64)
    if n_rules:
        all_syms = np.concatenate([final_arr, rules_arr.ravel()])
        rule_syms = all_syms[all_syms >= k] - k
        uses = np.bincount(rule_syms, minlength=n_rules)
    dl = (kt_code_bits(final_arr) + kt_code_bits(rules_arr.ravel()) + _universal_bits(n_rules)
          + _universal_bits(final_arr.size))
    return {"rules": rules_arr, "sequence": final_arr, "rule_uses": uses, "n_rules": int(n_rules),
            "description_bits": float(dl), "raw_bits": kt_code_bits(s), "input_length": n}


def _grammar_saving(seq: np.ndarray, k: int) -> float:
    g = repair(seq, k)
    return g["raw_bits"] - g["description_bits"]


def grammar_gain(seq: np.ndarray, k: int, rng: np.random.Generator, n_surrogates: int = 200,
                 surrogate: str = "doublet") -> dict:
    """Bits saved by the Re-Pair grammar over the order-0 code of the raw sequence, judged against surrogates.

    gain = raw_bits - description_bits, raw_bits the KT order-0 code length of `seq` (see `kt_code_bits`) and
    description_bits the two-part grammar code of `repair`. The same is computed on `n_surrogates` surrogates and
    p_value = (1 + #{surrogate gain >= observed gain}) / (n_surrogates + 1).

    surrogate: the Markov-1 null, by default realised exactly by "doublet" (`surrogates.doublet_shuffle`: uniform
    over sequences with the same adjacent-pair counts, so the test asks for structure beyond first-order
    transitions and its p-value is exact under any first-order Markov source). "markov1" uses the simulated
    `markov_surrogate` instead; because the transition table is re-estimated and then sampled, its surrogates carry
    extra spurious order-1 structure and the test is conservative (2.3 % false positives at the 5 % level on
    300 i.i.d. sequences of 500 symbols over 4 letters, against 5.0 % with "doublet"). "shuffle", "markov<m>"
    and "block<L>" are also accepted.

    Output dict: gain_bits, gain_per_symbol, raw_bits, description_bits, n_rules, null_gain_bits (array),
    null_mean, p_value, surrogate, n_surrogates.
    """
    s = _as_symbols(seq, k)
    g = repair(s, k)
    observed = g["raw_bits"] - g["description_bits"]
    make = make_surrogate(surrogate, k)
    null = np.array([_grammar_saving(make(s, rng), k) for _ in range(n_surrogates)])
    p = (1.0 + float(np.count_nonzero(null >= observed))) / (n_surrogates + 1.0)
    return {"gain_bits": float(observed), "gain_per_symbol": float(observed / max(s.size, 1)),
            "raw_bits": g["raw_bits"], "description_bits": g["description_bits"], "n_rules": g["n_rules"],
            "null_gain_bits": null, "null_mean": float(null.mean()) if null.size else float("nan"),
            "p_value": p, "surrogate": surrogate, "n_surrogates": int(n_surrogates)}


# ----------------------------------------------------------------------------------------------------------------
# transitional-probability dips
# ----------------------------------------------------------------------------------------------------------------

def _pair_tp(s: np.ndarray, k: int) -> np.ndarray:
    """Bidirectional transitional probability of each adjacent pair, c(ab) / max(c(a), c(b))."""
    uni = np.bincount(s, minlength=k)
    codes = s[:-1] * k + s[1:]
    keys, inverse, counts = np.unique(codes, return_inverse=True, return_counts=True)
    pair_count = counts[inverse.ravel()]
    return pair_count / np.maximum(uni[s[:-1]], uni[s[1:]])


def tp_segmentation(seq: np.ndarray, k: int, threshold: float | None = None) -> np.ndarray:
    """Chunk boundaries at dips of the bidirectional transitional probability (SeaQuest's TP-dip chunker).

    TP(a, b) = min(c(ab)/c(a), c(ab)/c(b)) = c(ab) / max(c(a), c(b)), with unigram counts c(.) over the whole
    sequence. A boundary is placed between positions i and i + 1 when the TP of that pair is a strict local
    minimum (lower than both neighbouring pairs; the sequence ends count as +infinity) or below `threshold`.
    Output: sorted int array of boundary positions including 0 and len(seq), so chunks are seq[b[j]:b[j+1]].

    On its own the boundary set says little: it places dips on i.i.d. noise too. Judge it with
    `tp_segmentation_null`.
    """
    s = _as_symbols(seq, k)
    n = s.size
    if n < 2:
        return np.array([0, n] if n else [0], dtype=np.int64)
    tp = _pair_tp(s, k)
    padded = np.concatenate([[np.inf], tp, [np.inf]])
    dip = (tp < padded[:-2]) & (tp < padded[2:])
    if threshold is not None:
        dip |= tp < threshold
    inner = np.nonzero(dip)[0] + 1
    return np.concatenate([[0], inner, [n]]).astype(np.int64)


def _chunk_stats(s: np.ndarray, bounds: np.ndarray) -> tuple[int, float]:
    """(number of interior boundaries, mean chunk frequency = chunk tokens / distinct chunk types)."""
    data = s.astype(np.int64)
    chunks = {data[a:b].tobytes() for a, b in zip(bounds[:-1].tolist(), bounds[1:].tolist())}
    n_tokens = bounds.size - 1
    return int(bounds.size - 2), float(n_tokens / max(len(chunks), 1))


def tp_segmentation_null(seq: np.ndarray, k: int, rng: np.random.Generator, n: int = 200,
                         threshold: float | None = None, surrogate: str = "doublet") -> dict:
    """The TP-dip chunker on `seq` and on `n` Markov-1 surrogates.

    The default surrogate ("doublet", the exact first-order shuffle) keeps every pair count, hence the whole
    transitional-probability table the chunker is built from, and only rearranges the sequence; so the test asks
    whether the chunks reflect anything beyond first-order transitions. "markov1" (simulated) is also accepted.
    Statistics: the number of interior boundaries and the mean chunk frequency (chunk tokens divided by distinct
    chunk types; recurring "words" raise it).
    Output dict: n_boundaries, mean_chunk_frequency (observed), null_n_boundaries, null_mean_chunk_frequency
    (arrays of length n), p_n_boundaries (two-sided: 2 * min of the two one-sided empirical p-values, capped at
    1), p_mean_chunk_frequency (one-sided, larger is more extreme), boundaries (the observed boundary array).
    """
    s = _as_symbols(seq, k)
    bounds = tp_segmentation(s, k, threshold)
    nb, mcf = _chunk_stats(s, bounds)
    make = make_surrogate(surrogate, k)
    null_nb = np.empty(n)
    null_mcf = np.empty(n)
    for r in range(n):
        sur = make(s, rng)
        null_nb[r], null_mcf[r] = _chunk_stats(sur, tp_segmentation(sur, k, threshold))
    p_hi = (1.0 + np.count_nonzero(null_nb >= nb)) / (n + 1.0)
    p_lo = (1.0 + np.count_nonzero(null_nb <= nb)) / (n + 1.0)
    return {"n_boundaries": nb, "mean_chunk_frequency": mcf, "null_n_boundaries": null_nb,
            "null_mean_chunk_frequency": null_mcf, "p_n_boundaries": float(min(1.0, 2.0 * min(p_hi, p_lo))),
            "p_mean_chunk_frequency": float((1.0 + np.count_nonzero(null_mcf >= mcf)) / (n + 1.0)),
            "boundaries": bounds}
