# New in Stargazer; replaces the null models scattered through SeaQuest (same author, commit 064d667):
# the within-sequence permutations of seaquest/structure/sequences.py and mi.py and the unigram resampling of
# seaquest/decipher/structure.py (`pareidolia_controls`). Changes: one place for surrogates, a Markov surrogate of
# any order (keeps short-range statistics, destroys long-range structure), block shuffles, and an empirical
# permutation p-value with the +1 correction in place of Gaussian z-scores from 4-20 draws.
"""Surrogate sequences and permutation p-values.

Every function takes an explicit `numpy.random.Generator`; use `stargazer.core.rng.rng_for` to make one.
"""
from __future__ import annotations

from bisect import bisect_right
from typing import Callable

import numpy as np

__all__ = ["block_shuffle", "doublet_shuffle", "make_surrogate", "markov_surrogate", "permutation_pvalue", "shuffle"]


def shuffle(seq: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A uniformly random permutation of `seq`: symbol frequencies kept exactly, all order destroyed."""
    return rng.permutation(np.asarray(seq))


def block_shuffle(seq: np.ndarray, block: int, rng: np.random.Generator) -> np.ndarray:
    """Cut `seq` into consecutive non-overlapping blocks of length `block` and put them in random order.

    Structure shorter than the block survives (except across block joins); longer structure is destroyed. A final
    partial block, if any, is kept as one more (shorter) block and permuted with the others, so every symbol is kept.
    """
    s = np.asarray(seq)
    if block < 1:
        raise ValueError("block must be a positive integer")
    starts = np.arange(0, s.size, block)
    pieces = [s[a:a + block] for a in starts]
    order = rng.permutation(len(pieces))
    return np.concatenate([pieces[i] for i in order]) if pieces else s.copy()


def markov_surrogate(seq: np.ndarray, k: int, order: int, rng: np.random.Generator) -> np.ndarray:
    """A sequence of the same length simulated from the order-`order` transition frequencies of `seq`.

    Transition probabilities are P(x | c) = (n(c, x) + 1/2) / (n(c) + k/2) (add-half smoothing); a context never
    seen in `seq` (reachable only through smoothing) therefore gets the uniform distribution. The first `order`
    symbols are a context copied from a uniformly chosen position of `seq`. Order 0 draws i.i.d. symbols from the
    smoothed symbol frequencies. Keeps the order-`order` statistics approximately (up to estimation noise) and
    destroys longer-range structure.

    The walk is sequential, so it runs in a Python loop (about 0.5 microseconds per symbol); contexts are stored in
    a dictionary keyed by the integer-encoded context, so large alphabets and orders need no dense table.
    """
    s = np.asarray(seq, dtype=np.int64).ravel()
    n = s.size
    if order < 0:
        raise ValueError("order must be non-negative")
    if n == 0:
        return s.copy()
    if order == 0:
        p = np.bincount(s, minlength=k) + 0.5
        return rng.choice(k, size=n, p=p / p.sum())
    if n <= order:
        return s.copy()
    if k ** (order + 1) >= 1 << 62:
        raise ValueError("k**(order + 1) must fit in a 64-bit integer")
    modulus = k ** order
    ctx = np.zeros(n - order + 1, dtype=np.int64)
    for j in range(order):
        ctx = ctx * k + s[j:n - order + 1 + j]
    joint = ctx[:-1] * k + s[order:]
    keys, counts = np.unique(joint, return_counts=True)
    ctx_keys = keys // k
    uniq_ctx, first = np.unique(ctx_keys, return_index=True)
    table = np.full((uniq_ctx.size, k), 0.5)
    table[np.searchsorted(uniq_ctx, ctx_keys), keys % k] += counts
    cdf = np.cumsum(table, axis=1)
    cdf /= cdf[:, -1:]
    cdf[:, -1] = 1.0
    cdfs = [row.tolist() for row in cdf]
    uniform = (np.arange(1, k + 1) / k).tolist()
    uniform[-1] = 1.0
    lookup = dict(zip(uniq_ctx.tolist(), range(uniq_ctx.size)))
    start = int(rng.integers(0, n - order + 1))
    out = s[start:start + order].tolist()
    code = int(ctx[start])
    u = rng.random(n - order).tolist()
    km1 = k - 1
    for x_u in u:
        row = lookup.get(code)
        x = bisect_right(cdfs[row] if row is not None else uniform, x_u)
        if x > km1:
            x = km1
        out.append(x)
        code = (code * k + x) % modulus
    return np.asarray(out, dtype=np.int64)


def doublet_shuffle(seq: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """A uniformly random sequence with exactly the same counts of adjacent pairs, first symbol and last symbol.

    This is the exact first-order shuffle (Altschul & Erickson 1985; Kandel et al. 1996): the sequence is an
    Eulerian path in the multigraph whose edges are the observed transitions, and a uniform random Eulerian path is
    drawn by choosing a random arborescence of "last exits" toward the final symbol (Wilson's algorithm,
    multiplicity-weighted, BEST theorem) and a random order for the other exits of each symbol.

    Why it exists: under any first-order Markov (or i.i.d.) source every sequence with the same transition counts
    and end symbols is equally likely, so the observed sequence and these surrogates are exchangeable and a
    permutation p-value against them is exact. A surrogate simulated from fitted transition frequencies
    (`markov_surrogate`) instead adds a second dose of estimation noise and makes tests conservative.
    """
    s = np.asarray(seq, dtype=np.int64).ravel()
    n = s.size
    if n < 3:
        return s.copy()
    src, dst = s[:-1], s[1:]
    by_src = np.argsort(src, kind="stable")
    sorted_src = src[by_src]
    vertices, starts, degree = np.unique(sorted_src, return_index=True, return_counts=True)
    where = {int(v): j for j, v in enumerate(vertices.tolist())}
    root = int(s[-1])
    last_edge: dict[int, int] = {}
    in_tree = {root}
    nxt_vertex: dict[int, int] = {}
    for v in vertices.tolist():
        u = v
        while u not in in_tree:
            j = where[u]
            e = int(by_src[starts[j] + int(rng.integers(degree[j]))])
            last_edge[u] = e
            nxt_vertex[u] = int(dst[e])
            u = nxt_vertex[u]
        u = v
        while u not in in_tree:
            in_tree.add(u)
            u = nxt_vertex[u]
    key = rng.random(n - 1)
    if last_edge:
        key[np.fromiter(last_edge.values(), dtype=np.int64)] = 2.0
    order = np.lexsort((key, src))
    targets = dst[order].tolist()
    pointer = dict(zip(vertices.tolist(), starts.tolist()))
    out = [int(s[0])]
    u = out[0]
    for _ in range(n - 1):
        p = pointer[u]
        pointer[u] = p + 1
        u = targets[p]
        out.append(u)
    return np.asarray(out, dtype=np.int64)


def make_surrogate(kind: str, k: int) -> Callable[[np.ndarray, np.random.Generator], np.ndarray]:
    """A surrogate generator by name: "shuffle", "doublet" (exact first-order shuffle), "markov<m>" (simulated,
    e.g. "markov1") or "block<L>" (e.g. "block16")."""
    if kind == "shuffle":
        return shuffle
    if kind == "doublet":
        return lambda s, rng: doublet_shuffle(s, k, rng)
    if kind.startswith("markov") and kind[6:].isdigit():
        m = int(kind[6:])
        return lambda s, rng: markov_surrogate(s, k, m, rng)
    if kind.startswith("block") and kind[5:].isdigit():
        b = int(kind[5:])
        return lambda s, rng: block_shuffle(s, b, rng)
    raise ValueError(f"unknown surrogate {kind!r}: use 'shuffle', 'markov<order>' or 'block<length>'")


def permutation_pvalue(statistic: Callable[[np.ndarray], float], seq: np.ndarray,
                       surrogate: Callable[[np.ndarray, np.random.Generator], np.ndarray], n: int,
                       rng: np.random.Generator, *, larger_is_more_extreme: bool = True
                       ) -> tuple[float, np.ndarray, float]:
    """Empirical surrogate test of `statistic(seq)`.

    statistic(seq) -> float; surrogate(seq, rng) -> a new sequence; n, number of surrogates.
    Returns (observed, null, p): the observed value, the n surrogate values, and the p-value
    (1 + #{null at least as extreme as observed}) / (n + 1), which is never zero and is valid (not anti-conservative)
    whenever the observed sequence and its surrogates are exchangeable under the null. NaN surrogate values count
    as not extreme. Use n >= 200 for decisions at the 1 % level (the smallest attainable p is 1 / (n + 1)).
    """
    if n < 1:
        raise ValueError("n must be positive")
    observed = float(statistic(seq))
    null = np.array([float(statistic(surrogate(seq, rng))) for _ in range(n)])
    with np.errstate(invalid="ignore"):
        extreme = null >= observed if larger_is_more_extreme else null <= observed
    p = (1.0 + float(np.count_nonzero(extreme))) / (n + 1.0)
    return observed, null, p
