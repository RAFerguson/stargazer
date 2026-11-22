"""From light curves to symbols: does a series sit on a few discrete levels, and which level is each point on?

The definition of "discrete" used here is strict, because the loose one is easy to fool. A sinusoid spends most
of its time near its two extremes, and a red-noise curve seen over a short span is lumpy; a two-component
mixture fitted to either says "two levels". Here a series is discrete only if **nearly every point lies within
its own error bar of one of k well-separated levels**: k >= 2 levels, adjacent levels at least `min_separation`
median error bars apart, every level holding real points, and at most `max_off_fraction` of the points further
than `tolerance` error bars from every level.

With the defaults (separation 8, tolerance 3) the +-3 error-bar zones around adjacent levels leave a gap of at
least 2 error bars between them, at least a quarter of the separation, and any smooth distribution of values
puts a corresponding share of its points in the gaps. So continuous curves fail the definition by construction,
whatever the null; the null (`discreteness_test`) is there for the cases where noise and sampling conspire.

What this cannot see: levels closer than `min_separation` error bars (reported as unresolved, not as absent);
states that differ in something other than the mean level (variance, colour, period); a discrete process
observed through a smoothing response slower than the sampling, whose points then sit between levels; and any
meaning of the symbols. A discrete series is not a message: an eclipsing binary is discrete.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import special

from ..core import EventList, Series, SymbolSequence
from ..nulls import drw
from ..trust.decoys import sequential_p

_MAX_BINS = 1024
_MIN_LEVEL_POINTS = 5
_MIN_LEVEL_FRACTION = 0.005


@dataclass
class LevelsResult:
    """The outcome of `find_levels` or `discreteness_test`.

    k: number of levels found; 1 means not discrete.
    levels: the level values, ascending (for k = 1, the inverse-variance weighted mean).
    off_level_fraction: fraction of points further than `tolerance` error bars from every level.
    separation: smallest gap between adjacent levels in units of the median error bar (nan for k = 1).
    occupancy: fraction of points on each level (nearest level).
    p_value: probability that the continuous null of `discreteness_test` looks at least as discrete (by the
        off-level fraction); nan when not tested.
    detected: k >= 2 and, when tested, p_value <= alpha.
    symbols: per-point level index and posterior probabilities, when detected.
    transitions: number of level changes along the series (nearest-level assignment).
    mean_dwell: mean number of consecutive points on one level.
    note: what was found and, when nothing was, why.
    """

    k: int
    levels: np.ndarray
    off_level_fraction: float
    separation: float
    occupancy: np.ndarray
    p_value: float
    detected: bool
    symbols: SymbolSequence | None
    transitions: int
    mean_dwell: float
    note: str


# ---------------------------------------------------------------------------------------------------------------
# Optimal one-dimensional k-means
# ---------------------------------------------------------------------------------------------------------------

def _binned(value: np.ndarray, weight: np.ndarray, width: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sums of weight, weight * value and weight * value^2 in equal bins of the given width (at most _MAX_BINS)."""
    lo, hi = float(value.min()), float(value.max())
    n_bins = int(min(_MAX_BINS, max(1, np.ceil((hi - lo) / width))))
    index = np.minimum(((value - lo) / max(hi - lo, 1e-300) * n_bins).astype(np.int64), n_bins - 1)
    w = np.bincount(index, weight, n_bins)
    s = np.bincount(index, weight * value, n_bins)
    q = np.bincount(index, weight * value * value, n_bins)
    return w, s, q


def _kmeans_levels(value: np.ndarray, weight: np.ndarray, k: int, width: float) -> np.ndarray:
    """Weighted k-means levels in one dimension, optimal on bins of `width` (dynamic programming), then
    polished by a few Lloyd iterations on the unbinned values. Returns k ascending levels (or fewer if the data
    cannot support k distinct clusters)."""
    w, s, q = _binned(value, weight, width)
    keep = w > 0
    w, s, q = w[keep], s[keep], q[keep]
    b = w.size
    k = min(k, b)
    cw = np.concatenate([[0.0], np.cumsum(w)])
    cs = np.concatenate([[0.0], np.cumsum(s)])
    cq = np.concatenate([[0.0], np.cumsum(q)])
    # cost[a, e]: weighted sum of squares of bins a..e-1 as one cluster (a < e)
    a = np.arange(b)[:, None]
    e = np.arange(1, b + 1)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        tw = cw[e] - cw[a]
        cost = (cq[e] - cq[a]) - (cs[e] - cs[a]) ** 2 / tw
    cost = np.where(e > a, np.maximum(cost, 0.0), np.inf)
    best = cost[0].copy()                        # best[e-1]: one cluster over bins 0..e-1
    back = np.zeros((k, b), dtype=np.int64)
    for j in range(1, k):
        # cluster j covers bins a..e-1, the first j clusters cover 0..a-1 (a >= j)
        prev = np.concatenate([[np.inf], best[:-1]])           # prev[a] = best over bins 0..a-1
        total = prev[:, None] + cost
        back[j] = np.argmin(total, axis=0)
        best = total[back[j], np.arange(b)]
    edges = [b]
    for j in range(k - 1, 0, -1):
        edges.append(int(back[j][edges[-1] - 1]))
    edges = edges[::-1]
    starts = [0] + edges[:-1]
    levels = np.array([(cs[e2] - cs[s2]) / (cw[e2] - cw[s2]) for s2, e2 in zip(starts, edges) if e2 > s2])
    for _ in range(5):                            # Lloyd polishing on the raw values
        assign = np.argmin(np.abs(value[:, None] - levels[None, :]), axis=1)
        sums = np.bincount(assign, weight * value, levels.size)
        wts = np.bincount(assign, weight, levels.size)
        levels = np.where(wts > 0, sums / np.where(wts > 0, wts, 1.0), levels)
    return np.unique(levels)


@dataclass
class _Candidate:
    k: int
    levels: np.ndarray
    off: float
    separation: float
    occupancy: np.ndarray
    counts: np.ndarray
    on_counts: np.ndarray


def _evaluate(value: np.ndarray, error: np.ndarray, levels: np.ndarray, tolerance: float) -> _Candidate:
    """Score one level set: refine each level on its own on-level points, then count off-level points."""
    weight = 1.0 / error ** 2
    for _ in range(2):   # re-centre each level on the points that sit on it (ignores points in transit)
        dist = np.abs(value[:, None] - levels[None, :]) / error[:, None]
        near = np.argmin(dist, axis=1)
        on = dist[np.arange(value.size), near] <= tolerance
        sums = np.bincount(near[on], (weight * value)[on], levels.size)
        wts = np.bincount(near[on], weight[on], levels.size)
        levels = np.where(wts > 0, sums / np.where(wts > 0, wts, 1.0), levels)
        order = np.argsort(levels)
        levels = levels[order]
    dist = np.abs(value[:, None] - levels[None, :]) / error[:, None]
    near = np.argmin(dist, axis=1)
    on = dist[np.arange(value.size), near] <= tolerance
    counts = np.bincount(near, minlength=levels.size)
    on_counts = np.bincount(near[on], minlength=levels.size)
    med = float(np.median(error))
    sep = float(np.min(np.diff(levels)) / med) if levels.size > 1 else float("nan")
    return _Candidate(levels.size, levels, float(1.0 - on.mean()), sep, counts / value.size, counts, on_counts)


def _dwell(value: np.ndarray, levels: np.ndarray) -> tuple[int, float]:
    if levels.size < 2 or value.size == 0:
        return 0, float(value.size)
    assign = np.argmin(np.abs(value[:, None] - levels[None, :]), axis=1)
    transitions = int(np.count_nonzero(np.diff(assign)))
    return transitions, value.size / (transitions + 1)


def _search(value: np.ndarray, error: np.ndarray, k_max: int, min_separation: float, tolerance: float,
            max_off_fraction: float, min_points: int, min_fraction: float
            ) -> tuple[_Candidate | None, _Candidate, list[_Candidate]]:
    """(the smallest qualifying candidate or None, the one-level candidate, all candidates for k >= 2)."""
    weight = 1.0 / error ** 2
    width = 0.25 * float(np.median(error))
    one = _evaluate(value, error, np.array([np.sum(weight * value) / np.sum(weight)]), tolerance)
    tried = []
    need = max(min_points, int(np.ceil(min_fraction * value.size)))
    for k in range(2, k_max + 1):
        levels = _kmeans_levels(value, weight, k, width)
        if levels.size < 2:
            break
        cand = _evaluate(value, error, levels, tolerance)
        tried.append(cand)
        if (cand.k == k and cand.separation >= min_separation and cand.off <= max_off_fraction
                and np.all(cand.on_counts >= need)):
            return cand, one, tried
    return None, one, tried


# ---------------------------------------------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------------------------------------------

def find_levels(value: np.ndarray, error: np.ndarray, k_max: int = 6, min_separation: float = 8.0,
                tolerance: float = 3.0, max_off_fraction: float = 0.10, *, min_points: int = _MIN_LEVEL_POINTS,
                min_fraction: float = _MIN_LEVEL_FRACTION) -> LevelsResult:
    """The smallest number of well-separated levels that nearly every point sits on, within its error bar.

    For k = 2..k_max the best level set is found by inverse-variance weighted k-means (optimal by dynamic
    programming on values binned at a quarter of the median error bar, then polished on the raw values), and
    each level is re-centred on the points that lie on it. A solution qualifies when adjacent levels are at least
    `min_separation` median error bars apart, at most `max_off_fraction` of the points lie further than
    `tolerance` of their own error bars from every level, and every level holds at least `min_points` points
    and `min_fraction` of the data on it. The smallest qualifying k is returned; if none qualifies, k = 1.

    The minimum level size: the default is 5 points and half a per cent of the data, not the 1% one might pick,
    because real rare states exist (an eclipse with a duty cycle of a few per cent, a short-lived flare plateau)
    and the size of a level is not what protects against false levels. That protection is the on-level
    requirement: a level is only counted from points that lie within their error bars of it, and stray outliers
    scattered over a range of values do not fall within 3 error bars of a common value. Five points is the
    least that shows a level's scatter is consistent with its error bars.

    No probability is computed (p_value is nan); see `discreteness_test`. The values are taken to be in time
    order for the transition and dwell counts.
    """
    value = np.asarray(value, float)
    error = np.asarray(error, float)
    if value.shape != error.shape or value.ndim != 1:
        raise ValueError("value and error must be one-dimensional arrays of equal length")
    if np.any(error <= 0) or not np.all(np.isfinite(value)) or not np.all(np.isfinite(error)):
        raise ValueError("errors must be positive and all numbers finite")
    nan = float("nan")
    if value.size < 2 * min_points:
        return LevelsResult(1, np.array([float(np.mean(value))]) if value.size else np.zeros(0), nan, nan,
                            np.ones(1), nan, False, None, 0, float(value.size),
                            "too few points to look for levels")
    chosen, one, tried = _search(value, error, k_max, min_separation, tolerance, max_off_fraction, min_points,
                                 min_fraction)
    if chosen is not None:
        transitions, dwell = _dwell(value, chosen.levels)
        note = (f"{chosen.k} levels, {chosen.separation:.1f} error bars apart at the closest, with "
                f"{100 * chosen.off:.1f}% of points off every level")
        return LevelsResult(chosen.k, chosen.levels, chosen.off, chosen.separation, chosen.occupancy, nan, True,
                            None, transitions, dwell, note)
    note = _why_not(one, tried, min_separation, max_off_fraction, min_points, min_fraction, value.size)
    return LevelsResult(1, one.levels, one.off, nan, np.ones(1), nan, False, None, 0, float(value.size), note)


def _why_not(one: _Candidate, tried: list[_Candidate], min_separation: float, max_off_fraction: float,
             min_points: int, min_fraction: float, n: int) -> str:
    """A plain-language reason why no level set qualified."""
    need = max(min_points, int(np.ceil(min_fraction * n)))
    unresolved = [c for c in tried if c.off <= max_off_fraction and c.separation < min_separation]
    if unresolved and one.off > max_off_fraction:
        c = unresolved[0]
        if c.k == 2:
            return (f"not discrete by this definition: two levels would hold the points within their error bars, "
                    f"but they are only {c.separation:.1f} error bars apart (fewer than {min_separation:g}); the "
                    f"levels are unresolved (possibly discrete, but not distinguishable from a continuous spread)")
        return (f"continuous or unresolved: covering the values needs {c.k} levels only {c.separation:.1f} error "
                f"bars apart (fewer than {min_separation:g}), which any spread of values allows")
    if one.off <= max_off_fraction:
        return ("one level: the points are consistent with a constant within their error bars "
                f"({100 * one.off:.1f}% off)")
    small = [c for c in tried if c.off <= max_off_fraction and c.separation >= min_separation
             and np.any(c.on_counts < need)]
    if small:
        return (f"{small[0].k} separated levels fit, but at least one holds fewer than {need} points on it; "
                "not counted as a level")
    best = min(tried, key=lambda c: c.off) if tried else one
    return (f"continuous: even the best level set ({best.k} levels) leaves {100 * best.off:.0f}% of the points "
            f"more than the tolerance from every level")


def _statistic(value: np.ndarray, error: np.ndarray, kw: dict) -> float:
    """The discreteness statistic: off-level fraction of the best qualifying level set, 1.0 when none."""
    chosen, _, _ = _search(value, error, kw.get("k_max", 6), kw.get("min_separation", 8.0),
                           kw.get("tolerance", 3.0), kw.get("max_off_fraction", 0.10),
                           kw.get("min_points", _MIN_LEVEL_POINTS), kw.get("min_fraction", _MIN_LEVEL_FRACTION))
    return 1.0 if chosen is None else chosen.off


def discreteness_test(series: Series, rng: np.random.Generator, *, null: str = "drw", n_decoys: int = 200,
                      alpha: float = 1e-2, **find_levels_kwargs) -> LevelsResult:
    """`find_levels` with a probability: how often does a continuous null look at least this discrete?

    Statistic: the off-level fraction of the smallest qualifying level set, or 1.0 when no level set qualifies.
    Smaller is more discrete.

    Nulls:
    - "drw": decoys are damped random walks drawn from the posterior of a fit to this series (`drw.fit`, then
      `DrwFit.draw` for (sigma, tau) and `drw.simulate` at the same times, with the same errors, around the
      fitted mean). This asks whether red noise, with this sampling, could sit on levels this well.
    - "smooth": a smoothed bootstrap of the values: each decoy value is a randomly chosen observed value plus
      Gaussian noise of three median error bars, keeping its own error bar. For series whose time order is
      irrelevant; it asks whether the values could come from a distribution smooth on the scale of a few
      error bars.

    Each decoy goes through the identical level search. p_value is the fraction of decoys whose statistic is
    at most the observed one (ties count against the series), by Besag and Clifford's sequential scheme; it can
    be no smaller than 1 / (n_decoys + 1). When the series itself has no qualifying level set its statistic is
    already the largest possible, p = 1 exactly, and no decoys are drawn.

    detected requires both a qualifying level set and p_value <= alpha. When detected, `symbols` holds the
    per-point posterior over levels from `symbolise`.
    """
    result = find_levels(series.value, series.error, **find_levels_kwargs)
    if result.k < 2:
        result.p_value = 1.0
        result.note += "; p = 1 (no level set qualifies, so no decoy can be less discrete)"
        return result
    observed = result.off_level_fraction
    time, error = series.time, series.error
    if null == "drw":
        fitted = drw.fit(series)

        def draw(m: int) -> np.ndarray:
            sigma, tau = fitted.draw(m, rng)
            sims = drw.simulate(time, error, sigma, tau, rng, mean=fitted.mean_map)
            return -np.array([_statistic(row, error, find_levels_kwargs) for row in sims])
    elif null == "smooth":
        bandwidth = 3.0 * float(np.median(error))
        value = series.value

        def draw(m: int) -> np.ndarray:
            out = np.empty(m)
            for i in range(m):
                pick = rng.integers(0, value.size, value.size)
                out[i] = -_statistic(value[pick] + bandwidth * rng.standard_normal(value.size), error[pick],
                                     find_levels_kwargs)
            return out
    else:
        raise ValueError("null must be 'drw' or 'smooth'")
    mc = sequential_p(-observed, draw, batch=min(25, n_decoys), max_draws=n_decoys)
    result.p_value = float(mc.p)
    result.detected = bool(mc.p <= alpha)
    if result.detected:
        result.symbols = symbolise(series, result.levels)
        result.note += f"; p = {mc.p:.3g} against the {null} null ({mc.draws} decoys)"
    else:
        result.note += f"; but the {null} null looks as discrete in {100 * mc.p:.0f}% of decoys: not detected"
    return result


def _occupancy_em(value: np.ndarray, error: np.ndarray, levels: np.ndarray, iterations: int = 20) -> np.ndarray:
    weights = np.full(levels.size, 1.0 / levels.size)
    log_like = -0.5 * ((value[:, None] - levels[None, :]) / error[:, None]) ** 2
    for _ in range(iterations):
        post = _posterior(log_like, weights)
        weights = post.mean(axis=0).clip(1e-12)
        weights /= weights.sum()
    return weights


def _posterior(log_like: np.ndarray, weights: np.ndarray) -> np.ndarray:
    log_post = log_like + np.log(weights)[None, :]
    log_post -= log_post.max(axis=1, keepdims=True)
    post = np.exp(log_post)
    return post / post.sum(axis=1, keepdims=True)


def symbolise(series: Series, levels: np.ndarray) -> SymbolSequence:
    """Each point's level, with its posterior probability over all levels.

    The probability of level j for point i is proportional to occupancy_j * N(value_i; level_j, error_i^2): the
    Gaussian error bar of the point times the fraction of time spent on each level, the occupancies estimated
    from the whole series by a few EM iterations with the levels held fixed. Symbol 0 is the lowest level (in
    magnitudes that is the brightest). A point between levels gets a split probability, so statistics computed
    from `sample()` draws carry the uncertainty. meta records the levels, occupancies, and `off_level`: how many
    points lie more than 3 error bars from every level (their symbols are guesses).
    """
    levels = np.sort(np.asarray(levels, float))
    value, error = series.value, series.error
    occupancy = _occupancy_em(value, error, levels)
    log_like = -0.5 * ((value[:, None] - levels[None, :]) / error[:, None]) ** 2
    post = _posterior(log_like, occupancy)
    off = int(np.sum(np.min(np.abs(value[:, None] - levels[None, :]) / error[:, None], axis=1) > 3.0))
    meta = {"levels": levels.tolist(), "occupancy": occupancy.tolist(), "off_level": off, "kind": series.kind,
            "source_id": series.source_id, "method": "levels"}
    return SymbolSequence(np.argmax(post, axis=1), levels.size, post, series.time.copy(), meta)


def quantise(series: Series, k: int = 4, method: str = "quantile") -> SymbolSequence:
    """A plain k-bin symboliser for continuous series. This is a lossy choice imposed by the analyst, not a
    property of the data: the symbols exist because the bins were drawn, and any structure found in them is
    structure of the continuous curve (for example its red noise), not evidence of discrete states.

    method "quantile" puts equal numbers of points in each bin; "uniform" uses equal-width bins between the
    smallest and largest value. The per-point probabilities are the Gaussian error bar's mass in each bin, so a
    point near an edge is uncertain. meta records the edges and that the symbolisation is lossy.
    """
    if k < 2:
        raise ValueError("k must be at least 2")
    value, error = series.value, series.error
    if method == "quantile":
        inner = np.quantile(value, np.arange(1, k) / k)
    elif method == "uniform":
        inner = np.linspace(value.min(), value.max(), k + 1)[1:-1]
    else:
        raise ValueError("method must be 'quantile' or 'uniform'")
    edges = np.concatenate([[-np.inf], inner, [np.inf]])
    cdf = special.ndtr((edges[None, :] - value[:, None]) / error[:, None])
    prob = np.diff(cdf, axis=1).clip(0.0)
    prob /= prob.sum(axis=1, keepdims=True)
    symbols = np.clip(np.searchsorted(inner, value, side="right"), 0, k - 1)
    meta = {"method": method, "edges": inner.tolist(), "lossy": True,
            "warning": "bins imposed on a continuous series; not evidence of discrete states",
            "source_id": series.source_id}
    return SymbolSequence(symbols, k, prob, series.time.copy(), meta)


def event_symbols(events: EventList, k_max: int = 6, *, min_separation: float = 8.0, tolerance: float = 3.0,
                  max_off_fraction: float = 0.10) -> SymbolSequence | None:
    """The level search on the logarithms of the intervals between events; None if they are not discrete.

    Intervals carry no error bars, so for each k the scatter is estimated robustly from the data: 1.4826 times
    the median absolute deviation of the log-intervals from their nearest k-means level, pooled over levels.
    This is a self-referential error bar, so the separation criterion does the work: a smooth distribution cut
    into k clusters has levels only about 3.5 scatters apart, far short of 8, whereas genuinely discrete
    intervals (pulse-position coding with small timing jitter) have levels many scatters apart. No probability
    is attached; the spacing test in `stargazer.detect.intervals` calibrates the related question of a common
    quantum. Symbol j is the j-th shortest interval level; time is the time of the event ending each interval;
    meta records the levels in time units and the scatter in log units.
    """
    d = np.diff(np.sort(np.asarray(events.time, float)))
    keep = d > 0
    log_d = np.log(d[keep])
    if log_d.size < 2 * _MIN_LEVEL_POINTS:
        return None
    ones = np.ones_like(log_d)
    for k in range(2, k_max + 1):
        levels = _kmeans_levels(log_d, ones, k, max(float(np.ptp(log_d)) / _MAX_BINS, 1e-12))
        if levels.size < k:
            break
        near = np.argmin(np.abs(log_d[:, None] - levels[None, :]), axis=1)
        scatter = max(1.4826 * float(np.median(np.abs(log_d - levels[near]))), 1e-9)
        error = np.full(log_d.size, scatter)
        cand = _evaluate(log_d, error, levels, tolerance)
        need = max(_MIN_LEVEL_POINTS, int(np.ceil(_MIN_LEVEL_FRACTION * log_d.size)))
        if (cand.k == k and cand.separation >= min_separation and cand.off <= max_off_fraction
                and np.all(cand.on_counts >= need)):
            log_like = -0.5 * ((log_d[:, None] - cand.levels[None, :]) / scatter) ** 2
            post = _posterior(log_like, _occupancy_em(log_d, error, cand.levels))
            meta = {"levels": np.exp(cand.levels).tolist(), "log_scatter": scatter, "off_level_fraction": cand.off,
                    "separation": cand.separation, "method": "log_interval_levels", "dropped_zero_intervals":
                    int(np.sum(~keep)), "source_id": events.source_id}
            return SymbolSequence(np.argmax(post, axis=1), k, post, np.sort(events.time)[1:][keep], meta)
    return None
