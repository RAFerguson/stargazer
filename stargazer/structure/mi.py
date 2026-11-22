# Adapted from SeaQuest (same author, commit 064d667), seaquest/structure/mi.py. Changes: pairs are counted with
# vectorised integer codes instead of Python lists of tuples; the shuffle null defaults to 200 draws and gives a
# 95 % band and per-lag empirical p-values (SeaQuest used 10-20 draws and no p-values); lags can be time-separation
# bins on an irregular time axis; `fit_decay` fits all three models by least squares on the mutual information
# itself and compares them by AIC on those residuals (SeaQuest compared log-space AICs with a linear-space one).
"""Mutual-information decay with lag (Sainburg et al. 2019; Lin & Tegmark 2017) against a shuffle null, and
exponential, power-law and composite fits to the decay. Exponential decay is what a finite-state Markov source
gives; a power-law tail is weak evidence of hierarchy, especially on short sequences."""
from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from .estimators import _as_symbols, _mi_curve, lag_pair_index, shuffled_mi_null

__all__ = ["fit_decay", "mutual_information_decay"]

HIERARCHY_CAVEAT = ("a power-law tail is weak evidence of hierarchy: on short sequences the mutual information at "
                    "long lags is at the noise floor, where the models are hard to tell apart, and mixtures of "
                    "Markov processes can mimic a power law over a limited range")


def mutual_information_decay(seq: np.ndarray, k: int, max_lag: int, n_shuffle: int, rng: np.random.Generator, *,
                             time: np.ndarray | None = None, lag_edges: np.ndarray | None = None) -> dict:
    """Shuffle-corrected mutual information (bits) as a function of lag, with a null band and p-values.

    Index mode: lags 1..max_lag, pairs (seq[i], seq[i + lag]). Time mode (`time` and `lag_edges` given; max_lag is
    ignored): one value per time-separation bin [lag_edges[b], lag_edges[b + 1]), in the units of `time`, over all
    pairs whose separation falls in the bin (see estimators.lagged_mutual_information).

    The null is n_shuffle random permutations of the symbols over the same positions (or times); use >= 200.
    Output dict:
      lags: integer lags, or bin centres in time mode; lag_edges (time mode, else None);
      mi_plugin: plug-in MI; mi: mi_plugin - shuffle_mean (bias-corrected; can be negative);
      shuffle_mean, shuffle_low, shuffle_high: mean and 2.5 / 97.5 percentiles of the plug-in MI under shuffling;
      p_value: per lag, (1 + #{shuffled >= observed}) / (n_shuffle + 1), not corrected for the number of lags;
      n_pairs: pairs per lag; n_shuffle.
    """
    s = _as_symbols(seq, k)
    if time is None and lag_edges is None:
        lags = np.arange(1, int(max_lag) + 1)
        pairings, order = lag_pair_index(s.size, lags)
        centres = lags
        edges = None
    else:
        pairings, order = lag_pair_index(s.size, None, time=time, lag_edges=lag_edges)
        edges = np.asarray(lag_edges, dtype=float)
        centres = 0.5 * (edges[:-1] + edges[1:])
    s = s[order]
    observed, _, n_pairs = _mi_curve(s, k, pairings)
    null = shuffled_mi_null(s, k, pairings, n_shuffle, rng)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(null, axis=0)
        low, high = np.nanpercentile(null, [2.5, 97.5], axis=0)
    p = (1.0 + (null >= observed).sum(axis=0)) / (n_shuffle + 1.0)
    p = np.where(np.isfinite(observed), p, np.nan)
    return {"lags": centres, "lag_edges": edges, "mi_plugin": observed, "mi": observed - mean, "shuffle_mean": mean,
            "shuffle_low": low, "shuffle_high": high, "p_value": p, "n_pairs": n_pairs, "n_shuffle": int(n_shuffle)}


def _exp(x: np.ndarray, a: float, b: float) -> np.ndarray:
    return a * np.exp(-b * x)


def _power(x: np.ndarray, a: float, e: float) -> np.ndarray:
    return a * x ** (-e)


def _composite(x: np.ndarray, a: float, b: float, c: float, e: float) -> np.ndarray:
    return a * np.exp(-b * x) + c * x ** (-e)


def _fit(f, x: np.ndarray, y: np.ndarray, p0: list[float]) -> tuple[np.ndarray, float] | None:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            warnings.simplefilter("ignore", RuntimeWarning)
            popt, _ = curve_fit(f, x, y, p0=p0, bounds=(0.0, np.inf), maxfev=20000)
    except (RuntimeError, ValueError):
        return None
    rss = float(np.sum((f(x, *popt) - y) ** 2))
    return popt, rss


def fit_decay(lags: np.ndarray, mi: np.ndarray) -> dict:
    """Fit exponential a*exp(-b*x), power law a*x^(-e) and composite a*exp(-b*x) + c*x^(-e) to an MI decay.

    All three are fitted by (bounded, non-negative-parameter) least squares on the mutual information itself, not
    on its logarithm, so all residuals are on one scale and the AICs are comparable; this also keeps lags whose
    bias-corrected MI is zero or negative, which a log-space fit must drop. AIC = n ln(RSS / n) + 2 p (Gaussian
    errors of equal variance), p = 2, 2 and 4 parameters. The composite needs at least 6 points.
    Starting values come from straight-line fits in log space on the positive points.

    Inputs: lags (positive), mi (bits); non-finite points are dropped.
    Output dict: exponential {a, b, rss, aic}, power_law {a, e, rss, aic}, composite {a, b, c, e, rss, aic} or None,
    preferred ("exponential", "power_law", "composite" or "insufficient" with fewer than 4 points), n, caveat.
    """
    x = np.asarray(lags, dtype=float).ravel()
    y = np.asarray(mi, dtype=float).ravel()
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 4:
        return {"preferred": "insufficient", "n": n, "exponential": None, "power_law": None, "composite": None,
                "caveat": HIERARCHY_CAVEAT}

    def aic(rss: float, n_par: int) -> float:
        return float(n * np.log(max(rss, 1e-300) / n) + 2 * n_par)

    pos = y > 0
    a0 = float(max(y.max(), 1e-12))
    b0, e0 = 0.5, 1.0
    if pos.sum() >= 2:
        ly = np.log(y[pos])
        b0 = float(max(-np.polyfit(x[pos], ly, 1)[0], 1e-3))
        e0 = float(max(-np.polyfit(np.log(x[pos]), ly, 1)[0], 1e-3))
    out: dict = {"n": n, "caveat": HIERARCHY_CAVEAT}
    models: dict[str, float] = {}
    fe = _fit(_exp, x, y, [a0 * np.exp(b0 * x.min()), b0])
    out["exponential"] = None
    if fe is not None:
        (a, b), rss = fe
        out["exponential"] = {"a": float(a), "b": float(b), "rss": rss, "aic": aic(rss, 2)}
        models["exponential"] = out["exponential"]["aic"]
    fp = _fit(_power, x, y, [a0 * x.min() ** e0, e0])
    out["power_law"] = None
    if fp is not None:
        (a, e), rss = fp
        out["power_law"] = {"a": float(a), "e": float(e), "rss": rss, "aic": aic(rss, 2)}
        models["power_law"] = out["power_law"]["aic"]
    out["composite"] = None
    if n >= 6:
        ea = out["exponential"] or {"a": a0, "b": b0}
        pa = out["power_law"] or {"a": a0, "e": e0}
        fc = _fit(_composite, x, y, [ea["a"], max(ea["b"], 1e-3), max(pa["a"] * 0.1, 1e-12), max(pa["e"], 1e-3)])
        if fc is not None:
            (a, b, c, e), rss = fc
            out["composite"] = {"a": float(a), "b": float(b), "c": float(c), "e": float(e), "rss": rss,
                                "aic": aic(rss, 4)}
            models["composite"] = out["composite"]["aic"]
    out["preferred"] = min(models, key=models.get) if models else "insufficient"
    out["aic"] = models
    return out
