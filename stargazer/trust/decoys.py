"""Decoys: calibrating a search by running the identical search on data known to contain nothing.

Every detector in Stargazer reports a statistic. What that statistic means is decided here, by comparing it
with the same statistic computed on decoys: series simulated from the fitted null model at the real sampling,
or sequences scrambled so as to keep only what the null allows. Because the decoys go through the whole
procedure, including every fit and every maximisation over trial parameters, the comparison corrects for the
look-elsewhere effect inside one object without any formula.

Across objects, the per-object probabilities are combined by false-discovery-rate control.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import stats


@dataclass
class MonteCarloP:
    """A Monte Carlo probability: the fraction of decoys at least as extreme as the real object."""

    p: float
    draws: int
    exceedances: int
    stopped_early: bool

    @property
    def sigma(self) -> float:
        """The same probability on the one-sided Gaussian scale."""
        return float(stats.norm.isf(self.p))

    @property
    def resolution(self) -> float:
        """The smallest probability this number of decoys could have reported."""
        return 1.0 / (self.draws + 1)

    def to_dict(self) -> dict:
        return {"p": self.p, "sigma": self.sigma, "draws": self.draws, "exceedances": self.exceedances,
                "stopped_early": self.stopped_early, "resolution": self.resolution}


def sequential_p(observed: float, draw_batch: Callable[[int], np.ndarray], *, batch: int = 50,
                 max_draws: int = 2000, stop_exceedances: int = 10) -> MonteCarloP:
    """Monte Carlo probability with early stopping (Besag and Clifford's sequential scheme).

    `draw_batch(n)` returns n decoy statistics. Drawing stops as soon as `stop_exceedances` decoys have
    matched or beaten the observed statistic, because by then the object is plainly unremarkable; the estimate
    is then exceedances / draws. Otherwise all `max_draws` are used and the estimate is
    (exceedances + 1) / (draws + 1), which can never be zero. Unremarkable objects cost a few dozen decoys and
    remarkable ones get the full budget.
    """
    draws = exceedances = 0
    while draws < max_draws:
        values = np.asarray(draw_batch(min(batch, max_draws - draws)), dtype=float)
        hits = values >= observed
        if exceedances + int(hits.sum()) >= stop_exceedances:
            # stop at the exact draw where the count was reached, so the estimate is unbiased
            needed = stop_exceedances - exceedances
            position = int(np.flatnonzero(hits)[needed - 1]) + 1
            return MonteCarloP(stop_exceedances / (draws + position), draws + position, stop_exceedances, True)
        draws += values.size
        exceedances += int(hits.sum())
    return MonteCarloP((exceedances + 1) / (draws + 1), draws, exceedances, False)


def benjamini_hochberg(p: np.ndarray, q: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """False-discovery-rate control across objects.

    Returns (rejected, adjusted): a boolean array of the objects declared discoveries at rate q, and each
    object's adjusted probability (the smallest q at which it would be declared).
    """
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return np.zeros(0, dtype=bool), np.zeros(0)
    order = np.argsort(p)
    scaled = p[order] * n / np.arange(1, n + 1)
    adjusted_sorted = np.minimum.accumulate(scaled[::-1])[::-1].clip(max=1.0)
    adjusted = np.empty(n)
    adjusted[order] = adjusted_sorted
    return adjusted <= q, adjusted


def uniformity(p: np.ndarray) -> dict:
    """How far a set of probabilities is from uniform: the check that a null population is calibrated.

    Returns the observed and expected counts at or below 0.05 and 0.01, each with the number of binomial
    standard deviations by which the observed count departs from the expected one. Probabilities from
    sequential Monte Carlo are discrete, so a test that assumes continuous values (Kolmogorov–Smirnov) would
    reject a perfectly calibrated set; counts at thresholds are the right check.
    """
    p = np.asarray(p, dtype=float)
    out: dict = {"n": int(p.size)}
    for level in (0.05, 0.01):
        count, expected = int((p <= level).sum()), level * p.size
        spread = np.sqrt(p.size * level * (1 - level)) if p.size else float("nan")
        out[f"below_{level}"] = count
        out[f"expected_below_{level}"] = expected
        out[f"z_below_{level}"] = float((count - expected) / spread) if spread else float("nan")
    return out
