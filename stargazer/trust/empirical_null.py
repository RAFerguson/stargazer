"""An empirical null: letting the bulk of the population say what "ordinary" looks like.

Every probability in a search is computed against a model of the uninteresting case, and the model is never
exactly right. If it is slightly wrong in the same way for every object, the probabilities of the whole
population shift together: there are too many small ones, and a naive reading finds "discoveries" that are only
the model's error. Efron (2004, "Large-scale simultaneous hypothesis testing: the choice of a null hypothesis")
pointed out that when thousands of objects are tested and nearly all are uninteresting, the population measures
its own null. Convert each probability to a z-score; the centre of the z distribution should be a standard
normal curve; if it is instead centred at delta with width s, use that curve as the null.

Here the centre and width are fitted to the fractions of objects below a few moderate probability thresholds
(see `fit`), and the correction is only ever allowed to make results less significant: a width below one or a
centre below zero is ignored.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.stats import norm


@dataclass
class EmpiricalNull:
    centre: float            # fitted centre of the null z distribution (0 if the theoretical null holds)
    width: float             # fitted width (1 if the theoretical null holds)
    centre_used: float       # after the conservative floor
    width_used: float
    n: int

    def recalibrate(self, p: np.ndarray) -> np.ndarray:
        """Probabilities re-expressed against the empirical null."""
        z = to_z(p)
        return norm.sf((z - self.centre_used) / self.width_used)

    def to_dict(self) -> dict:
        return asdict(self)


def to_z(p: np.ndarray) -> np.ndarray:
    """z-scores for one-sided probabilities: small p gives large positive z."""
    return norm.isf(np.clip(np.asarray(p, float), 1e-300, 1 - 1e-12))


THRESHOLDS = (0.5, 0.25, 0.1, 0.05)


def fit(p: np.ndarray) -> EmpiricalNull:
    """Fit the null's centre and width to the fractions of probabilities at or below a few thresholds.

    If z-scores of uninteresting objects follow a normal curve with centre c and width w, the fraction with
    probability at or below a is the upper tail of that curve beyond the standard z-score of a. Turning each
    observed fraction back into a z-score u gives a straight line, z(a) = c + w u, which is fitted by weighted
    least squares. Fractions at thresholds are used, not quantiles of the z-scores, because sequential Monte
    Carlo probabilities are coarse where they are large and fine where they are small, and the thresholds
    chosen are ones at which those probabilities are exactly calibrated. With fewer than 200 objects the
    theoretical null (centre 0, width 1) is kept.
    """
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    n = p.size
    if n < 200:
        return EmpiricalNull(0.0, 1.0, 0.0, 1.0, int(n))
    levels = np.array(THRESHOLDS)
    fraction = np.clip(np.array([(p <= a).mean() for a in levels]), 0.5 / n, 1 - 0.5 / n)
    u = norm.isf(fraction)
    z = norm.isf(levels)
    weight = n * norm.pdf(u) ** 2 / (fraction * (1 - fraction))
    mean_u, mean_z = np.average(u, weights=weight), np.average(z, weights=weight)
    width = float(np.sum(weight * (u - mean_u) * (z - mean_z)) / np.sum(weight * (u - mean_u) ** 2))
    centre = float(mean_z - width * mean_u)
    return EmpiricalNull(centre, width, max(centre, 0.0), max(width, 1.0), int(n))
