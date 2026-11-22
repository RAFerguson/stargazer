"""Data-quality steps applied before any search: honest error bars and removal of wild points.

Both steps are fixed rules, applied identically to every object, and recorded in the run configuration.
"""
from __future__ import annotations

import numpy as np

from ..core import Series


def pair_residuals(series: Series, max_gap_days: float = 0.5) -> np.ndarray:
    """Differences between consecutive measurements taken within `max_gap_days`, in units of their combined
    quoted error. A quasar does not change in a few hours, so these should scatter with unit variance; if they
    scatter more, the quoted errors are too small by that factor."""
    gap = np.diff(series.time)
    close = gap <= max_gap_days
    combined = np.sqrt(series.error[1:] ** 2 + series.error[:-1] ** 2)
    return (np.diff(series.value) / combined)[close]


def error_scale(residuals: np.ndarray, min_pairs: int = 20) -> float | None:
    """The factor by which quoted errors understate the scatter of same-night pairs (robust), or None when
    there are too few pairs to say."""
    if residuals.size < min_pairs:
        return None
    return float(1.4826 * np.median(np.abs(residuals - np.median(residuals))))


def rescale_errors(series: Series, scale: float) -> Series:
    """Errors multiplied by max(scale, 1): errors are only ever enlarged."""
    factor = max(float(scale), 1.0)
    out = Series(series.time, series.value, series.error * factor, series.band, series.source_id, series.kind,
                 series.time_standard, dict(series.meta), series.provenance)
    out.meta["error_scale"] = factor
    return out


def clip_outliers(series: Series, n_sigma: float = 6.0) -> Series:
    """Drop points further than `n_sigma` robust scatters from the median, where the scatter of a point is its
    own error added in quadrature to the series' robust scatter. Only wild points go: with the default, a
    genuine flare of a few tenths of a magnitude stays."""
    centre = np.median(series.value)
    scatter = 1.4826 * np.median(np.abs(series.value - centre))
    keep = np.abs(series.value - centre) <= n_sigma * np.sqrt(scatter ** 2 + series.error ** 2)
    out = series.select(keep)
    out.meta["n_clipped"] = int(np.sum(~keep))
    return out
