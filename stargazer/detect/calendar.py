"""The calendar: periods that the Earth, the Moon and a ground-based survey imprint on every light curve.

Ground-based photometry carries systematic errors that follow the lunar month (bright sky near full Moon) and
the year (airmass, season, reference-image effects). In an engineering check on 47 ZTF quasars in a field that
is not part of any registered search, the stacked periodogram showed a 29.5-day line in 60% of the objects, at
a mean statistic of 30 where 2 is expected, with a second line at 27.3 days (the lunar month beating with the
year) and a weaker one at 32.1 days. That is the Moon, not the quasars. The defences:

* the search is confined to periods long enough that neither a trial period nor any harmonic of the periodic
  model reaches the lunar band (`SearchConfig.min_period_days`, 40 days by default);
* trial periods near the year and its fractions are never chosen as candidates, nor are periods whose second
  or third harmonic would fall on the year or the half-year;
* a vetting gate rejects any candidate that ends up near a calendar period all the same.

Fitting and subtracting sinusoids at the lunar and annual frequencies was tried and rejected. It halved the
lunar line but did not remove it (the line's strength varies from year to year), and it raised the statistic
at every neighbouring frequency where a sinusoid had been subtracted, because the subtraction is made in one
metric and the search in another. Staying away is safer than cleaning.

The cost is stated plainly: a real signal at one of these periods, or shorter than the lunar band, cannot be
found from the ground this way.
"""
from __future__ import annotations

import numpy as np

YEAR = 365.2422
SYNODIC_MONTH = 29.530589
SIDEREAL_MONTH = 27.321662


def calendar_periods() -> dict[str, float]:
    """Periods a ground-based observing pattern can imprint, in days, by name."""
    out = {"one year": YEAR}
    out.update({f"one year / {k}": YEAR / k for k in (2, 3, 4)})
    out["synodic month"] = SYNODIC_MONTH
    out["sidereal month"] = SIDEREAL_MONTH
    for sign, name in ((+1, "plus"), (-1, "minus")):
        out[f"synodic month {name} one cycle per year"] = 1.0 / (1.0 / SYNODIC_MONTH + sign / YEAR)
    return out


def mask(freq: np.ndarray, span_days: float, tolerance: float = 0.5, harmonics: int = 3) -> np.ndarray:
    """True for trial frequencies to keep away from: within `tolerance` resolution elements (1 / span) of a
    calendar period, or with a harmonic (up to `harmonics`) that close to the year or the half-year."""
    freq = np.asarray(freq, float)
    out = np.zeros(freq.size, dtype=bool)
    for period in calendar_periods().values():
        out |= np.abs(freq - 1.0 / period) * span_days <= tolerance
    for period in (YEAR, YEAR / 2):
        for k in range(2, harmonics + 1):
            out |= np.abs(k * freq - 1.0 / period) * span_days <= tolerance
    return out
