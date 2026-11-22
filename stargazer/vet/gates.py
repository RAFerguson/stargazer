"""Vetting gates for a periodicity candidate.

A small probability is where vetting starts, not where it ends. Each gate asks one question that a real,
celestial, strictly periodic signal should survive and that a known kind of impostor should not. A gate
returns what it measured, a verdict, and one plain sentence. Gates never see whether an object is real, a
decoy twin or a plant (see `trust.blind`), so the false-alarm rate of the gates and the statistic together is
measured, not assumed, and no single threshold here has to be exactly right.

| Gate      | Question                                                        | Impostor it removes                    |
|-----------|-----------------------------------------------------------------|----------------------------------------|
| calendar  | Is the period one the Earth, Moon or survey would impose?       | seasonal and lunar systematics         |
| window    | Does the sampling pattern itself repeat at this period?         | aliases of the observing pattern       |
| seasons   | Does the signal survive removing any one year of data?          | one flare, one bad season              |
| bands     | Do the bands agree in phase?                                    | a fault in one filter or one chip      |
| holdout   | Do the early data alone find it, and predict the late data?     | red noise, which does not keep time    |
| coherence | Does the search do worse when the time axis is warped?          | noise the model got wrong, flicker     |
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from ..core import Series
from ..detect.calendar import calendar_periods
from ..detect.periodicity import (PeriodicityResult, SearchConfig, evaluate, fit_sinusoid, frequency_grid,
                                   harmonic_design, resolve_period)
from ..nulls import drw

@dataclass
class GateResult:
    gate: str
    passed: bool | None            # None: the gate could not be evaluated on these data
    sentence: str
    detail: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


def calendar_gate(period_days: float, span_days: float, tolerance: float = 0.5) -> GateResult:
    """Fail if the period is within `tolerance` frequency-resolution elements (1 / span) of a calendar period."""
    freq = 1.0 / period_days
    nearest, distance = "", np.inf
    for name, period in calendar_periods().items():
        d = abs(freq - 1.0 / period) * span_days
        if d < distance:
            nearest, distance = name, d
    passed = distance > tolerance
    sentence = (f"The period of {period_days:.1f} days is {distance:.1f} resolution elements from the nearest "
                f"calendar period ({nearest}); " + ("that is clear of it." if passed else "too close to tell apart."))
    return GateResult("calendar", passed, sentence, {"nearest": nearest, "distance_resolution_elements": distance,
                                                     "tolerance": tolerance})


def spectral_window(time: np.ndarray, error: np.ndarray, freq: np.ndarray) -> np.ndarray:
    """|sum of weights x exp(2 pi i f t)|^2 / (sum of weights)^2: how strongly the sampling itself repeats at f."""
    weight = 1.0 / error ** 2
    phase = 2.0 * np.pi * np.outer(freq, time)
    return np.abs((weight * np.exp(1j * phase)).sum(axis=1)) ** 2 / weight.sum() ** 2


def window_gate(series: Series, period_days: float, limit: float = 0.25) -> GateResult:
    """Fail if the sampling pattern has more than `limit` of its power at the candidate period."""
    freq = 1.0 / period_days
    resolution = 1.0 / series.span
    local = freq + resolution * np.linspace(-0.5, 0.5, 11)
    value = float(spectral_window(series.time, series.error, local).max())
    passed = value <= limit
    sentence = (f"The observing pattern has {100 * value:.0f}% of its power at this period; "
                + ("the sampling does not explain the signal." if passed else "the sampling alone could produce it."))
    return GateResult("window", passed, sentence, {"window_power": value, "limit": limit})


def _phase_distance(a: float, b: float) -> float:
    return abs((a - b + 0.5) % 1.0 - 0.5)


def season_gate(series: Series, period_days: float, sigma: float, tau: float, *, block_days: float = 365.2422,
                max_phase_shift: float = 0.15, max_share: float = 0.6) -> GateResult:
    """Remove each year of data in turn and refit. Pass if the phase never moves by more than
    `max_phase_shift` cycles, the amplitude never changes by more than a factor of two, and no single year
    carries more than `max_share` of the evidence. Years are counted from the first observation, so each block
    holds one observing season whether or not the season has gaps."""
    label = np.floor((series.time - series.time[0]) / block_days).astype(int)
    seasons = [s for s in np.unique(label) if np.sum(label == s) >= 5]
    if len(seasons) < 3:
        return GateResult("seasons", None, "Fewer than three years of data: the test cannot be made.",
                          {"n_seasons": len(seasons)})
    full = fit_sinusoid(series, period_days, sigma, tau)
    shifts, ratios, shares = [], [], []
    for s in seasons:
        rest = fit_sinusoid(series.select(label != s), period_days, sigma, tau)
        shifts.append(_phase_distance(rest["phase_cycles"], full["phase_cycles"]))
        ratios.append(rest["amplitude"] / full["amplitude"] if full["amplitude"] > 0 else np.inf)
        shares.append(1.0 - rest["statistic"] / full["statistic"] if full["statistic"] > 0 else 1.0)
    worst_shift, worst_share = float(max(shifts)), float(max(shares))
    ratio_low, ratio_high = float(min(ratios)), float(max(ratios))
    passed = worst_shift <= max_phase_shift and worst_share <= max_share and 0.5 <= ratio_low and ratio_high <= 2.0
    sentence = (f"Removing any one of {len(seasons)} years moves the phase by at most {worst_shift:.2f} cycles and "
                f"the amplitude by a factor of {ratio_low:.2f} to {ratio_high:.2f}; the most important year "
                f"carries {100 * worst_share:.0f}% of the evidence. "
                + ("The signal does not hang on one year." if passed else "The signal depends on one year."))
    return GateResult("seasons", passed, sentence, {"n_seasons": len(seasons), "max_phase_shift": worst_shift,
                                                    "amplitude_ratio_range": [ratio_low, ratio_high],
                                                    "max_share_of_evidence": worst_share})


def band_gate(bands: list[Series], period_days: float, config: SearchConfig = SearchConfig(), *,
              max_phase_difference: float = 0.15, min_statistic: float = 6.0) -> GateResult:
    """Fit the candidate period to each band on its own. Pass if every band in which the sinusoid is visible
    (fit improvement of at least `min_statistic`) agrees in phase within `max_phase_difference` cycles.

    Agreement between bands is necessary, not sufficient: a quasar's own flicker is common to all bands. A
    disagreement, though, points to a fault in one filter, one chip or one reference image.
    """
    fits = []
    for band in bands:
        binned = band.binned(config.bin_days) if config.bin_days > 0 else band
        if len(binned) < 20:
            continue
        null = drw.fit(binned)
        one = fit_sinusoid(binned, period_days, null.sigma_map, null.tau_map)
        fits.append({"band": band.band, "source_id": band.source_id, "n": len(binned), **one})
    visible = [f for f in fits if f["statistic"] >= min_statistic]
    detail = {"bands": [{k: f[k] for k in ("band", "n", "statistic", "amplitude", "phase_cycles")} for f in fits]}
    if len(visible) < 2:
        return GateResult("bands", None, "The signal is visible in fewer than two bands: no comparison possible.",
                          detail)
    worst = max(_phase_distance(a["phase_cycles"], b["phase_cycles"]) for a in visible for b in visible)
    passed = worst <= max_phase_difference
    detail["max_phase_difference"] = worst
    sentence = (f"The {len(visible)} bands in which the signal is visible agree in phase to within {worst:.2f} "
                "cycles." if passed else
                f"The bands disagree in phase by up to {worst:.2f} cycles, which a real source would not do.")
    return GateResult("bands", passed, sentence, detail)


def holdout_gate(series: Series, period_days: float, *, train_fraction: float = 0.7,
                 config: SearchConfig = SearchConfig()) -> GateResult:
    """Search the early data alone, then predict the late data.

    The first `train_fraction` of the time span is searched on its own, as if the rest did not exist: noise
    fit, periodogram, joint fit of the periodic model, choice between a period and its harmonics. The gate
    passes only if (a) that search picks the same period, to within one frequency-resolution element of the
    early data, and (b) the early periodic model, extended unchanged, makes the late data more probable than
    no periodic model does. Both sides of (b) use the noise fitted to the early data alongside the periodic
    model, so the comparison is only about whether the extended model helps. Nothing about the late data
    enters the model, so this is a real prediction. Noise that does not keep time fails (b) more often than
    not.
    """
    cut = series.time[0] + train_fraction * series.span
    early_mask = series.time <= cut
    early = series.select(early_mask)
    late_n = len(series) - len(early)
    if len(early) < config.min_points or late_n < 15:
        return GateResult("holdout", None, "Too few points to hold some back.", {"n_early": len(early),
                                                                                "n_late": late_n})
    freq = frequency_grid(early.span, config)
    resolution = 1.0 / early.span
    if freq.size < 10 or 1.0 / period_days < freq[0] - resolution:
        return GateResult("holdout", None, "The early data cover too few cycles of this period to find it "
                          "independently.", {"early_span_days": early.span})
    sigma_axis, tau_axis = drw.default_grids(early.time, early.value, early.error, n_sigma=config.n_sigma,
                                             n_tau=config.n_tau)
    found = evaluate(early.time, early.error, early.value, freq, sigma_axis, tau_axis, config)
    early_freq, coef, sigma, tau = resolve_period(early, freq, found, sigma_axis, tau_axis, config)
    early_period = 1.0 / early_freq
    same = abs(early_freq - 1.0 / period_days) <= resolution
    design = harmonic_design(series.time[None, :], np.array([early_freq]), config.harmonics,
                             1.0 / config.min_period_days)
    model = design[0, :, 1:] @ coef[1:]
    pair = (np.array([sigma]), np.array([tau]))

    def late_log_like(value: np.ndarray) -> float:
        whole, _ = drw.log_likelihood_grid(series.time, series.error, value, *pair)
        first, _ = drw.log_likelihood_grid(series.time[early_mask], series.error[early_mask], value[early_mask],
                                           *pair)
        return float(whole[0] - first[0])

    gain = late_log_like(series.value - model) - late_log_like(series.value)
    passed = bool(same and gain > 0)
    share = f"the first {100 * train_fraction:.0f}% of the time span"
    if not same:
        sentence = (f"Searched alone, {share} prefers a period of {early_period:.1f} days, not "
                    f"{period_days:.1f}: the early data do not find this signal independently.")
    else:
        sentence = (f"Searched alone, {share} finds the same period ({early_period:.1f} days), and its periodic "
                    f"model predicts the remaining {late_n} points " + ("better" if gain > 0 else "worse")
                    + f" than no periodic model (log-likelihood gain {gain:+.1f}).")
    return GateResult("holdout", passed, sentence, {"log_likelihood_gain": gain, "early_period_days": early_period,
                                                    "same_period": bool(same), "n_early": len(early),
                                                    "n_late": late_n, "train_fraction": train_fraction})


def coherence_gate(result: PeriodicityResult, limit: float = 0.05) -> GateResult:
    """Pass if searches of the same data in warped time rarely do as well as the search in true time."""
    p = result.coherence_p
    if not np.isfinite(p):
        return GateResult("coherence", None, "The warped-time test was not run.", {})
    passed = p <= limit
    sentence = (f"With the time axis warped at random, the search does as well as in true time {100 * p:.1f}% of "
                "the time; " + ("the variation keeps time." if passed else "the variation does not clearly keep time."))
    return GateResult("coherence", passed, sentence, {"coherence_p": p, "limit": limit})


def run_gates(merged: Series, bands: list[Series], result: PeriodicityResult,
              config: SearchConfig = SearchConfig()) -> list[GateResult]:
    """All gates for one candidate. `merged` is the series that was searched (before binning); `bands` are its
    separate bands; `result` is the search result."""
    series = merged.binned(config.bin_days) if config.bin_days > 0 else merged
    sigma, tau = result.noise["sigma_with_signal"], result.noise["tau_with_signal_days"]
    period = result.period_days
    return [calendar_gate(period, series.span, config.calendar_tolerance),
            window_gate(series, period),
            season_gate(series, period, sigma, tau),
            band_gate(bands, period, config),
            holdout_gate(series, period, config=config),
            coherence_gate(result)]


def verdict(gates: list[GateResult]) -> tuple[bool, list[str]]:
    """A candidate is vetted if no gate failed. Gates that could not be evaluated are listed, not counted as
    passes: the dossier says which checks were not possible."""
    failed = [g.gate for g in gates if g.passed is False]
    return not failed, failed
