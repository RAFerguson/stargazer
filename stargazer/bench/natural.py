"""The natural collection: signals structured by natural processes. A tool must not call these artificial.

Each item is structured relative to white noise or independent symbols, and each names in `truth["fair_null"]`
the natural model against which it should NOT stand out. Content tiers stay at or below 2: continuous red-noise
and periodic light curves reach tier 0 (no discrete states); the eclipsing binary and the pulsar's on/off states
reach tier 1 (two discrete states recovered); none carries structure beyond its natural model.
"""
from __future__ import annotations

import math

import numpy as np

from stargazer.core import BitStream, ContentTier, EventList, Series

from .controls import markov_bits
from .zoo import ZooItem, expect, register

MJD0 = 58_000.0
YEAR = 365.25


def survey_times(rng: np.random.Generator, n_seasons: int = 7, season_days: float = 240.0,
                 nights: tuple[int, int] = (30, 45)) -> np.ndarray:
    """Realistic ground-survey cadence: yearly seasons of about 8 months, 30 to 45 distinct nights per season
    chosen at random, one visit per night at a random time of night (MJD)."""
    times = []
    for s in range(n_seasons):
        n = int(rng.integers(nights[0], nights[1] + 1))
        night = np.sort(rng.choice(int(season_days), size=n, replace=False))
        times.append(MJD0 + s * YEAR + night + rng.uniform(0.15, 0.40, n))
    return np.concatenate(times)


def simulate_drw(time: np.ndarray, sigma: float, tau: float, rng: np.random.Generator) -> np.ndarray:
    """An exact damped-random-walk (Ornstein-Uhlenbeck) realisation at the given times, zero mean, stationary
    standard deviation sigma and timescale tau: x[i] = a x[i-1] + sqrt(sigma^2 (1 - a^2)) z, a = exp(-dt/tau)."""
    time = np.asarray(time, dtype=float)
    x = np.empty(time.size)
    z = rng.standard_normal(time.size)
    x[0] = sigma * z[0]
    a = np.exp(-np.diff(time) / tau)
    scale = sigma * np.sqrt(1.0 - a ** 2)
    for i in range(1, time.size):
        x[i] = a[i - 1] * x[i - 1] + scale[i - 1] * z[i]
    return x


def _survey_series(name: str, time: np.ndarray, model: np.ndarray, rng: np.random.Generator) -> Series:
    error = rng.uniform(0.01, 0.03, time.size)
    value = model + error * rng.standard_normal(time.size)
    return Series(time, value, error, band="g", source_id=name, kind="mag", time_standard="MJD")


def _seasons(size: str) -> int:
    return 7 if size == "small" else 12


@register("natural/drw_quasar")
def drw_quasar(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A quasar-like damped-random-walk light curve at seasonal survey sampling, around magnitude 17 with errors
    of 0.01 to 0.03 mag. truth: sigma, tau, mean_mag, the noise-free curve (`model`), fair_null "drw"."""
    name = "natural/drw_quasar"
    sigma, tau, mean = 0.15, 300.0, 17.0
    time = survey_times(rng, _seasons(size))
    model = mean + simulate_drw(time, sigma, tau, rng)
    truth = {"sigma": sigma, "tau": tau, "mean_mag": mean, "model": model, "error_range": (0.01, 0.03),
             "fair_null": "drw"}
    return ZooItem(name, "natural", _survey_series(name, time, model, rng), truth,
                   expect(True, ContentTier.NONE, "natural"),
                   "A damped-random-walk quasar light curve (sigma 0.15 mag, tau 300 d) at survey sampling.")


@register("natural/variable_star_p180")
def variable_star(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A strictly periodic sinusoidal variable on white noise at survey sampling. truth: period, amplitude
    (semi-amplitude in mag), phase (radians, model = mean + A sin(2 pi t / P + phase)), model, fair_null
    "periodic_plus_noise"."""
    name = "natural/variable_star_p180"
    period, amplitude, mean = 180.0, 0.2, 17.0
    phase = float(rng.uniform(0, 2 * np.pi))
    time = survey_times(rng, _seasons(size))
    model = mean + amplitude * np.sin(2 * np.pi * time / period + phase)
    truth = {"period": period, "amplitude": amplitude, "phase": phase, "mean_mag": mean, "model": model,
             "fair_null": "periodic_plus_noise"}
    return ZooItem(name, "natural", _survey_series(name, time, model, rng), truth,
                   expect(True, ContentTier.NONE, "natural"),
                   "A strictly periodic variable star (P = 180 d, 0.2 mag) on white noise at survey sampling.")


@register("natural/drw_plus_weak_sine")
def drw_plus_weak_sine(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """The hard case: a damped random walk (sigma 0.15, tau 300 d) plus a weak sinusoid (0.05 mag, 400 d).
    It does stand out from a pure-DRW null, weakly, and is still natural. truth: period, amplitude, phase,
    sigma, tau, model, fair_null "drw_plus_periodic"."""
    name = "natural/drw_plus_weak_sine"
    sigma, tau, period, amplitude, mean = 0.15, 300.0, 400.0, 0.05, 17.0
    phase = float(rng.uniform(0, 2 * np.pi))
    time = survey_times(rng, _seasons(size))
    model = mean + simulate_drw(time, sigma, tau, rng) + amplitude * np.sin(2 * np.pi * time / period + phase)
    truth = {"period": period, "amplitude": amplitude, "phase": phase, "sigma": sigma, "tau": tau, "mean_mag": mean,
             "model": model, "fair_null": "drw_plus_periodic", "stands_out_from_drw": True}
    return ZooItem(name, "natural", _survey_series(name, time, model, rng), truth,
                   expect(True, ContentTier.NONE, "natural"),
                   "A damped random walk plus a weak 400-day sinusoid: periodicity hidden in red noise.")


@register("natural/eclipsing_binary")
def eclipsing_binary(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """An eclipsing-binary-like light curve in normalised flux at even 30-minute sampling: flat at 1 with box
    dips (period 3.7 d, depth 2%, duration 3 h) plus white noise. truth: period, depth, duration (days), epoch
    (centre of the first eclipse), noise sigma, in-eclipse mask, model, fair_null "periodic_plus_noise"."""
    name = "natural/eclipsing_binary"
    period, depth, duration, noise = 3.7, 0.02, 3.0 / 24.0, 0.002
    days = 90.0 if size == "small" else 360.0
    cadence = 30.0 / (24 * 60)
    time = MJD0 + np.arange(int(round(days / cadence))) * cadence
    epoch = MJD0 + float(rng.uniform(0, period))
    phase = ((time - epoch) / period + 0.5) % 1.0 - 0.5
    in_eclipse = np.abs(phase * period) < duration / 2
    model = np.where(in_eclipse, 1.0 - depth, 1.0)
    value = model + noise * rng.standard_normal(time.size)
    series = Series(time, value, np.full(time.size, noise), band="TESS", source_id=name, kind="flux",
                    time_standard="MJD")
    truth = {"period": period, "depth": depth, "duration": duration, "epoch": epoch, "noise_sigma": noise,
             "cadence": cadence, "in_eclipse": in_eclipse, "model": model, "fair_null": "periodic_plus_noise"}
    return ZooItem(name, "natural", series, truth, expect(True, ContentTier.SYMBOLS, "natural"),
                   "Box-shaped periodic eclipses (P = 3.7 d, 2% deep, 3 h long) at 30-minute sampling.")


PULSAR_TRANSITION = np.array([[0.7, 0.3], [0.1, 0.9]])  # rows: from (off, on); columns: to (off, on)


def pulsar_states(rng: np.random.Generator, n_periods: int) -> np.ndarray:
    """The pulsar's nulling sequence: 1 = pulse on, 0 = null, a two-state Markov chain with P(on->on) = 0.9 and
    P(off->off) = 0.7, started from its stationary distribution."""
    return markov_bits(rng, n_periods, PULSAR_TRANSITION)


def _pulsar_truth(states: np.ndarray, period: float) -> dict:
    return {"period": period, "on_off": states, "transition_matrix": PULSAR_TRANSITION.copy(),
            "state_order": ("off", "on"), "p_stay_on": 0.9, "p_stay_off": 0.7,
            "p_on_stationary": 0.3 / 0.4, "n_periods": int(states.size)}


@register("natural/pulsar_nulling")
def pulsar_nulling(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A nulling pulsar as an EventList: pulses expected every 0.714 s at t0 + k P with 2 ms Gaussian timing
    jitter and lognormal amplitudes, absent when the Markov nulling chain is off. truth: period, t0, jitter, the
    on/off sequence, the transition matrix, the amplitude law, the period index of each pulse, fair_null
    "periodic_with_markov_nulling"."""
    name = "natural/pulsar_nulling"
    period, jitter, t0 = 0.714, 0.002, 1.0
    n_periods = 5_000 if size == "small" else 20_000
    states = pulsar_states(rng, n_periods)
    index = np.flatnonzero(states)
    time = t0 + index * period + jitter * rng.normal(size=index.size)
    amplitude = rng.lognormal(mean=0.0, sigma=0.5, size=index.size)
    events = EventList(time, amplitude=amplitude, exposure=[(0.0, t0 + n_periods * period)], source_id=name,
                       time_unit="s")
    truth = _pulsar_truth(states, period)
    truth.update({"t0": t0, "jitter": jitter, "amplitude_law": "lognormal(mu=0, sigma=0.5)", "pulse_index": index,
                  "fair_null": "periodic_with_markov_nulling"})
    return ZooItem(name, "natural", events, truth, expect(True, ContentTier.SYMBOLS, "natural"),
                   "A 0.714 s pulsar with timing jitter, lognormal amplitudes and Markov nulling.")


@register("natural/pulsar_nulling_bits")
def pulsar_nulling_bits(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A nulling pulsar's on/off sequence (1 = pulse, 0 = null) as a BitStream, drawn from the same Markov chain
    as natural/pulsar_nulling (with this item's own seed). truth: the transition matrix and sequence, fair_null
    "markov1"."""
    name = "natural/pulsar_nulling_bits"
    states = pulsar_states(rng, 5_000 if size == "small" else 20_000)
    truth = _pulsar_truth(states, 0.714)
    truth["fair_null"] = "markov1"
    # A payload that is already a bit stream has no discreteness left to discover, so the justified tier is 0:
    # the on/off sequence is a two-state chain and nothing more.
    return ZooItem(name, "natural", BitStream(states), truth, expect(True, ContentTier.NONE, "natural"),
                   "A pulsar's nulling on/off sequence: Markov-1, structured only against independent bits.")


def _merge_slots(chosen: np.ndarray, start: float, slot: float) -> list[tuple[float, float]]:
    """Merge chosen equal-width slots (boolean array) into maximal (start, stop) intervals."""
    edges = np.diff(np.concatenate([[0], chosen.astype(int), [0]]))
    begins, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    return [(start + b * slot, start + e * slot) for b, e in zip(begins, ends)]


def _intersect(a: list[tuple[float, float]], b: list[tuple[float, float]]) -> list[tuple[float, float]]:
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


@register("natural/repeating_burster")
def repeating_burster(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A repeating-burster-like EventList (in days): a clustered Weibull renewal process (shape 0.6) running only
    inside periodic activity windows (period 16.35 d, 5-day window) and only seen during exposure (random half-day
    blocks covering 60% of the time). The renewal clock runs on observable active time, so every burst lies inside
    both a window and an exposure interval. truth: period, window, window_start (phase of the first window),
    shape, scale, exposure fraction, active intervals, fair_null "windowed_weibull_renewal"."""
    name = "natural/repeating_burster"
    period, window, shape, n_bursts = 16.35, 5.0, 0.6, 400 if size == "small" else 2_000
    duration = 1_000.0 if size == "small" else 5_000.0
    slot, exposure_fraction = 0.5, 0.6
    n_slots = int(duration / slot)
    chosen = np.zeros(n_slots, dtype=bool)
    chosen[rng.choice(n_slots, size=int(round(exposure_fraction * n_slots)), replace=False)] = True
    exposure = _merge_slots(chosen, MJD0, slot)
    window_start = MJD0 + float(rng.uniform(0, period))
    first = window_start - period * math.ceil((window_start - MJD0) / period)
    windows = [(max(s, MJD0), min(s + window, MJD0 + duration))
               for s in np.arange(first, MJD0 + duration, period) if s + window > MJD0]
    active = _intersect(windows, exposure)
    lengths = np.array([hi - lo for lo, hi in active])
    total = float(lengths.sum())
    # Renewal intervals of unit scale, rescaled so that n_bursts fit in the observable active time.
    gaps = rng.weibull(shape, n_bursts + 1)
    scale = total / gaps.sum()
    clock = np.cumsum(gaps[:-1]) * scale
    cum = np.concatenate([[0.0], np.cumsum(lengths)])
    k = np.clip(np.searchsorted(cum, clock, side="right") - 1, 0, len(active) - 1)
    starts = np.array([lo for lo, _ in active])
    time = starts[k] + (clock - cum[k])
    events = EventList(time, exposure=exposure, source_id=name, time_unit="d")
    truth = {"period": period, "window": window, "window_start": window_start, "shape": shape, "scale": scale,
             "n_bursts": n_bursts, "duration": duration, "exposure_fraction": exposure_fraction,
             "active_intervals": active, "fair_null": "windowed_weibull_renewal"}
    return ZooItem(name, "natural", events, truth, expect(True, ContentTier.NONE, "natural"),
                   "Clustered Weibull bursts inside 5-day windows of a 16.35-day activity cycle, with 60% exposure.")
