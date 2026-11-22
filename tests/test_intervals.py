"""Tests for stargazer.detect.intervals: the clock test, the spacing test and slot bits, on the benchmark items,
on synthetic slotted and pulse-position-coded events, and calibration on Poisson and bursty event lists."""
from __future__ import annotations

import numpy as np
import pytest

from stargazer.bench import build_zoo
from stargazer.core import BitStream, EventList
from stargazer.core.rng import rng_for
from stargazer.detect.intervals import (ClockResult, SpacingResult, clock_test, jitter_surrogates,
                                        rayleigh_spectrum, rotation_surrogates, spacing_test, to_bits)
from stargazer.trust.decoys import MonteCarloP

SEED = 20261122


@pytest.fixture(scope="module")
def zoo() -> dict:
    return {item.name: item for item in build_zoo(seed=0, collections=("natural", "human", "control"))}


# ---------------------------------------------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------------------------------------------

def poisson_times(rng: np.random.Generator, n: int) -> np.ndarray:
    return np.cumsum(rng.exponential(1.0, n))


def switching_times(rng: np.random.Generator, n: int, rates: tuple[float, float] = (0.3, 3.0)) -> np.ndarray:
    """A doubly stochastic Poisson process whose rate alternates between two values, each spell lasting an
    exponential time of about 50 mean intervals (at the average rate)."""
    mean_rate = 0.5 * (rates[0] + rates[1])
    out: list[float] = []
    t, high = 0.0, bool(rng.integers(2))
    while len(out) < n:
        spell = rng.exponential(50.0 / mean_rate)
        rate = rates[1] if high else rates[0]
        out.extend(t + np.sort(rng.uniform(0.0, spell, rng.poisson(rate * spell))))
        t += spell
        high = not high
    return np.array(out[:n])


def weibull_times(rng: np.random.Generator, n: int, shape: float = 0.6) -> np.ndarray:
    """A clustered renewal process: Weibull intervals with shape below one."""
    return np.cumsum(rng.weibull(shape, n))


def ppm_times(rng: np.random.Generator, n: int, levels: int, jitter: float) -> np.ndarray:
    """Pulse-position coding: each interval is a whole number 1..levels of unit quanta, plus timing jitter that
    accumulates (there is no global clock)."""
    return np.cumsum(1.0 + rng.integers(0, levels, n) + jitter * rng.standard_normal(n))


def p_of(result: ClockResult | SpacingResult) -> float:
    return result.p_value.p if isinstance(result.p_value, MonteCarloP) else float(result.p_value)


# ---------------------------------------------------------------------------------------------------------------
# building blocks
# ---------------------------------------------------------------------------------------------------------------

def test_rayleigh_spectrum_peaks_at_clock_and_matches_exact_sum() -> None:
    rng = rng_for(SEED, "rayleigh")
    ticks = np.flatnonzero(rng.random(800) < 0.5)
    times = 3.0 + 2.5 * ticks + 0.01 * rng.standard_normal(ticks.size)
    freq, power = rayleigh_spectrum(times, 1.0, 10.0)
    assert np.all(np.diff(freq) > 0)
    best = 1.0 / freq[np.argmax(power)]
    assert best == pytest.approx(2.5, rel=1e-3) or best == pytest.approx(2.5 / 2, rel=1e-3) or \
        best == pytest.approx(2.5 / 3, rel=1e-3)
    # the binned FFT agrees with the exact sum to within the binning loss
    f0 = freq[np.argmax(power)]
    exact = np.abs(np.exp(2j * np.pi * f0 * (times - times[0])).sum()) ** 2 * 2 / times.size
    assert power.max() == pytest.approx(exact, rel=0.08)
    # random times: power is chi-squared with two degrees of freedom, mean 2
    _, noise = rayleigh_spectrum(np.sort(rng.uniform(0, 1000, 1000)), 0.1, 1.0)
    assert 1.8 < noise.mean() < 2.2


@pytest.mark.parametrize("make", [jitter_surrogates, rotation_surrogates])
def test_surrogates_keep_window_counts_and_destroy_a_clock(make) -> None:
    rng = rng_for(SEED, "surrogates", make.__name__)
    times = np.arange(2000) * 1.0 + 0.3 + 0.002 * rng.standard_normal(2000)
    fake = make(times, 10.0, rng, 100)
    assert fake.shape == (100, times.size)
    assert np.all(np.diff(fake, axis=1) >= 0)
    # the clock is gone: Rayleigh power at the true period (about 4000 in the data) drops, on average, to the
    # level the ten events of one window carry by themselves (about 20, rotation) or to noise (about 2, jitter)
    original = np.abs(np.exp(2j * np.pi * times).sum()) ** 2 * 2 / times.size
    z = np.abs(np.exp(2j * np.pi * fake).sum(axis=1)) ** 2 * 2 / times.size
    assert z.mean() < 0.01 * original
    fake = fake[:4]
    # slow rate changes survive: counts in 100-unit blocks are kept to within a window's worth of events
    counts = np.histogram(times, bins=np.arange(0, 2001, 100))[0]
    for row in fake:
        assert np.all(np.abs(np.histogram(row, bins=np.arange(0, 2001, 100))[0] - counts) <= 20)


def test_rotation_surrogates_keep_short_intervals() -> None:
    rng = rng_for(SEED, "rotation")
    times = weibull_times(rng, 3000)
    fake = rotation_surrogates(times, 20.0, rng, 3)
    short = np.sort(np.diff(times))[:1500]
    for row in fake:
        assert np.allclose(np.sort(np.diff(row))[:1400], short[:1400], rtol=0.2, atol=1e-3)


# ---------------------------------------------------------------------------------------------------------------
# clock test on the benchmark
# ---------------------------------------------------------------------------------------------------------------

def test_ook_pulses_clock_and_bits(zoo) -> None:
    item = zoo["human/ook_pulses"]
    result = clock_test(item.payload, rng_for(SEED, "ook"))
    assert result.detected
    assert result.period == pytest.approx(item.truth["slot_length"], rel=1e-4)
    assert result.residual_rms < 0.02
    bits = to_bits(item.payload, result.period, result.phase)
    assert isinstance(bits, BitStream)
    assert np.array_equal(bits.bits, item.truth["bits"])
    assert bits.meta["n_far"] == 0
    assert result.occupancy == pytest.approx(item.truth["bits"].mean(), abs=1e-9)


def test_pulsar_is_a_clock_and_bits_reproduce_nulling(zoo) -> None:
    item = zoo["natural/pulsar_nulling"]
    result = clock_test(item.payload, rng_for(SEED, "pulsar"))
    assert result.detected                   # a pulsar is a clock; content analysis judges its on/off sequence
    assert result.period == pytest.approx(item.truth["period"], rel=1e-5)
    assert result.residual_rms < 0.01
    bits = to_bits(item.payload, result.period, result.phase)
    assert np.array_equal(bits.bits, item.truth["on_off"])
    assert bits.meta["slot0_time"] == pytest.approx(item.truth["t0"], abs=0.01)


def test_repeating_burster_has_no_clock_and_no_quantum(zoo) -> None:
    events = zoo["natural/repeating_burster"].payload
    assert not clock_test(events, rng_for(SEED, "burster")).detected
    assert not spacing_test(events, rng_for(SEED, "burster spacing")).detected


def test_poisson_control_has_no_clock_and_no_quantum(zoo) -> None:
    events = zoo["control/poisson_events"].payload
    clock = clock_test(events, rng_for(SEED, "poisson"))
    spacing = spacing_test(events, rng_for(SEED, "poisson spacing"))
    assert not clock.detected and p_of(clock) > 0.05
    assert not spacing.detected and p_of(spacing) > 0.05


def test_harmonic_choice_prefers_the_slot_length() -> None:
    # sparse occupancy: the strongest power can fall on a sub-multiple; the reported period must be the slot
    rng = rng_for(SEED, "harmonic")
    ticks = np.flatnonzero(rng.random(3000) < 0.3)
    events = EventList(10.0 + 0.4 * ticks + 0.001 * rng.standard_normal(ticks.size))
    result = clock_test(events, rng, n_surrogates=50, alpha=0.05)
    assert result.detected
    assert result.period == pytest.approx(0.4, rel=1e-4)
    assert result.occupancy == pytest.approx(0.3, abs=0.03)


def test_to_bits_marks_exposure_gaps_unknown() -> None:
    bits_true = np.array([1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1], dtype=np.uint8)
    keep = bits_true.astype(bool)
    keep[6:9] = False                         # slots 6-8 fall in a gap
    times = np.flatnonzero(keep) * 2.0 + 1.0
    events = EventList(times, exposure=[(0.0, 12.0), (18.0, 32.0)])
    bits = to_bits(events, 2.0, 1.0)
    assert len(bits) == 16
    assert bits.known is not None and not bits.known[6:9].any() and bits.known[:6].all()
    assert np.array_equal(bits.bits[bits.known], bits_true[bits.known])


# ---------------------------------------------------------------------------------------------------------------
# spacing test
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("levels", [2, 4, 16])
def test_spacing_finds_pulse_position_quantum(levels: int) -> None:
    rng = rng_for(SEED, "ppm", levels)
    events = EventList(ppm_times(rng, 400, levels, 0.01))
    result = spacing_test(events, rng)
    assert result.detected
    assert result.quantum == pytest.approx(1.0, rel=2e-3)
    assert result.n_levels == levels
    assert np.array_equal(result.levels, np.arange(1, levels + 1))


def test_spacing_survives_where_the_clock_is_lost() -> None:
    # jitter accumulates, so there is no global clock, but every interval is still a multiple of the quantum
    rng = rng_for(SEED, "ppm drift")
    events = EventList(ppm_times(rng, 400, 8, 0.1))
    assert spacing_test(events, rng).detected
    assert not clock_test(events, rng, n_surrogates=100).detected


# ---------------------------------------------------------------------------------------------------------------
# calibration
# ---------------------------------------------------------------------------------------------------------------

def _rates(make, n_lists: int, n_events: int, tag: str, n_surrogates: int) -> tuple[float, float]:
    clock_p, spacing_p = [], []
    for i in range(n_lists):
        rng = rng_for(SEED, "calibration", tag, i)
        events = EventList(make(rng, n_events))
        clock_p.append(p_of(clock_test(events, rng, n_surrogates=n_surrogates)))
        spacing_p.append(p_of(spacing_test(events, rng, n_surrogates=n_surrogates)))
    return float(np.mean(np.array(clock_p) <= 0.05)), float(np.mean(np.array(spacing_p) <= 0.05))


@pytest.mark.parametrize("make,tag", [(poisson_times, "poisson"), (switching_times, "switching"),
                                      (weibull_times, "weibull")])
def test_false_alarm_rate_quick(make, tag: str) -> None:
    # 40 short lists each: expected 2 below 0.05; allow up to 6 (binomial tail probability about 1.4%)
    clock, spacing = _rates(make, 40, 250, tag, 100)
    assert clock <= 6 / 40 and spacing <= 6 / 40


@pytest.mark.slow
@pytest.mark.parametrize("make,tag,n_lists", [(poisson_times, "poisson", 200), (switching_times, "switching", 100),
                                              (weibull_times, "weibull", 100)])
def test_false_alarm_rate_full(make, tag: str, n_lists: int) -> None:
    clock, spacing = _rates(make, n_lists, 500, tag + " full", 200)
    print(f"{tag}: fraction with p <= 0.05: clock {clock:.3f}, spacing {spacing:.3f} ({n_lists} lists)")
    limit = 0.05 + 2.5 * np.sqrt(0.05 * 0.95 / n_lists)
    assert clock <= limit and spacing <= limit
