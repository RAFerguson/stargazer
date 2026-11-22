"""Tests for stargazer.decode.symbolise: discrete levels on the benchmark light curves, keyed multi-level
series, the traps (sinusoid, sawtooth, red noise), unresolved levels, symbol probabilities, the lossy
quantiser, interval levels, and calibration against damped random walks at survey sampling."""
from __future__ import annotations

import numpy as np
import pytest

from stargazer.bench import build_zoo
from stargazer.bench.natural import simulate_drw, survey_times
from stargazer.core import EventList, Series, SymbolSequence
from stargazer.core.rng import rng_for
from stargazer.decode.symbolise import (LevelsResult, discreteness_test, event_symbols, find_levels, quantise,
                                        symbolise)

SEED = 20261122


@pytest.fixture(scope="module")
def zoo() -> dict:
    return {item.name: item for item in build_zoo(seed=0, collections=("natural", "control"))}


def keyed_series(rng: np.random.Generator, k: int, n: int = 600, step: float = 0.1, error: float = 0.01,
                 dwell: int = 12) -> tuple[Series, np.ndarray]:
    """A series keyed between k levels step apart, holding each randomly chosen level for `dwell` points."""
    states = np.repeat(rng.integers(0, k, n // dwell + 1), dwell)[:n]
    err = np.full(n, error)
    time = np.sort(rng.uniform(0, 100, n))
    return Series(time, states * step + err * rng.standard_normal(n), err, kind="flux"), states


def evenly(value: np.ndarray, error: float) -> Series:
    n = value.size
    return Series(np.arange(n) * 0.1, value, np.full(n, error), kind="flux")


# ---------------------------------------------------------------------------------------------------------------
# benchmark
# ---------------------------------------------------------------------------------------------------------------

def test_eclipsing_binary_two_levels_matching_truth(zoo) -> None:
    item = zoo["natural/eclipsing_binary"]
    result = discreteness_test(item.payload, rng_for(SEED, "eclipse"))
    assert isinstance(result, LevelsResult)
    assert result.detected and result.k == 2
    assert result.levels == pytest.approx([1.0 - item.truth["depth"], 1.0], abs=5e-4)
    assert result.p_value <= 0.01
    assert result.symbols is not None
    agree = np.mean((result.symbols.symbols == 0) == item.truth["in_eclipse"])
    assert agree >= 0.99
    assert result.occupancy[0] == pytest.approx(item.truth["in_eclipse"].mean(), abs=0.005)
    assert result.transitions == pytest.approx(2 * 24, abs=4)        # about 24 eclipses in 90 days


@pytest.mark.parametrize("name", ["natural/drw_quasar", "natural/variable_star_p180", "natural/drw_plus_weak_sine",
                                  "control/white_noise_series"])
def test_continuous_benchmarks_are_not_discrete(zoo, name: str) -> None:
    result = discreteness_test(zoo[name].payload, rng_for(SEED, name))
    assert result.k == 1 and not result.detected
    assert result.symbols is None
    assert result.p_value == 1.0


def test_white_noise_sits_on_one_level(zoo) -> None:
    series = zoo["control/white_noise_series"].payload
    result = find_levels(series.value, series.error)
    assert result.k == 1 and result.off_level_fraction < 0.02
    assert "one level" in result.note
    assert np.isnan(result.p_value)


# ---------------------------------------------------------------------------------------------------------------
# discrete series
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("k", [2, 3, 4])
def test_keyed_series_found_with_right_k(k: int) -> None:
    rng = rng_for(SEED, "keyed", k)
    series, states = keyed_series(rng, k)
    result = discreteness_test(series, rng)
    assert result.detected and result.k == k
    assert result.levels == pytest.approx(np.arange(k) * 0.1, abs=0.005)
    assert np.mean(result.symbols.symbols == states) >= 0.99
    assert result.separation == pytest.approx(10.0, abs=0.5)


def test_rare_state_is_kept() -> None:
    # a low state occupying 2% of the time, 10 error bars below the main level
    rng = rng_for(SEED, "rare")
    n = 3000
    states = np.zeros(n, dtype=int)
    for start in rng.choice(n - 10, 6, replace=False):
        states[start:start + 10] = 1
    err = np.full(n, 0.002)
    series = evenly(1.0 - 0.02 * states + err * rng.standard_normal(n), 0.002)
    result = find_levels(series.value, series.error)
    assert result.k == 2 and result.detected
    assert result.occupancy[0] == pytest.approx(states.mean(), abs=0.003)


def test_heteroscedastic_errors() -> None:
    rng = rng_for(SEED, "hetero")
    n = 800
    states = np.repeat(rng.integers(0, 2, n // 20), 20)
    err = rng.uniform(0.005, 0.015, n)
    series = Series(np.arange(n) * 1.0, 0.2 * states + err * rng.standard_normal(n), err)
    result = discreteness_test(series, rng, null="smooth", n_decoys=100)
    assert result.detected and result.k == 2
    assert np.mean(result.symbols.symbols == states) >= 0.995


# ---------------------------------------------------------------------------------------------------------------
# the traps
# ---------------------------------------------------------------------------------------------------------------

def test_sinusoid_with_small_errors_is_not_discrete() -> None:
    # its values pile up at the two extremes: a naive two-component mixture calls this two-level
    rng = rng_for(SEED, "sine")
    t = np.sort(rng.uniform(0, 500, 800))
    err = np.full(t.size, 0.002)
    series = Series(t, 0.2 * np.sin(2 * np.pi * t / 37.3) + err * rng.standard_normal(t.size), err)
    result = discreteness_test(series, rng)
    assert result.k == 1 and not result.detected
    assert result.off_level_fraction > 0.5
    assert "continuous" in result.note


def test_sawtooth_is_not_discrete() -> None:
    rng = rng_for(SEED, "sawtooth")
    t = np.sort(rng.uniform(0, 500, 800))
    err = np.full(t.size, 0.002)
    series = Series(t, 0.2 * ((t / 37.3) % 1.0) + err * rng.standard_normal(t.size), err)
    assert not discreteness_test(series, rng).detected


def test_red_noise_with_tiny_errors_is_not_discrete() -> None:
    rng = rng_for(SEED, "red")
    t = np.arange(600) * 1.0
    err = np.full(t.size, 1e-4)
    series = Series(t, simulate_drw(t, 0.2, 100.0, rng) + err * rng.standard_normal(t.size), err)
    result = discreteness_test(series, rng)
    assert result.k == 1 and not result.detected


def test_short_span_red_noise_lumps_are_not_levels() -> None:
    # a slow random walk seen over a short span with moderate errors looks lumpy
    rng = rng_for(SEED, "lumpy")
    t = np.arange(400) * 1.0
    err = np.full(t.size, 0.01)
    series = Series(t, simulate_drw(t, 0.1, 5000.0, rng) + err * rng.standard_normal(t.size), err)
    assert not discreteness_test(series, rng).detected


def test_drw_null_rejects_levels_made_by_sampling() -> None:
    # red noise seen in five short visits: each visit is flat to within its errors and the visits are far apart
    # in value, so the strict definition alone is fooled about half the time; the red-noise decoys at the same
    # times are fooled just as often, so nothing is detected
    qualified = detected = 0
    for i in range(40):
        rng = rng_for(SEED, "visits", i)
        t = np.concatenate([e * 90.0 + np.sort(rng.uniform(0, 0.1, 10)) for e in range(5)])
        err = np.full(t.size, 0.005)
        series = Series(t, simulate_drw(t, 0.2, 200.0, rng) + err * rng.standard_normal(t.size), err)
        result = discreteness_test(series, rng, n_decoys=100)
        qualified += result.k >= 2
        detected += result.detected
    assert qualified >= 10
    assert detected == 0


def test_two_levels_four_error_bars_apart_are_unresolved() -> None:
    rng = rng_for(SEED, "unresolved")
    series, _ = keyed_series(rng, 2, step=0.04, error=0.01)
    result = discreteness_test(series, rng)
    assert result.k == 1 and not result.detected
    assert "unresolved" in result.note


def test_separation_threshold_is_respected() -> None:
    rng = rng_for(SEED, "threshold")
    series, _ = keyed_series(rng, 2, step=0.06, error=0.01)     # 6 error bars apart
    assert find_levels(series.value, series.error).k == 1
    assert find_levels(series.value, series.error, min_separation=5.0).k == 2


# ---------------------------------------------------------------------------------------------------------------
# symbols
# ---------------------------------------------------------------------------------------------------------------

def test_symbolise_probabilities() -> None:
    value = np.array([0.0, 0.1, 0.05, 0.001, 0.099])
    series = evenly(value, 0.01)
    symbols = symbolise(series, np.array([0.1, 0.0]))
    assert isinstance(symbols, SymbolSequence)
    assert symbols.alphabet_size == 2
    assert list(symbols.symbols[[0, 1, 3, 4]]) == [0, 1, 0, 1]
    assert np.allclose(symbols.probabilities.sum(axis=1), 1.0)
    assert 0.2 < symbols.probabilities[2, 0] < 0.8                   # half way between: uncertain
    assert symbols.probabilities[0, 0] > 0.999
    assert symbols.meta["levels"] == [0.0, 0.1]
    assert symbols.meta["off_level"] == 1


def test_quantise_is_labelled_lossy() -> None:
    rng = rng_for(SEED, "quantise")
    t = np.arange(1000) * 1.0
    series = Series(t, simulate_drw(t, 0.2, 50.0, rng), np.full(t.size, 0.01))
    symbols = quantise(series, 4)
    assert symbols.meta["lossy"] is True
    assert np.bincount(symbols.symbols, minlength=4) == pytest.approx([250] * 4, abs=2)
    assert np.allclose(symbols.probabilities.sum(axis=1), 1.0)
    uniform = quantise(series, 3, method="uniform")
    assert uniform.alphabet_size == 3 and len(uniform.meta["edges"]) == 2


def test_event_symbols_on_pulse_position_code_and_poisson() -> None:
    rng = rng_for(SEED, "event symbols")
    code = rng.integers(0, 4, 500)
    events = EventList(np.cumsum(1.0 + code + 0.005 * rng.standard_normal(500)))
    symbols = event_symbols(events)
    assert symbols is not None and symbols.alphabet_size == 4
    assert np.array_equal(symbols.symbols, code[1:])
    assert symbols.meta["levels"] == pytest.approx([1, 2, 3, 4], rel=0.01)
    assert event_symbols(EventList(np.cumsum(rng.exponential(1.0, 500)))) is None
    assert event_symbols(EventList(np.cumsum(rng.weibull(0.6, 500)))) is None


# ---------------------------------------------------------------------------------------------------------------
# calibration
# ---------------------------------------------------------------------------------------------------------------

def test_drw_survey_light_curves_are_rarely_discrete() -> None:
    detected = 0
    for i in range(200):
        rng = rng_for(SEED, "drw calibration", i)
        t = survey_times(rng)
        sigma = rng.uniform(0.05, 0.3)
        tau = float(np.exp(rng.uniform(np.log(20.0), np.log(3000.0))))
        err = rng.uniform(0.003, 0.01, t.size)
        series = Series(t, 17.0 + simulate_drw(t, sigma, tau, rng) + err * rng.standard_normal(t.size), err)
        detected += discreteness_test(series, rng, alpha=0.01).detected
    print(f"DRW calibration: {detected} of 200 detected at alpha = 0.01")
    assert detected <= 6
