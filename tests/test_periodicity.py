"""The periodicity search: recovery of planted signals, behaviour on pure noise, and the supporting pieces."""
import numpy as np
import pytest
from scipy.stats import chi2

from stargazer.bench.natural import simulate_drw, survey_times
from stargazer.core import Series, rng_for
from stargazer.detect import periodicity as per
from stargazer.trust import inject

FAST = per.SearchConfig(max_decoys=150, coherence_templates=60)


def _quasar(seed, sigma=0.12, tau=200.0, error=0.03):
    rng = rng_for(seed, "test-quasar")
    time = np.unique(np.floor(np.concatenate([survey_times(rng), survey_times(rng)]))) + 0.3
    value = simulate_drw(time, sigma, tau, rng) + error * rng.standard_normal(time.size)
    return Series(time, value, np.full(time.size, error), "r", f"test:{seed}")


def test_frequency_grid_respects_cycles_and_minimum_period():
    grid = per.frequency_grid(2400.0, per.SearchConfig(min_cycles=3, min_period_days=20, oversample=5))
    assert grid[0] == pytest.approx(3 / 2400.0) and grid[-1] < 1 / 20.0
    assert np.allclose(np.diff(grid), 1 / (5 * 2400.0))
    assert per.frequency_grid(2400.0)[-1] < 1 / 40.0          # the default stays clear of the lunar month


def test_calendar_periods_are_never_candidates():
    from stargazer.detect import calendar

    freq = per.frequency_grid(2400.0)
    masked = calendar.mask(freq, 2400.0, 0.5, 3)
    for period in (365.24, 182.6, 730.5, 121.7):
        assert masked[np.argmin(np.abs(freq - 1 / period))]
    assert not masked[np.argmin(np.abs(freq - 1 / 230.0))] and masked.mean() < 0.25
    host = _quasar(31)
    yearly = inject.sinusoid(host, 365.24, 0.3, 0.1)
    assert abs(per.search(yearly, rng_for(9, "year"), FAST).period_days - 365.24) > 20


def test_no_model_harmonic_is_shorter_than_the_minimum_period():
    time = np.linspace(0, 1000, 50)[None, :]
    design = per.harmonic_design(time, np.array([1 / 100.0]), 3, max_freq=1 / 40.0)
    assert np.any(design[0, :, 1] != 0) and np.any(design[0, :, 2] != 0) and np.all(design[0, :, 3] == 0)
    assert np.all(design[0, :, 6] == 0)


def test_harmonic_score_is_the_chi_squared_tail():
    freq = np.arange(0.00125, 0.05, 1 / 12000.0)
    scan = np.random.default_rng(0).chisquare(2, size=(1, freq.size))
    score = per.harmonic_score(scan, freq, 3)
    step = freq[1] - freq[0]
    for j in (0, 40, 300, 550):
        terms = [scan[0, i] for i in (int(round((k * freq[j] - freq[0]) / step)) for k in (1, 2, 3)) if i < freq.size]
        assert score[0, j] == pytest.approx(-chi2.logsf(sum(terms), 2 * len(terms)), rel=1e-9)


def test_merge_bands_puts_bands_on_one_scale():
    a, b = _quasar(1), _quasar(2)
    b = b.with_values(20.0 + 3.0 * b.value)
    merged = per.merge_bands([a, b])
    assert len(merged) == len(a) + len(b) and np.all(np.diff(merged.time) >= 0)
    assert abs(np.median(merged.value)) < 0.5 and 0.5 < np.std(merged.value) < 2.0


@pytest.mark.parametrize("shape,amplitude,period", [("sinusoid", 0.15, 230.0), ("sawtooth", 0.25, 150.0),
                                                    ("flare", 0.35, 210.0)])
def test_planted_signals_are_found_at_the_right_period(shape, amplitude, period):
    host = _quasar(11)
    planted = {"sinusoid": inject.sinusoid, "sawtooth": inject.sawtooth,
               "flare": lambda s, p, a, ph: inject.periodic_flares(s, p, a, 0.15, ph)}[shape](host, period, amplitude, 0.3)
    result = per.search(planted, rng_for(5, shape), FAST)
    assert result.status == "searched"
    assert result.period_days == pytest.approx(period, rel=0.05)
    assert result.p <= 0.02 and result.decoy.exceedances == 0
    assert result.coherence_p <= 0.1
    assert result.statistic > result.first_pass * 0.5
    if shape != "sinusoid":
        assert result.harmonic_fraction > 0.05


def test_pure_red_noise_is_rarely_significant():
    p = [per.search(_quasar(100 + i), rng_for(6, i), FAST).p for i in range(16)]
    assert sum(value <= 0.05 for value in p) <= 4          # expected 0.8; five or more has probability < 1e-3
    assert min(p) > 1e-3


def test_white_noise_test_would_have_been_fooled():
    results = [per.search(_quasar(200 + i), rng_for(7, i), FAST) for i in range(6)]
    fooled = sum(r.classical["white_noise_fap"] < 1e-3 for r in results)
    assert fooled >= 5                                      # red noise looks periodic to a white-noise test
    assert all(r.p > 1e-3 for r in results)


def test_too_little_data_is_reported_not_searched():
    series = _quasar(3).select(np.arange(len(_quasar(3))) < 20)
    assert per.search(series, rng_for(8), FAST).status.startswith("too few points")


def test_tail_probability_extrapolates_sensibly():
    rng = np.random.default_rng(1)
    maxima = rng.chisquare(2, size=(2000, 100)).max(axis=1)
    p, info = per.tail_probability(40.0, maxima)
    exact = 1 - (1 - np.exp(-20.0)) ** 100
    assert 0.2 * exact < p < 5 * exact
    assert 50 < info["effective_frequencies"] < 200
