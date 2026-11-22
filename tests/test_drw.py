"""The damped-random-walk engine: likelihood, simulation, whitening and least squares."""
import numpy as np
import pytest

from stargazer.core import Series
from stargazer.nulls import drw


def _sampling(rng, n=90, span=900.0):
    time = np.sort(rng.uniform(0, span, n))
    return time, rng.uniform(0.02, 0.05, n)


def _dense_log_like(time, error, value, sigma, tau):
    """Reference: the same marginal likelihood from the full covariance matrix."""
    cov = sigma ** 2 * np.exp(-np.abs(time[:, None] - time[None, :]) / tau) + np.diag(error ** 2)
    inv = np.linalg.inv(cov)
    ones = np.ones_like(time)
    a, b, c = ones @ inv @ ones, ones @ inv @ value, value @ inv @ value
    _, logdet = np.linalg.slogdet(cov)
    return -0.5 * (c - b * b / a) - 0.5 * logdet - 0.5 * np.log(a) - 0.5 * (time.size - 1) * np.log(2 * np.pi)


def test_likelihood_matches_dense_reference():
    rng = np.random.default_rng(1)
    time, error = _sampling(rng)
    value = drw.simulate(time, error, 0.2, 120.0, rng)[0]
    for sigma, tau in ((0.2, 120.0), (0.05, 10.0), (0.6, 4000.0)):
        fast, _ = drw.log_likelihood_grid(time, error, value, np.array([sigma]), np.array([tau]))
        assert fast[0] == pytest.approx(_dense_log_like(time, error, value, sigma, tau), abs=1e-6)


def test_per_row_parameters_match_shared_parameters():
    rng = np.random.default_rng(2)
    time, error = _sampling(rng)
    values = drw.simulate(time, error, np.full(3, 0.2), np.full(3, 100.0), rng)
    sigma, tau = np.array([0.1, 0.2, 0.3]), np.array([30.0, 100.0, 300.0])
    shared, _ = drw.log_likelihood_grid(time, error, values, sigma, tau)
    per_row, _ = drw.log_likelihood_grid(time, error, values, np.tile(sigma, (3, 1)), np.tile(tau, (3, 1)))
    assert np.allclose(shared, per_row)


def test_simulation_has_the_requested_variance_and_memory():
    rng = np.random.default_rng(3)
    time = np.arange(0.0, 400.0, 1.0)
    sims = drw.simulate(time, np.full(time.size, 1e-3), np.full(400, 0.3), np.full(400, 50.0), rng)
    assert sims.var(axis=0).mean() == pytest.approx(0.09, rel=0.08)
    lag10 = np.mean(sims[:, :-10] * sims[:, 10:])
    assert lag10 == pytest.approx(0.09 * np.exp(-10 / 50.0), rel=0.12)


def test_whitening_gives_unit_normal_innovations():
    rng = np.random.default_rng(4)
    time, error = _sampling(rng, n=600, span=3000.0)
    value = drw.simulate(time, error, 0.2, 150.0, rng)[0]
    z, _ = drw.whiten(Series(time, value, error), 0.2, 150.0)
    assert abs(z.mean()) < 0.15 and z.std() == pytest.approx(1.0, abs=0.1)


def test_white_noise_limit_is_lomb_scargle():
    from astropy.timeseries import LombScargle

    rng = np.random.default_rng(5)
    time, error = _sampling(rng, n=120)
    value = 0.1 * np.sin(2 * np.pi * time / 70.0) + error * rng.standard_normal(time.size)
    freq = np.linspace(0.005, 0.05, 40)
    stat, _, _ = drw.sinusoid_scan(time, error, value, freq, 1e-9, 1.0)
    reference = LombScargle(time, value, error).power(freq, normalization="psd")
    assert np.allclose(stat[0], 2 * reference, rtol=1e-4, atol=1e-6)


def test_least_squares_agrees_with_the_template_scan():
    rng = np.random.default_rng(6)
    time, error = _sampling(rng, n=150)
    value = drw.simulate(time, error, 0.2, 100.0, rng)[0]
    freq = 1 / 80.0
    stat, a, b = drw.sinusoid_scan(time, error, value, np.array([freq]), 0.2, 100.0)
    design = np.stack([np.ones_like(time), np.sin(2 * np.pi * freq * time), np.cos(2 * np.pi * freq * time)], 1)[None]
    coef, reduction = drw.gls(time, error, value[None], design, np.array([0.2]), np.array([100.0]))
    _, level_only = drw.gls(time, error, value[None], design[:, :, :1], np.array([0.2]), np.array([100.0]))
    assert coef[0, 1] == pytest.approx(a[0, 0], rel=1e-6) and coef[0, 2] == pytest.approx(b[0, 0], rel=1e-6)
    assert reduction[0] - level_only[0] == pytest.approx(stat[0, 0], rel=1e-6)


def test_refined_maximum_is_close_to_a_fine_grid():
    rng = np.random.default_rng(7)
    time, error = _sampling(rng, n=200, span=2000.0)
    values = drw.simulate(time, error, np.full(12, 0.15), np.full(12, 200.0), rng)
    coarse = drw.default_grids(time, values[0], error, n_sigma=16, n_tau=24)
    fine = drw.default_grids(time, values[0], error, n_sigma=80, n_tau=110)
    best_coarse, sigma, tau = drw.max_log_like(time, error, values, *coarse)
    best_fine, _, _ = drw.max_log_like(time, error, values, *fine, zooms=3)
    assert np.max(np.abs(best_fine - best_coarse)) < 0.3
    direct, _ = drw.log_likelihood_grid(time, error, values, sigma[:, None], tau[:, None])
    assert np.max(np.abs(direct[:, 0] - best_coarse)) < 0.1    # the value reported belongs to the point reported


def test_forecast_reverts_to_the_mean():
    rng = np.random.default_rng(8)
    time, error = _sampling(rng, n=100)
    series = Series(time, 5.0 + drw.simulate(time, error, 0.2, 100.0, rng)[0], error)
    mean, sd = drw.forecast(series, 0.2, 100.0, time[-1] + np.array([1.0, 5000.0]))
    assert sd[0] < sd[1] and sd[1] == pytest.approx(0.2, rel=0.02)
    assert abs(mean[1] - 5.0) < 0.2
