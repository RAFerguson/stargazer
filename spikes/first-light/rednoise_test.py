#!/usr/bin/env python3
"""Is the strongest periodogram peak in a quasar light curve significant against red noise?

Quasars flicker with correlated ("red") noise, which piles power into the lowest
frequencies of any periodogram. A test that assumes white noise therefore reports
absurdly small false-alarm odds for an ordinary quasar. This script calibrates the peak
against the right null instead:

1. Fit a damped random walk (DRW) plus the reported photometric errors to the light
   curve, on a grid, keeping the whole likelihood surface rather than one best fit.
2. Draw DRW parameters from that surface and simulate light curves at the real
   observation times with the real error bars.
3. Run the same Lomb-Scargle search on every simulation and ask how often the highest
   simulated peak beats the observed one. Taking the highest peak over the whole
   frequency grid builds the look-elsewhere correction into the answer.

    python rednoise_test.py data/pg1302_ztf_dr24.json --claimed-period 1884
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from astropy.timeseries import LombScargle
from scipy.linalg import cho_factor, cho_solve

SIGMA_GRID = np.exp(np.linspace(np.log(0.02), np.log(1.0), 40))   # DRW amplitude, magnitudes
TAU_GRID = np.exp(np.linspace(np.log(5.0), np.log(2.0e4), 50))    # DRW timescale, days
MIN_PERIOD_DAYS = 30.0
OVERSAMPLING = 10


def nightly_means(t, mag, err):
    """Inverse-variance weighted mean per night."""
    night = np.floor(t).astype(int)
    rows = []
    for n in np.unique(night):
        pick = night == n
        weight = 1.0 / err[pick] ** 2
        rows.append((np.sum(weight * t[pick]) / weight.sum(), np.sum(weight * mag[pick]) / weight.sum(),
                     1.0 / np.sqrt(weight.sum())))
    return (np.array(column) for column in zip(*rows))


def covariance(t, err, sigma, tau):
    lag = np.abs(t[:, None] - t[None, :])
    return sigma ** 2 * np.exp(-lag / tau) + np.diag(err ** 2)


def drw_log_likelihood(t, mag, err, sigma, tau):
    """Gaussian log-likelihood of a DRW, with the mean magnitude profiled out."""
    try:
        factor = cho_factor(covariance(t, err, sigma, tau), lower=True)
    except np.linalg.LinAlgError:
        return -np.inf
    ones = np.ones_like(mag)
    mean = (ones @ cho_solve(factor, mag)) / (ones @ cho_solve(factor, ones))
    residual = mag - mean
    return -0.5 * residual @ cho_solve(factor, residual) - np.sum(np.log(np.diag(factor[0])))


def analyse(t, mag, err, claimed_period, simulations, rng):
    t, mag, err = nightly_means(t, mag, err)
    span = t.max() - t.min()
    freq = np.arange(1.0 / span, 1.0 / MIN_PERIOD_DAYS, 1.0 / (OVERSAMPLING * span))
    periodogram = LombScargle(t, mag, err)
    power = periodogram.power(freq)
    peak = int(np.argmax(power))
    white_noise_fap = float(periodogram.false_alarm_probability(
        power[peak], method="baluev", minimum_frequency=freq[0], maximum_frequency=freq[-1]))

    surface = np.array([[drw_log_likelihood(t, mag, err, s, tau) for tau in TAU_GRID] for s in SIGMA_GRID])
    weight = np.exp(surface - surface.max())
    weight /= weight.sum()
    best_sigma, best_tau = np.unravel_index(np.argmax(surface), surface.shape)

    at_claim = int(np.argmin(np.abs(freq - 1.0 / claimed_period))) if claimed_period else peak
    simulated = np.empty((simulations, len(freq)))
    for row, draw in enumerate(rng.choice(weight.size, size=simulations, p=weight.ravel())):
        i, j = np.unravel_index(draw, surface.shape)
        fake = np.linalg.cholesky(covariance(t, err, SIGMA_GRID[i], TAU_GRID[j])) @ rng.standard_normal(len(t))
        simulated[row] = LombScargle(t, fake, err).power(freq)
    highest = simulated.max(axis=1)

    result = {
        "nights": int(len(t)),
        "span_days": float(span),
        "first_mjd": float(t.min()),
        "last_mjd": float(t.max()),
        "peak_period_days": float(1.0 / freq[peak]),
        "peak_power": float(power[peak]),
        "cycles_covered_at_peak": float(span * freq[peak]),
        "white_noise_fap": white_noise_fap,
        "red_noise_fap": float(np.mean(highest >= power[peak])),
        "drw_best_sigma_mag": float(SIGMA_GRID[best_sigma]),
        "drw_best_tau_days": float(TAU_GRID[best_tau]),
        "claimed_period_days": claimed_period,
        "power_at_claimed": float(power[at_claim]),
        "red_noise_fap_at_claimed": float(np.mean(simulated[:, at_claim] >= power[at_claim])),
        "simulations": simulations,
    }
    curves = {"t": t.tolist(), "mag": mag.tolist(), "err": err.tolist(), "freq": freq.tolist(),
              "power": power.tolist(), "null_p99": np.percentile(simulated, 99, axis=0).tolist()}
    return result, curves


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("lightcurve", type=pathlib.Path)
    parser.add_argument("--claimed-period", type=float, default=None, help="published period to test, days")
    parser.add_argument("--simulations", type=int, default=2000)
    parser.add_argument("--min-epochs", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261114)
    args = parser.parse_args()

    record = json.loads(args.lightcurve.read_text())
    rng = np.random.default_rng(args.seed)
    output = {"name": record["name"], "source": record["source"], "seed": args.seed, "bands": {}}
    for curve in record["curves"]:
        t, mag, err, flags = (np.asarray(curve[key], dtype=float) for key in ("hmjd", "mag", "magerr", "catflags"))
        good = flags == 0
        if good.sum() < args.min_epochs:
            continue
        result, curves = analyse(t[good], mag[good], err[good], args.claimed_period, args.simulations, rng)
        output["bands"][curve["band"]] = {"result": result, "curves": curves}
        print(f"{curve['band']} band: {result['nights']} nights over {result['span_days']:.0f} d; "
              f"peak at {result['peak_period_days']:.0f} d ({result['cycles_covered_at_peak']:.1f} cycles); "
              f"white-noise FAP {result['white_noise_fap']:.1e}; red-noise FAP {result['red_noise_fap']:.3f}; "
              f"red-noise FAP at claimed period {result['red_noise_fap_at_claimed']:.3f}")

    target = args.lightcurve.with_name(args.lightcurve.stem.replace("_ztf_dr24", "") + "_rednoise.json")
    target.write_text(json.dumps(output))
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
