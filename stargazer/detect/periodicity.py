"""Periodicity in irregular, noisy, red-noise-dominated light curves.

The search has three parts, and the order matters.

1. **A statistic that knows the noise is red and does not assume a sine wave.** A first pass fits a sinusoid
   at every trial frequency by generalised least squares under the damped-random-walk covariance fitted to the
   series (`nulls.drw.sinusoid_scan`). The strongest peaks then get a full joint fit: a Fourier series with
   several harmonics, with the noise parameters fitted again once that periodic model has been removed. The
   statistic is twice the gain in maximised log-likelihood of "periodic model plus red noise" over "red noise
   alone". The second step matters because a noise model fitted to data that contain a periodic signal absorbs
   it: a sawtooth looks like a random walk, and the first pass alone is nearly blind to one.
2. **Calibration by decoys.** The same procedure, fits included, is run on series simulated from the fitted
   null at the real observation times with the real error bars. The reported probability is the fraction of
   decoys whose statistic matches or beats the real one. It already contains the look-elsewhere correction for
   the frequency grid and the cost of estimating the noise from the data.
3. **A second calibration that needs no noise simulation.** A clock keeps time; noise does not. The whole
   search is repeated on the real data with the time axis warped: stretched and squeezed at random, as a clock
   whose cycles wander in length would run. Warping destroys a strict period and leaves noise statistically as
   it was, so the warped searches show what this object's own data yield when no strict clock is present. The
   fraction of warped searches that match the real one is a probability that does not depend on having
   simulated the noise correctly.

Only periods short enough to repeat `min_cycles` times within the data are searched. A period the data covers
once is a trend, not a period.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from ..core import Series
from ..nulls import drw
from ..trust.decoys import MonteCarloP, sequential_p
from . import calendar


@dataclass(frozen=True)
class SearchConfig:
    """Everything that defines a periodicity search. Hash it and commit the hash before running."""

    min_cycles: float = 3.0          # a trial period must fit at least this many times in the data
    min_period_days: float = 40.0    # no trial period and no harmonic of the model is shorter than this
    oversample: float = 5.0          # frequency grid points per natural resolution element 1 / span
    bin_days: float = 1.0            # points closer than this are averaged first (0 disables)
    harmonics: int = 3               # the periodic model is a Fourier series with this many harmonics
    candidates: int = 2              # peaks of the first-pass periodogram that get the full joint fit
    harmonic_margin: float = 8.0     # a multiple of the best frequency is preferred if its fit is this close
    calendar_tolerance: float = 0.5  # no candidate within this many resolution elements of a calendar period
    batch: int = 50
    max_decoys: int = 2000
    stop_exceedances: int = 10
    coherence_templates: int = 200   # number of warped-time searches
    coherence_jitter: float = 0.2    # fractional scatter of cycle lengths in the warped clocks
    coherence_below_p: float = 0.05  # run the warped-time test only when the decoy probability is below this
    n_sigma: int = 16
    n_tau: int = 24
    min_points: int = 40


@dataclass
class PeriodicityResult:
    source_id: str
    band: str
    n_points: int
    span_days: float
    status: str                          # "searched" or the reason the search did not run
    n_frequencies: int = 0
    period_days: float = float("nan")
    statistic: float = float("nan")      # twice the gain in maximised log-likelihood from the periodic model
    first_pass: float = float("nan")     # the sinusoid-only statistic at the same period, under the null's noise
    amplitude: float = float("nan")      # of the fundamental, in the unit of the series passed in
    phase_cycles: float = float("nan")   # phase of the fundamental at time zero, in cycles
    harmonic_fraction: float = float("nan")   # share of the periodic model's power in harmonics above the first
    coefficients: list = field(default_factory=list)   # level, then sine terms, then cosine terms
    noise: dict = field(default_factory=dict)          # fitted noise without and with the periodic model
    cycles_covered: float = float("nan")
    decoy: MonteCarloP | None = None     # probability against the red-noise null, trials over frequency included
    coherence_p: float = float("nan")    # fraction of warped-time searches that do as well (no noise simulation)
    p_tail: float = float("nan")         # extrapolated from a fit to the decoys' distribution; indicative only
    tail: dict = field(default_factory=dict)
    null: dict = field(default_factory=dict)
    classical: dict = field(default_factory=dict)   # what a white-noise test would have said, for the record
    flags: list[str] = field(default_factory=list)

    @property
    def p(self) -> float:
        return self.decoy.p if self.decoy else float("nan")

    def to_dict(self) -> dict:
        out = asdict(self)
        out["decoy"] = self.decoy.to_dict() if self.decoy else None
        return out


def frequency_grid(span_days: float, config: SearchConfig = SearchConfig()) -> np.ndarray:
    """Trial frequencies (per day) from `min_cycles` cycles across the data up to 1 / min_period_days."""
    f_min = config.min_cycles / span_days
    f_max = 1.0 / config.min_period_days
    if f_max <= f_min:
        return np.zeros(0)
    return np.arange(f_min, f_max, 1.0 / (config.oversample * span_days))


def band_scale(band: Series) -> tuple[float, float]:
    """The centre (inverse-variance weighted mean) and robust scatter (1.4826 x the median absolute deviation)
    by which `merge_bands` standardises one band."""
    weight = 1.0 / band.error ** 2
    centre = float(np.sum(weight * band.value) / weight.sum())
    scatter = 1.4826 * float(np.median(np.abs(band.value - np.median(band.value))))
    scatter = scatter if scatter > 0 else float(np.std(band.value)) or 1.0
    return centre, scatter


def merge_bands(bands: list[Series]) -> Series:
    """One series from several bands of the same source.

    Each band has its weighted mean removed and is divided by its own robust scatter, then the bands are
    interleaved in time. This puts bands of different amplitude on one scale, so the search uses every
    measurement, and it keeps the null honest: a quasar's flicker is the same realisation in every band, so
    agreement between bands is not independent evidence, and decoys for the merged series are single
    realisations sampled at all the bands' times.
    """
    if not bands:
        raise ValueError("no bands to merge")
    times, values, errors, names = [], [], [], []
    for band in bands:
        centre, scatter = band_scale(band)
        times.append(band.time)
        values.append((band.value - centre) / scatter)
        errors.append(band.error / scatter)
        names.append(band.band)
    meta = {"merged_bands": names, "standardised": True}
    return Series(np.concatenate(times), np.concatenate(values), np.concatenate(errors), "+".join(names),
                  bands[0].source_id, "standardised", bands[0].time_standard, meta)


def _classical(series: Series, freq: np.ndarray) -> dict:
    """The floating-mean Lomb–Scargle peak and its white-noise false-alarm probability, for contrast only."""
    from astropy.timeseries import LombScargle

    ls = LombScargle(series.time, series.value, series.error)
    power = ls.power(freq)
    k = int(np.argmax(power))
    fap = float(ls.false_alarm_probability(power[k], method="baluev", minimum_frequency=freq[0],
                                           maximum_frequency=freq[-1]))
    return {"period_days": float(1.0 / freq[k]), "power": float(power[k]), "white_noise_fap": fap}


def harmonic_design(template_time: np.ndarray, freq: np.ndarray, harmonics: int,
                    max_freq: float | None = None) -> np.ndarray:
    """Design matrices of a Fourier series: a constant, then sin and cos of each harmonic.

    template_time has shape (R, N) and freq shape (R,): row r is a series of period 1 / freq[r] evaluated on
    its own time axis (the true one, or a warped one). Returns (R, N, 1 + 2 x harmonics). Harmonics above
    `max_freq` are left out (their columns are zero and their coefficients come out zero), so that no part of
    the model reaches periods shorter than the search allows.
    """
    order = np.arange(1, harmonics + 1)
    angle = 2.0 * np.pi * template_time[:, :, None] * freq[:, None, None] * order[None, None, :]
    keep = np.ones((freq.size, 1, harmonics)) if max_freq is None else \
        (order[None, None, :] * freq[:, None, None] <= max_freq * (1 + 1e-9)).astype(float)
    return np.concatenate([np.ones(template_time.shape + (1,)), np.sin(angle) * keep, np.cos(angle) * keep], axis=2)


def harmonic_score(scan: np.ndarray, freq: np.ndarray, harmonics: int) -> np.ndarray:
    """A score for each trial frequency that adds up the first-pass statistics at its harmonics.

    A pulse or a sawtooth puts its power into several harmonics, and under red noise a harmonic can stand
    higher above the noise than the fundamental does. Choosing candidates by the single highest peak would then
    pick a harmonic and report half or a third of the true period. For each frequency f this sums the
    first-pass statistic at f, 2f, ... (as far as the grid reaches) and converts the sum to minus the log of
    its chi-squared tail probability, so that sums of different numbers of terms can be compared. It is used
    only to choose which frequencies get the full joint fit.
    """
    scan = np.clip(np.atleast_2d(scan), 0.0, None)
    step = freq[1] - freq[0]
    total = scan.copy()
    terms = np.ones(freq.size, dtype=int)
    for k in range(2, harmonics + 1):
        index = np.rint((k * freq - freq[0]) / step).astype(int)
        ok = index < freq.size
        total[:, ok] += scan[:, index[ok]]
        terms[ok] += 1
    x = 0.5 * total
    series, term = np.ones_like(x), np.ones_like(x)
    score = np.empty_like(x)
    for k in range(1, harmonics + 1):
        mine = terms == k
        score[:, mine] = x[:, mine] - np.log(series[:, mine])
        term = term * x / k
        series = series + term
    return score


def _candidate_score(scan: np.ndarray, freq: np.ndarray, time: np.ndarray, config: SearchConfig) -> np.ndarray:
    """The harmonic-sum score, with calendar periods ruled out as candidates."""
    score = harmonic_score(scan, freq, config.harmonics)
    if config.calendar_tolerance > 0:
        score[:, calendar.mask(freq, float(time[-1] - time[0]), config.calendar_tolerance,
                               config.harmonics)] = -np.inf
    return score


def top_peaks(statistic: np.ndarray, count: int, exclusion: int) -> np.ndarray:
    """Indices of the `count` highest peaks in each row of `statistic` (M, F), no two within `exclusion`
    grid points of each other. Shape (M, count)."""
    work = np.array(statistic, dtype=float)
    grid = np.arange(work.shape[1])[None, :]
    found = []
    for _ in range(count):
        best = np.argmax(work, axis=1)
        found.append(best)
        work[np.abs(grid - best[:, None]) <= exclusion] = -np.inf
    return np.stack(found, axis=1)


def joint_statistic(time: np.ndarray, error: np.ndarray, values: np.ndarray, template_time: np.ndarray,
                    freq: np.ndarray, sigma0: np.ndarray, tau0: np.ndarray, log_like0: np.ndarray,
                    sigma_axis: np.ndarray, tau_axis: np.ndarray, harmonics: int, max_freq: float | None = None
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Twice the gain in maximised log-likelihood from a periodic model, for each row.

    For row r: fit the Fourier series of period 1 / freq[r] to values[r] by generalised least squares under
    the null's noise (sigma0[r], tau0[r]); remove it; maximise the red-noise likelihood of what is left over
    the noise grid; compare with the null's maximised log-likelihood log_like0[r]. Returns (statistic,
    coefficients, sigma, tau) where sigma and tau are the noise parameters fitted alongside the periodic model.
    """
    design = harmonic_design(template_time, freq, harmonics, max_freq)
    coef, _ = drw.gls(time, error, values, design, sigma0, tau0)
    residual = values - np.einsum("rnc,rc->rn", design[:, :, 1:], coef[:, 1:])
    log_like1, sigma1, tau1 = drw.max_log_like(time, error, residual, sigma_axis, tau_axis)
    return 2.0 * (log_like1 - log_like0), coef, sigma1, tau1


@dataclass
class Evaluation:
    """What the search statistic and its ingredients were for each of M series."""

    statistic: np.ndarray       # (M,) best joint statistic over the candidate peaks
    freq: np.ndarray            # (M,) frequency of that candidate
    coefficients: np.ndarray    # (M, 1 + 2 x harmonics)
    first_pass: np.ndarray      # (M,) sinusoid-only statistic at that frequency
    log_like0: np.ndarray       # (M,) maximised log-likelihood of red noise alone
    sigma0: np.ndarray
    tau0: np.ndarray
    sigma1: np.ndarray          # noise fitted alongside the periodic model
    tau1: np.ndarray


def evaluate(time: np.ndarray, error: np.ndarray, values: np.ndarray, freq: np.ndarray, sigma_axis: np.ndarray,
             tau_axis: np.ndarray, config: SearchConfig) -> Evaluation:
    """The complete statistic for each row of `values` (M, N): noise fit, first-pass periodogram, joint fit of
    the strongest peaks. Real data and decoys both go through here."""
    values = np.atleast_2d(np.asarray(values, float))
    log_like0, sigma0, tau0 = drw.max_log_like(time, error, values, sigma_axis, tau_axis)
    scan, _, _ = drw.sinusoid_scan(time, error, values, freq, sigma0, tau0)
    peaks = top_peaks(_candidate_score(scan, freq, time, config), config.candidates,
                      int(np.ceil(config.oversample)))
    m, c = peaks.shape
    rows = np.repeat(np.arange(m), c)
    stat, coef, sigma1, tau1 = joint_statistic(
        time, error, values[rows], np.broadcast_to(time, (rows.size, time.size)), freq[peaks.ravel()],
        sigma0[rows], tau0[rows], log_like0[rows], sigma_axis, tau_axis, config.harmonics, 1.0 / config.min_period_days)
    pick = np.argmax(stat.reshape(m, c), axis=1)
    chosen = np.arange(m) * c + pick
    index = peaks[np.arange(m), pick]
    return Evaluation(stat[chosen], freq[index], coef[chosen], scan[np.arange(m), index], log_like0, sigma0, tau0,
                      sigma1[chosen], tau1[chosen])


def resolve_period(series: Series, freq: np.ndarray, real: Evaluation, sigma_axis: np.ndarray,
                   tau_axis: np.ndarray, config: SearchConfig) -> tuple[float, np.ndarray, float, float]:
    """Choose between a period and its harmonics: the shortest period that explains the data.

    A Fourier series at twice the true period contains the true signal as its second harmonic, so it fits a
    plain sinusoid as well as the right model does. Starting from the best-fitting frequency, this moves to
    twice or three times that frequency whenever the model there fits nearly as well (its statistic is within
    `config.harmonic_margin` of the best), and repeats. A pulse or a sawtooth stays where it is, because its
    odd harmonics are lost at twice the frequency. Only the reported period changes; the search statistic and
    its calibration do not. Returns (frequency, coefficients, sigma, tau).
    """
    best = (float(real.freq[0]), real.coefficients[0], float(real.sigma1[0]), float(real.tau1[0]))
    moved = True
    while moved:
        moved = False
        for factor in (3, 2):
            trial = factor * best[0]
            if trial > freq[-1]:
                continue
            stat, coef, sigma1, tau1 = joint_statistic(
                series.time, series.error, series.value[None, :], series.time[None, :], np.array([trial]),
                real.sigma0[:1], real.tau0[:1], real.log_like0[:1], sigma_axis, tau_axis, config.harmonics, 1.0 / config.min_period_days)
            if stat[0] >= real.statistic[0] - config.harmonic_margin:
                best = (float(trial), coef[0], float(sigma1[0]), float(tau1[0]))
                moved = True
                break
    return best


def search(series: Series, rng: np.random.Generator, config: SearchConfig = SearchConfig()) -> PeriodicityResult:
    """Search one series for a period and calibrate the best one against red noise."""
    if config.bin_days > 0:
        series = series.binned(config.bin_days)
    result = PeriodicityResult(series.source_id, series.band, len(series), series.span, "searched")
    if len(series) < config.min_points:
        result.status = f"too few points ({len(series)} < {config.min_points})"
        return result
    freq = frequency_grid(series.span, config)
    if freq.size < 10:
        result.status = "data too short for the requested number of cycles"
        return result
    sigma_axis, tau_axis = drw.default_grids(series.time, series.value, series.error, n_sigma=config.n_sigma,
                                             n_tau=config.n_tau)
    fit = drw.fit(series, sigma_grid=sigma_axis, tau_grid=tau_axis)
    real = evaluate(series.time, series.error, series.value, freq, sigma_axis, tau_axis, config)
    observed = float(real.statistic[0])
    drawn: list[np.ndarray] = []

    def decoys(n: int) -> np.ndarray:
        sigma, tau = fit.draw(n, rng)
        fake = drw.simulate(series.time, series.error, sigma, tau, rng)
        drawn.append(evaluate(series.time, series.error, fake, freq, sigma_axis, tau_axis, config).statistic)
        return drawn[-1]

    k = config.harmonics
    best_freq, coef, sigma1, tau1 = resolve_period(series, freq, real, sigma_axis, tau_axis, config)
    if best_freq != real.freq[0]:
        result.flags.append(f"a model at {1.0 / real.freq[0]:.1f} days fitted best; the period reported is its "
                            f"harmonic at {1.0 / best_freq:.1f} days, which fits nearly as well")
    power = coef[1:1 + k] ** 2 + coef[1 + k:] ** 2
    result.n_frequencies = int(freq.size)
    result.period_days = float(1.0 / best_freq)
    result.statistic = observed
    result.first_pass = float(real.first_pass[0])
    result.amplitude = float(np.sqrt(power[0]))
    result.phase_cycles = float(np.arctan2(coef[1 + k], coef[1]) / (2.0 * np.pi) % 1.0)
    result.harmonic_fraction = float(power[1:].sum() / power.sum()) if power.sum() > 0 else float("nan")
    result.coefficients = [float(c) for c in coef]
    result.cycles_covered = float(series.span * best_freq)
    result.noise = {"sigma_null": float(real.sigma0[0]), "tau_null_days": float(real.tau0[0]),
                    "sigma_with_signal": float(sigma1), "tau_with_signal_days": float(tau1)}
    result.decoy = sequential_p(observed, decoys, batch=config.batch, max_draws=config.max_decoys,
                                stop_exceedances=config.stop_exceedances)
    result.null = fit.summary()
    result.classical = _classical(series, freq)
    maxima = np.concatenate(drawn)
    if maxima.size >= 200:
        result.p_tail, result.tail = tail_probability(observed, maxima)
    if not fit.tau_constrained:
        result.flags.append("red-noise timescale not constrained by this baseline")
    if best_freq in (freq[0], freq[-1]):
        result.flags.append("best period at the edge of the search range")
    if result.decoy.p <= config.coherence_below_p:
        result.coherence_p = coherence_test(series, freq, result.period_days, observed, real, sigma_axis, tau_axis,
                                            rng, config)
    return result


def tail_probability(observed: float, null_maxima: np.ndarray) -> tuple[float, dict]:
    """An extrapolated probability for statistics beyond the reach of the decoys drawn.

    The largest of F independent chi-squared(2) variables scaled by c exceeds s with probability
    1 - (1 - exp(-s / 2c))^F. F (the effective number of independent frequencies) and c are fitted to the
    decoys' maxima by maximum likelihood, and the formula is evaluated at the observed statistic.

    This is an extrapolation of a fitted curve, not a measured frequency, and the curve is not the true
    distribution of the joint statistic. Checked against 40,000 simulated nulls for one realistic sampling
    pattern, fits to 1,000 decoys overstated the probability every time: by a factor of 4 at a true 0.01, 10
    at 0.001 and 24 at 0.0001. It errs towards calling things less significant, which is the safe direction.
    It is used only where no simulated null matched the observed statistic, and never to report a probability
    larger than the simulation's own resolution (see `pipelines.quasar.probability_used`).
    """
    from scipy.optimize import minimize

    m = np.asarray(null_maxima, float)
    m = m[m > 0]

    def negative_log_like(theta: np.ndarray) -> float:
        c, f = np.exp(theta)
        u = np.exp(-m / (2 * c))
        return -float(np.sum(np.log(f) + (f - 1) * np.log1p(-u) - m / (2 * c) - np.log(2 * c)))

    best = minimize(negative_log_like, np.log([1.0, max(np.exp(np.median(m) / 2) * 0.7, 2.0)]), method="Nelder-Mead")
    c, f = np.exp(best.x)
    p = float(-np.expm1(f * np.log1p(-np.exp(-observed / (2 * c)))))
    return max(p, 1e-300), {"scale": float(c), "effective_frequencies": float(f), "n_decoys": int(m.size)}


def _jittered_phase(time: np.ndarray, period: float, jitter: float, rng: np.random.Generator,
                    n_templates: int) -> np.ndarray:
    """Phases (radians) at the observation times for templates whose cycles vary in length.

    Each template's successive cycle lengths are drawn independently around `period` with fractional scatter
    `jitter` (log-normal, mean-preserving), from a random starting point, so every template has the same number
    of cycles on average as the periodic one but does not keep time. Shape (N, n_templates).
    """
    start, stop = float(time[0]), float(time[-1])
    n_cycles = int(np.ceil((stop - start) / period * 1.8)) + 4
    lengths = period * np.exp(jitter * rng.standard_normal((n_templates, n_cycles)) - 0.5 * jitter ** 2)
    edges = start - period * rng.random((n_templates, 1)) + np.concatenate(
        [np.zeros((n_templates, 1)), np.cumsum(lengths, axis=1)], axis=1)
    phase = np.empty((time.size, n_templates))
    for k in range(n_templates):
        cycle = np.searchsorted(edges[k], time, side="right") - 1
        cycle = np.clip(cycle, 0, n_cycles - 1)
        phase[:, k] = 2.0 * np.pi * (cycle + (time - edges[k, cycle]) / lengths[k, cycle])
    return phase


def coherence_test(series: Series, freq: np.ndarray, period_days: float, observed: float, real: Evaluation,
                   sigma_axis: np.ndarray, tau_axis: np.ndarray, rng: np.random.Generator,
                   config: SearchConfig = SearchConfig(), chunk: int = 20) -> float:
    """The fraction of warped-time searches whose statistic matches the real search's (`observed`).

    Each warp replaces time by the reading of a clock whose cycles, nominally `period_days` long, vary in
    length by the fraction `config.coherence_jitter` from one to the next. The search is repeated with every
    periodic template built on warped time (first pass, then the joint fit of the strongest peaks), while the
    noise model stays on true time. Smooth-spectrum noise has no phase to lose, so its warped searches do about
    as well as the real one and the fraction is large. A strict clock loses its coherence under warping and the
    fraction is small.

    Unlike the decoy probability this uses no simulated noise: the comparison set is made from the object's own
    data. It therefore still means something when the noise model is wrong, provided the true noise spectrum
    is smooth. It cannot separate a strict clock from an oscillation that stays coherent for longer than the
    data last, and with few cycles in the data it has little power. The smallest value it can return is
    1 / (number of warps + 1).
    """
    n_warps = config.coherence_templates
    phase = _jittered_phase(series.time, period_days, config.coherence_jitter, rng, n_warps)   # (N, n_warps)
    warped_time = phase * period_days / (2.0 * np.pi)
    exclusion = int(np.ceil(config.oversample))
    best = np.empty(n_warps)
    value = series.value[None, :]
    for start in range(0, n_warps, chunk):
        part = warped_time[:, start:start + chunk]                              # (N, c)
        c = part.shape[1]
        angle = (2.0 * np.pi * part[:, :, None] * freq[None, None, :]).reshape(len(series), -1)
        scan, _, _ = drw.template_scan(series.time, series.error, value, np.sin(angle), np.cos(angle),
                                       real.sigma0[0], real.tau0[0])
        peaks = top_peaks(_candidate_score(scan[0].reshape(c, freq.size), freq, series.time, config),
                          config.candidates, exclusion)                                     # (c, candidates)
        rows = np.repeat(np.arange(c), peaks.shape[1])
        stat, _, _, _ = joint_statistic(
            series.time, series.error, np.broadcast_to(series.value, (rows.size, len(series))), part.T[rows],
            freq[peaks.ravel()], np.full(rows.size, real.sigma0[0]), np.full(rows.size, real.tau0[0]),
            np.full(rows.size, real.log_like0[0]), sigma_axis, tau_axis, config.harmonics, 1.0 / config.min_period_days)
        best[start:start + c] = stat.reshape(c, -1).max(axis=1)
    return float((1 + np.sum(best >= observed)) / (n_warps + 1))


def fit_sinusoid(series: Series, period_days: float, sigma: float, tau: float) -> dict:
    """Amplitude, phase and fit improvement of a sinusoid of fixed period under the red-noise covariance."""
    statistic, amp_s, amp_c = drw.sinusoid_scan(series.time, series.error, series.value,
                                                np.array([1.0 / period_days]), sigma, tau)
    return {"statistic": float(statistic[0, 0]), "amplitude": float(np.hypot(amp_s[0, 0], amp_c[0, 0])),
            "phase_cycles": float(np.arctan2(amp_c[0, 0], amp_s[0, 0]) / (2.0 * np.pi) % 1.0),
            "sin": float(amp_s[0, 0]), "cos": float(amp_c[0, 0])}
