"""The damped random walk: the default red-noise null model for quasars and other flickering sources.

A damped random walk (an Ornstein–Uhlenbeck process) wanders at random and forgets its past over a timescale
tau. Its stationary standard deviation is sigma, and two points dt apart are correlated by exp(-dt / tau).
Observed through measurement noise it is a linear-Gaussian state-space model, so its exact likelihood, its
simulation and its forecasts all cost O(N) per series with a Kalman filter, at any sampling.

Everything here is vectorised over a leading batch of parameter pairs or of simulated series, because the
searches built on it fit thousands of simulated series to calibrate one real one.

Conventions: time in days, values in the unit of the series, sigma in the same unit, tau in days.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core import Series

LOG_2PI = float(np.log(2.0 * np.pi))


def default_grids(time: np.ndarray, value: np.ndarray, error: np.ndarray, *, n_sigma: int = 24,
                  n_tau: int = 32) -> tuple[np.ndarray, np.ndarray]:
    """Log-spaced grids for sigma and tau that bracket anything the data can constrain.

    sigma runs from a fifth of the typical error to five times the scatter of the values; tau from half the
    median sampling interval to twenty times the length of the data. The same rule is applied to real and
    simulated series, which is what matters for calibration.
    """
    scatter = float(np.std(value)) if len(value) > 1 else float(np.median(error))
    floor = float(np.median(error))
    sigma = np.geomspace(max(0.2 * floor, 1e-6), max(5.0 * scatter, 2.0 * floor), n_sigma)
    span = float(time[-1] - time[0])
    cadence = float(np.median(np.diff(time))) if len(time) > 2 else 1.0
    tau = np.geomspace(max(0.5 * cadence, 1e-3), max(20.0 * span, 10.0), n_tau)
    return sigma, tau


def mesh(sigma_axis: np.ndarray, tau_axis: np.ndarray, span: float) -> tuple[np.ndarray, np.ndarray]:
    """The grid of (sigma, tau) pairs, flattened with sigma varying slowest.

    For timescales longer than the data, only the ratio sigma^2 / tau (how fast the walk diffuses) is measured,
    and the likelihood forms a ridge along sigma proportional to sqrt(tau). The grid follows the ridge: at
    timescale tau the sigma axis is multiplied by sqrt(max(1, tau / span)). This keeps the ridge on the grid
    and makes the surface close to separable along the two grid directions.
    """
    stretch = np.sqrt(np.maximum(1.0, np.asarray(tau_axis, float) / span))
    sigma = (np.asarray(sigma_axis, float)[:, None] * stretch[None, :]).ravel()
    tau = np.broadcast_to(np.asarray(tau_axis, float)[None, :], (len(sigma_axis), len(tau_axis))).ravel()
    return sigma, tau


def log_likelihood_grid(time: np.ndarray, error: np.ndarray, value: np.ndarray, sigma: np.ndarray,
                        tau: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Marginal log-likelihood of a damped random walk plus noise, with the mean integrated out.

    value has shape (N,) or (M, N): one series or a batch sharing the same times and errors. sigma and tau are
    arrays of equal shape: (G,), one parameter pair per entry, tried on every series; or (M, G), a separate
    set of G pairs for each series. Returns (log_like, mean), each of shape (G,) or (M, G); mean is the
    generalised-least-squares estimate of the constant level at each pair.
    """
    time = np.asarray(time, float)
    error = np.asarray(error, float)
    y = np.atleast_2d(np.asarray(value, float))            # (M, N)
    sigma = np.asarray(sigma, float)
    tau = np.asarray(tau, float)
    if sigma.ndim == 1:
        sigma, tau = sigma[None, :], tau[None, :]          # (1, G): shared by all series
    var = sigma ** 2
    n = time.size
    dt = np.concatenate([[np.inf], np.diff(time)])
    shape = (y.shape[0], var.shape[1])
    m_one = np.zeros_like(var)                              # filtered state for the constant column
    m_y = np.zeros(shape)                                   # filtered state for each series
    p = var.copy()
    acc_a = np.zeros_like(var)
    acc_d = np.zeros_like(var)
    acc_b = np.zeros(shape)
    acc_c = np.zeros(shape)
    for i in range(n):
        a = np.exp(-dt[i] / tau)
        p_pred = a * a * p + var * (1.0 - a * a)
        s = p_pred + error[i] ** 2
        gain = p_pred / s
        v_one = 1.0 - a * m_one
        v_y = y[:, i, None] - a * m_y
        acc_a += v_one * v_one / s
        acc_b += v_one * v_y / s
        acc_c += v_y * v_y / s
        acc_d += np.log(s)
        m_one = a * m_one + gain * v_one
        m_y = a * m_y + gain * v_y
        p = (1.0 - gain) * p_pred
    mean = acc_b / acc_a
    log_like = -0.5 * (acc_c - acc_b ** 2 / acc_a) - 0.5 * acc_d - 0.5 * np.log(acc_a) - 0.5 * (n - 1) * LOG_2PI
    if np.ndim(value) == 1:
        return log_like[0], mean[0]
    return log_like, mean


@dataclass
class DrwFit:
    """The likelihood surface of a damped random walk for one series, on a grid.

    weight is the posterior over the grid under priors flat in log sigma and log tau (times any extra prior
    given to `fit`). sigma_map and tau_map are the grid point of highest posterior.
    """

    sigma: np.ndarray        # (G,) flattened grid
    tau: np.ndarray          # (G,)
    log_like: np.ndarray     # (G,)
    weight: np.ndarray       # (G,) sums to one
    mean: np.ndarray         # (G,)
    sigma_map: float
    tau_map: float
    mean_map: float
    span: float
    n: int
    log_step: tuple[float, float] = (0.0, 0.0)   # grid spacing in log sigma and log tau

    def draw(self, n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        """n parameter pairs drawn from the posterior: a grid cell by its weight, then a point inside the cell.

        Spreading the draw over the cell matters when the data pin the parameters down more tightly than the
        grid spacing; without it every simulation would use the same one or two grid points.
        """
        index = rng.choice(self.weight.size, size=n, p=self.weight)
        tau = self.tau[index] * np.exp(self.log_step[1] * rng.uniform(-0.5, 0.5, n))
        ridge = np.sqrt(np.maximum(1.0, tau / self.span) / np.maximum(1.0, self.tau[index] / self.span))
        sigma = self.sigma[index] * np.exp(self.log_step[0] * rng.uniform(-0.5, 0.5, n)) * ridge
        return sigma, tau

    def interval(self, name: str, level: float = 0.68) -> tuple[float, float]:
        """Central credible interval of "sigma" or "tau" from the grid posterior."""
        values = getattr(self, name)
        order = np.argsort(values)
        cumulative = np.cumsum(self.weight[order])
        low = values[order][np.searchsorted(cumulative, (1 - level) / 2)]
        high = values[order][min(np.searchsorted(cumulative, 1 - (1 - level) / 2), values.size - 1)]
        return float(low), float(high)

    @property
    def tau_constrained(self) -> bool:
        """False when the timescale is not pinned down: the upper end of its 95% interval exceeds a tenth of
        the data's length, the regime in which published work finds the fitted timescale biased."""
        return self.interval("tau", 0.95)[1] < 0.1 * self.span

    def summary(self) -> dict:
        return {"sigma_map": self.sigma_map, "tau_map_days": self.tau_map, "mean_map": self.mean_map,
                "sigma_68": self.interval("sigma"), "tau_68_days": self.interval("tau"),
                "tau_constrained": self.tau_constrained, "n": self.n, "span_days": self.span}


def fit(series: Series, *, sigma_grid: np.ndarray | None = None, tau_grid: np.ndarray | None = None,
        log_prior=None) -> DrwFit:
    """Fit a damped random walk to one series on a grid and return the whole surface."""
    if sigma_grid is None or tau_grid is None:
        default_sigma, default_tau = default_grids(series.time, series.value, series.error)
        sigma_grid = default_sigma if sigma_grid is None else sigma_grid
        tau_grid = default_tau if tau_grid is None else tau_grid
    sg, tg = mesh(sigma_grid, tau_grid, series.span)
    log_like, mean = log_likelihood_grid(series.time, series.error, series.value, sg, tg)
    log_post = log_like + (log_prior(sg, tg) if log_prior is not None else 0.0)
    weight = np.exp(log_post - log_post.max())
    weight /= weight.sum()
    best = int(np.argmax(log_post))
    steps = (float(np.log(sigma_grid[1] / sigma_grid[0])) if len(sigma_grid) > 1 else 0.0,
             float(np.log(tau_grid[1] / tau_grid[0])) if len(tau_grid) > 1 else 0.0)
    return DrwFit(sg, tg, log_like, weight, mean, float(sg[best]), float(tg[best]), float(mean[best]),
                  series.span, len(series), steps)


def map_parameters(time: np.ndarray, error: np.ndarray, values: np.ndarray, sigma: np.ndarray,
                   tau: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Best grid point for each row of `values` (M, N): the same estimate `fit` reports, in batch."""
    log_like, _ = log_likelihood_grid(time, error, values, sigma, tau)
    best = np.argmax(np.atleast_2d(log_like), axis=1)
    return sigma[best], tau[best]


def simulate(time: np.ndarray, error: np.ndarray, sigma: np.ndarray, tau: np.ndarray, rng: np.random.Generator,
             *, mean: float | np.ndarray = 0.0, noise: bool = True) -> np.ndarray:
    """Exact simulations of a damped random walk at the given times, one row per (sigma, tau) pair.

    Measurement noise with the given errors is added unless noise is False. Returns shape (M, N).
    """
    time = np.asarray(time, float)
    sigma = np.atleast_1d(np.asarray(sigma, float))
    tau = np.atleast_1d(np.asarray(tau, float))
    m, n = sigma.size, time.size
    out = np.empty((m, n))
    x = sigma * rng.standard_normal(m)
    out[:, 0] = x
    for i, dt in enumerate(np.diff(time), start=1):
        a = np.exp(-dt / tau)
        x = a * x + sigma * np.sqrt(1.0 - a * a) * rng.standard_normal(m)
        out[:, i] = x
    if noise:
        out += np.asarray(error, float)[None, :] * rng.standard_normal((m, n))
    return out + np.asarray(mean, float).reshape(-1, 1) if np.ndim(mean) else out + mean


def whiten(series: Series, sigma: float, tau: float) -> tuple[np.ndarray, float]:
    """Standardised one-step-ahead innovations of the series under the model, and the fitted mean.

    If the model is right the innovations are independent standard normal numbers whatever the sampling,
    which makes them the natural input for content tests: any structure left in them is structure the red
    noise does not explain.
    """
    _, mean = log_likelihood_grid(series.time, series.error, series.value, np.array([sigma]), np.array([tau]))
    mean = float(np.atleast_1d(mean)[0])
    var = sigma ** 2
    dt = np.concatenate([[np.inf], np.diff(series.time)])
    m, p = 0.0, var
    z = np.empty(len(series))
    for i in range(len(series)):
        a = np.exp(-dt[i] / tau)
        p_pred = a * a * p + var * (1.0 - a * a)
        s = p_pred + series.error[i] ** 2
        v = (series.value[i] - mean) - a * m
        z[i] = v / np.sqrt(s)
        gain = p_pred / s
        m = a * m + gain * v
        p = (1.0 - gain) * p_pred
    return z, mean


def forecast(series: Series, sigma: float, tau: float, new_time: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean and standard deviation of the process at times after the last observation, given the data."""
    _, mean = log_likelihood_grid(series.time, series.error, series.value, np.array([sigma]), np.array([tau]))
    mean = float(np.atleast_1d(mean)[0])
    var = sigma ** 2
    dt = np.concatenate([[np.inf], np.diff(series.time)])
    m, p = 0.0, var
    for i in range(len(series)):
        a = np.exp(-dt[i] / tau)
        p_pred = a * a * p + var * (1.0 - a * a)
        s = p_pred + series.error[i] ** 2
        gain = p_pred / s
        m = a * m + gain * ((series.value[i] - mean) - a * m)
        p = (1.0 - gain) * p_pred
    ahead = np.asarray(new_time, float) - series.time[-1]
    if np.any(ahead < 0):
        raise ValueError("forecast times must not precede the last observation")
    a = np.exp(-ahead / tau)
    return mean + a * m, np.sqrt(a * a * p + var * (1.0 - a * a))


def sinusoid_scan(time: np.ndarray, error: np.ndarray, values: np.ndarray, freq: np.ndarray, sigma: np.ndarray,
                  tau: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The periodogram a sinusoid deserves when the noise is red.

    For each row of `values` (M, N), with its own (sigma, tau), and each trial frequency, fit
    level + A sin + B cos by generalised least squares under damped-random-walk-plus-noise covariance, and
    return how much twice the log-likelihood improves over the level alone. Under the null each entry is
    close to chi-squared with two degrees of freedom; in white noise this reduces to the floating-mean
    Lomb–Scargle periodogram.

    Returns (statistic, a, b), each of shape (M, F), where a and b are the fitted sine and cosine amplitudes.
    """
    phase = 2.0 * np.pi * np.outer(np.asarray(time, float), np.asarray(freq, float))    # (N, F)
    return template_scan(time, error, values, np.sin(phase), np.cos(phase), sigma, tau)


def template_scan(time: np.ndarray, error: np.ndarray, values: np.ndarray, sin_t: np.ndarray, cos_t: np.ndarray,
                  sigma: np.ndarray, tau: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`sinusoid_scan` for arbitrary pairs of template columns.

    sin_t and cos_t have shape (N, F): F pairs of templates evaluated at the observation times. Each pair is
    fitted, with a constant level, by generalised least squares under the red-noise covariance. Used with
    sines and cosines of a strictly periodic phase this is the periodogram; used with a phase whose cycles
    have been stretched and squeezed at random it gives the response of the same data to templates that are
    not periodic, which is what `detect.periodicity.coherence_test` compares against.
    """
    time = np.asarray(time, float)
    error = np.asarray(error, float)
    y = np.atleast_2d(np.asarray(values, float))            # (M, N)
    sigma = np.broadcast_to(np.atleast_1d(np.asarray(sigma, float)), (y.shape[0],))
    tau = np.broadcast_to(np.atleast_1d(np.asarray(tau, float)), (y.shape[0],))
    var = sigma ** 2
    m_rows, n = y.shape
    dt = np.concatenate([[np.inf], np.diff(time)])
    shape = (m_rows, sin_t.shape[1])
    m_one, m_y, p = np.zeros(m_rows), np.zeros(m_rows), var.copy()
    m_s, m_c = np.zeros(shape), np.zeros(shape)
    a11, b1 = np.zeros(m_rows), np.zeros(m_rows)
    a1s, a1c, ass, asc, acc, bs, bc = (np.zeros(shape) for _ in range(7))
    for i in range(n):
        a = np.exp(-dt[i] / tau)
        p_pred = a * a * p + var * (1.0 - a * a)
        s = p_pred + error[i] ** 2
        gain = p_pred / s
        w = 1.0 / s
        v_one = 1.0 - a * m_one
        v_y = y[:, i] - a * m_y
        v_s = sin_t[i][None, :] - a[:, None] * m_s
        v_c = cos_t[i][None, :] - a[:, None] * m_c
        a11 += w * v_one * v_one
        b1 += w * v_one * v_y
        wv_one, wv_y = (w * v_one)[:, None], (w * v_y)[:, None]
        a1s += wv_one * v_s
        a1c += wv_one * v_c
        ass += w[:, None] * v_s * v_s
        asc += w[:, None] * v_s * v_c
        acc += w[:, None] * v_c * v_c
        bs += wv_y * v_s
        bc += wv_y * v_c
        m_one = a * m_one + gain * v_one
        m_y = a * m_y + gain * v_y
        m_s = a[:, None] * m_s + gain[:, None] * v_s
        m_c = a[:, None] * m_c + gain[:, None] * v_c
        p = (1.0 - gain) * p_pred
    # Eliminate the level, then solve the two-by-two system for the sine and cosine amplitudes.
    inv = 1.0 / a11[:, None]
    tss = ass - a1s * a1s * inv
    tsc = asc - a1s * a1c * inv
    tcc = acc - a1c * a1c * inv
    ts = bs - a1s * b1[:, None] * inv
    tc = bc - a1c * b1[:, None] * inv
    det = tss * tcc - tsc * tsc
    det = np.where(det > 1e-12 * tss * tcc, det, np.inf)    # degenerate frequencies carry no evidence
    amp_s = (ts * tcc - tc * tsc) / det
    amp_c = (tc * tss - ts * tsc) / det
    statistic = ts * amp_s + tc * amp_c
    return statistic, amp_s, amp_c


def whiten_columns(time: np.ndarray, error: np.ndarray, columns: np.ndarray, sigma: np.ndarray,
                   tau: np.ndarray) -> np.ndarray:
    """Pass each column through the whitening filter of the noise model.

    columns has shape (R, N, C): for each of R rows, C series sampled at the N observation times; row r is
    whitened with its own (sigma[r], tau[r]). The result has the same shape. A column that is pure noise from
    the model comes out as independent standard normal numbers, so generalised least squares under the
    red-noise covariance becomes ordinary least squares on whitened columns.
    """
    time = np.asarray(time, float)
    error = np.asarray(error, float)
    columns = np.asarray(columns, float)
    rows, n, _ = columns.shape
    sigma = np.broadcast_to(np.atleast_1d(np.asarray(sigma, float)), (rows,))
    tau = np.broadcast_to(np.atleast_1d(np.asarray(tau, float)), (rows,))
    var = sigma ** 2
    dt = np.concatenate([[np.inf], np.diff(time)])
    state = np.zeros((rows, columns.shape[2]))
    p = var.copy()
    out = np.empty_like(columns)
    for i in range(n):
        a = np.exp(-dt[i] / tau)
        p_pred = a * a * p + var * (1.0 - a * a)
        s = p_pred + error[i] ** 2
        gain = p_pred / s
        v = columns[:, i, :] - a[:, None] * state
        out[:, i, :] = v / np.sqrt(s)[:, None]
        state = a[:, None] * state + gain[:, None] * v
        p = (1.0 - gain) * p_pred
    return out


def gls(time: np.ndarray, error: np.ndarray, values: np.ndarray, design: np.ndarray, sigma: np.ndarray,
        tau: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Generalised least squares of each row of `values` (R, N) on its own design matrix (R, N, C) under the
    red-noise covariance of that row. Returns (coefficients (R, C), reduction in chi-squared (R,))."""
    values = np.atleast_2d(np.asarray(values, float))
    white = whiten_columns(time, error, np.concatenate([design, values[:, :, None]], axis=2), sigma, tau)
    wx, wy = white[:, :, :-1], white[:, :, -1]
    normal = np.einsum("rnc,rnd->rcd", wx, wx)
    rhs = np.einsum("rnc,rn->rc", wx, wy)
    ridge = 1e-10 * np.einsum("rcc->r", normal)[:, None, None] * np.eye(design.shape[2])[None]
    coef = np.linalg.solve(normal + ridge, rhs[:, :, None])[:, :, 0]
    return coef, np.einsum("rc,rc->r", rhs, coef)


def max_log_like(time: np.ndarray, error: np.ndarray, values: np.ndarray, sigma_axis: np.ndarray,
                 tau_axis: np.ndarray, zooms: int = 2) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The maximum of the marginal log-likelihood over (sigma, tau) for each row of `values` (R, N).

    The likelihood is evaluated on the grid `mesh(sigma_axis, tau_axis, span)`. Around each row's best grid
    point a five-by-five grid of half the spacing is evaluated, `zooms` times over, and the last one is
    finished with a parabola along each direction. Returns (log-likelihood, sigma, tau), each of shape (R,).
    Real data and every simulation go through this same function.
    """
    values = np.atleast_2d(np.asarray(values, float))
    span = float(time[-1] - time[0])
    sg, tg = mesh(sigma_axis, tau_axis, span)
    surface, _ = log_likelihood_grid(time, error, values, sg, tg)
    rows = np.arange(values.shape[0])
    flat = np.argmax(surface, axis=1)
    i, j = np.unravel_index(flat, (sigma_axis.size, tau_axis.size))
    best = surface[rows, flat]
    u, w = np.log(sigma_axis[i]), np.log(tau_axis[j])       # position along the two grid directions
    step_u = np.log(sigma_axis[1] / sigma_axis[0]) if sigma_axis.size > 1 else 0.0
    step_w = np.log(tau_axis[1] / tau_axis[0]) if tau_axis.size > 1 else 0.0
    offsets = np.arange(-2, 3)
    local = None
    for _ in range(zooms):
        step_u, step_w = step_u / 2.0, step_w / 2.0
        grid_u = np.broadcast_to(u[:, None, None] + step_u * offsets[None, :, None], (rows.size, 5, 5))
        grid_w = np.broadcast_to(w[:, None, None] + step_w * offsets[None, None, :], (rows.size, 5, 5))
        tau = np.exp(grid_w).reshape(rows.size, 25)
        sigma = np.exp(grid_u).reshape(rows.size, 25) * np.sqrt(np.maximum(1.0, tau / span))
        local, _ = log_likelihood_grid(time, error, values, sigma, tau)
        local = local.reshape(rows.size, 5, 5)
        flat = np.argmax(local.reshape(rows.size, 25), axis=1)
        a, b = np.unravel_index(flat, (5, 5))
        best = local[rows, a, b]
        u, w = u + step_u * offsets[a], w + step_w * offsets[b]
    if local is not None:
        def parabola(index: np.ndarray, take) -> tuple[np.ndarray, np.ndarray]:
            inside = (index > 0) & (index < 4)
            y0, y2 = take(np.clip(index - 1, 0, 4)), take(np.clip(index + 1, 0, 4))
            curvature = y0 - 2.0 * best + y2
            ok = inside & (curvature < -1e-12)
            safe = np.where(ok, curvature, -1.0)
            shift = np.where(ok, np.clip(0.5 * (y0 - y2) / safe, -0.5, 0.5), 0.0)
            return shift, np.where(ok, np.clip(-(y0 - y2) ** 2 / (8.0 * safe), 0.0, None), 0.0)

        shift_u, gain_u = parabola(a, lambda k: local[rows, k, b])
        shift_w, gain_w = parabola(b, lambda k: local[rows, a, k])
        best = best + gain_u + gain_w
        u, w = u + shift_u * step_u, w + shift_w * step_w
    tau = np.exp(w)
    return best, np.exp(u) * np.sqrt(np.maximum(1.0, tau / span)), tau
