"""Hidden clocks and discrete spacings in the arrival times of events (pulses, bursts, photons).

Two questions are asked of an `EventList`, each answered by a statistic and a probability from surrogates that
destroy only the property being tested.

1. **Is there a clock?** (`clock_test`) Do the events fall on the ticks of one regular clock, with some ticks
   possibly empty? On-off keying, a pulsar and any slotted transmitter look like this. The statistic is the
   largest Rayleigh power over a grid of trial periods. The null keeps every event inside its own window of a
   fixed grid of windows, which keeps every change of rate slower than the window. By default the contents of
   each window are rotated as a block (`rotation_surrogates`), which also keeps clustering inside the window and
   destroys only phase coherence between windows; plain `jitter_surrogates`, which move each event
   independently, are available and are shown in `clock_test` to over-detect clustered sources.
2. **Are the spacings quantised?** (`spacing_test`) Are the intervals between consecutive events whole
   multiples of one quantum, as in pulse-position coding? The statistic is the largest phase coherence of the
   intervals over a grid of trial quanta. The null is a *smoothed bootstrap*: intervals drawn from a kernel
   density estimate of log-interval, which keeps the shape of the interval distribution and erases any comb.

`to_bits` then turns a clock into a bit stream of slot occupancy.

A clock is not a message. A pulsar is a clock and is correctly detected here; whether its on/off sequence
carries anything more than a two-state chain is a question for the content analysis downstream.

What these tests cannot see: a clock whose period drifts by more than about a tenth of a period over the data
(the statistic assumes strict periodicity; no period derivative is searched); periodic *modulation of the rate*
with many events per cycle (activity windows, orbital modulation), which the null deliberately keeps (see
`clock_test`); and information carried in amplitudes or in the order of intervals rather than in their values.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.fft as sfft

from ..core import BitStream, EventList
from ..trust.decoys import MonteCarloP, sequential_p

_MAX_FFT = 2 ** 22          # largest transform used; the frequency grid is coarsened rather than exceed it
_MIN_EVENTS = 20


# ---------------------------------------------------------------------------------------------------------------
# Rayleigh power
# ---------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class _Grid:
    """The binning and frequency grid shared by the data and its surrogates (so both see the same procedure)."""

    t0: float
    dt: float           # bin width
    length: int         # transform length; times are wrapped modulo length * dt, which is exact on the grid
    k_lo: int           # first and last frequency index searched
    k_hi: int

    @property
    def freq(self) -> np.ndarray:
        return np.arange(self.k_lo, self.k_hi + 1) / (self.length * self.dt)


def _make_grid(times: np.ndarray, min_period: float, max_period: float, oversample: float) -> _Grid:
    span = float(times[-1] - times[0])
    dt = min_period / 8.0
    length = sfft.next_fast_len(int(np.ceil(oversample * max(span, max_period) / dt)) + 1, real=True)
    if length > _MAX_FFT:   # keep memory and time bounded: coarsen the frequency spacing instead of the bins
        length = _MAX_FFT
    df = 1.0 / (length * dt)
    k_lo = max(1, int(np.ceil(1.0 / (max_period * df))))
    k_hi = min(length // 2, int(np.floor(1.0 / (min_period * df))))
    return _Grid(float(times[0]), dt, length, k_lo, max(k_hi, k_lo))


def _fft_power(times: np.ndarray, grid: _Grid) -> np.ndarray:
    """Rayleigh power on the grid for one set of times (N,) or a batch (M, N), via binned FFT."""
    t = np.atleast_2d(times)
    m, n = t.shape
    index = np.rint((t - grid.t0) / grid.dt).astype(np.int64) % grid.length
    out = np.empty((m, grid.k_hi - grid.k_lo + 1))
    for row in range(m):
        counts = np.bincount(index[row], minlength=grid.length).astype(np.float32)
        spec = sfft.rfft(counts)[grid.k_lo:grid.k_hi + 1]
        out[row] = (spec.real.astype(float) ** 2 + spec.imag.astype(float) ** 2) * 2.0 / n
    return out if np.ndim(times) == 2 else out[0]


def _exact_power(times: np.ndarray, freq: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact Rayleigh power |sum exp(2 pi i f t)|^2 * 2 / N and the phase angle of the sum, at each frequency."""
    t = times - times[0]
    total = np.zeros(np.size(freq), dtype=complex)
    for chunk in range(0, t.size, 4096):
        total += np.exp(2j * np.pi * np.outer(np.atleast_1d(freq), t[chunk:chunk + 4096])).sum(axis=1)
    return np.abs(total) ** 2 * 2.0 / t.size, np.angle(total)


def _refine(times: np.ndarray, f0: float, half_width: float, steps: int = 2) -> float:
    """Locally maximise the exact Rayleigh power around f0 by successively finer grids."""
    f = f0
    for _ in range(steps):
        trial = f + np.linspace(-half_width, half_width, 41)
        trial = trial[trial > 0]
        f = float(trial[int(np.argmax(_exact_power(times, trial)[0]))])
        half_width /= 10.0
    return f


def rayleigh_spectrum(times: np.ndarray, min_period: float, max_period: float,
                      oversample: float = 4) -> tuple[np.ndarray, np.ndarray]:
    """Rayleigh power of event times over trial periods from min_period to max_period.

    The Rayleigh power at frequency f is Z = 2 |sum_j exp(2 pi i f t_j)|^2 / N. For events at independent
    uniform random times it is close to chi-squared with two degrees of freedom (mean 2) at every frequency; a
    fraction R of events locked to one phase of period 1/f raises it to about 2 N R^2. It is computed by
    binning the times in bins of min_period / 8 and taking one FFT, zero-padded so that the frequency spacing is
    1 / (oversample * span). Binning to the nearest bin costs at most 3% of the amplitude at the shortest
    period. Returns (frequencies, power), frequencies ascending; periods are 1 / frequencies.
    """
    times = np.sort(np.asarray(times, float))
    grid = _make_grid(times, min_period, max_period, oversample)
    return grid.freq, _fft_power(times, grid)


# ---------------------------------------------------------------------------------------------------------------
# Surrogates
# ---------------------------------------------------------------------------------------------------------------

def jitter_surrogates(times: np.ndarray, width: float, rng: np.random.Generator, n: int) -> np.ndarray:
    """n surrogate event lists, shape (n, N), sorted along each row.

    The time axis is cut into windows of length `width` on a fixed grid whose origin is drawn afresh for each
    surrogate, and every event is moved to a uniform random place inside its own window. The number of events in
    every window is kept, so every change of rate slower than `width` survives; all timing finer than `width`,
    including any clock with a period well below it, is destroyed. A clock of period P keeps at most a fraction
    of about P / (pi * width) of its coherent amplitude.
    """
    times = np.asarray(times, float)
    origin = times[0] - width * rng.random((n, 1))
    window = np.floor((times[None, :] - origin) / width)
    out = origin + (window + rng.random((n, times.size))) * width
    out.sort(axis=1)
    return out


def rotation_surrogates(times: np.ndarray, width: float, rng: np.random.Generator, n: int) -> np.ndarray:
    """n surrogate event lists, shape (n, N), sorted along each row.

    The same fixed grid of windows as `jitter_surrogates` (random origin per surrogate), but instead of moving
    each event independently, the contents of each window are rotated cyclically by one common random amount,
    uniform in [0, width). Every event stays in its own window, so every change of rate slower than `width` is
    kept, exactly as with jitter; in addition every interval inside a window is kept (except the one pair that
    straddles the rotation point), so clustering, bursts and dead time on scales shorter than the window survive.
    What is destroyed is the phase relation between different windows, which is what a clock is: a clock of
    period P << width keeps only the coherence that the few events of one window carry.
    """
    times = np.asarray(times, float)
    origin = times[0] - width * rng.random((n, 1))
    window = np.floor((times[None, :] - origin) / width)
    n_windows = int(window.max()) + 1
    shift = width * rng.random((n, n_windows))
    rows = np.arange(n)[:, None]
    start = origin + window * width
    out = start + np.mod(times[None, :] - start + shift[rows, window.astype(np.int64)], width)
    out.sort(axis=1)
    return out


class _LiveTime:
    """Maps between real time and live time (time counted only while the instrument was looking).

    Surrogates are drawn in live time and mapped back, so that no surrogate event falls in an exposure gap and
    the sharp edges of the exposure, which carry power at all frequencies, are the same in data and null.
    """

    def __init__(self, exposure: list[tuple[float, float]]) -> None:
        spans = sorted((float(a), float(b)) for a, b in exposure if b > a)
        merged: list[list[float]] = []
        for a, b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        self.start = np.array([a for a, _ in merged])
        self.length = np.array([b - a for a, b in merged])
        self.cum = np.concatenate([[0.0], np.cumsum(self.length)[:-1]])

    def to_live(self, t: np.ndarray) -> np.ndarray:
        i = np.clip(np.searchsorted(self.start, t, side="right") - 1, 0, self.start.size - 1)
        return self.cum[i] + np.clip(t - self.start[i], 0.0, self.length[i])

    def to_real(self, u: np.ndarray) -> np.ndarray:
        i = np.clip(np.searchsorted(self.cum, u, side="right") - 1, 0, self.cum.size - 1)
        return self.start[i] + (u - self.cum[i])


# ---------------------------------------------------------------------------------------------------------------
# Clock test
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class ClockResult:
    """The outcome of `clock_test`.

    period, phase: the clock found; ticks at phase + k * period (phase is a time, in [0, period)).
    power: exact Rayleigh power at that period (about 2 for no clock, about 2 N for a perfect one).
    p_value: probability that the surrogates (window rotation by default, see `clock_test`) reach the observed
        maximum power over the trial periods; it already accounts for searching many periods.
    n_events, n_cycles: number of events, and span / period.
    occupancy: fraction of clock ticks within the data that carry at least one event.
    residual_rms: root-mean-square distance of the events from their nearest tick, in units of the period
        (0.289 for phases spread uniformly).
    jitter_width: the window width W of the surrogates (rotation or jitter).
    detected: p_value <= alpha.
    harmonics_note: which multiple of the best-power period was reported, and why.
    """

    period: float
    phase: float
    power: float
    p_value: MonteCarloP | float
    n_events: int
    n_cycles: float
    occupancy: float
    residual_rms: float
    jitter_width: float
    detected: bool
    harmonics_note: str
    min_period: float = float("nan")
    max_period: float = float("nan")


def _phase_stats(times: np.ndarray, period: float) -> tuple[float, float, float]:
    """(mean resultant length, phase as a time in [0, period), circular scatter in time units)."""
    angle = 2.0 * np.pi * (times - times[0]) / period
    total = np.exp(1j * angle).mean()
    r = float(np.abs(total))
    phase = (times[0] + period * np.angle(total) / (2.0 * np.pi)) % period
    scatter = period / (2.0 * np.pi) * np.sqrt(-2.0 * np.log(max(r, 1e-300)))
    return r, float(phase), float(scatter)


def default_period_range(events: EventList) -> tuple[float, float]:
    """The default trial periods: from a tenth of the median interval to twice it.

    A clock whose ticks carry at most one event has a period no longer than the typical spacing of events, and
    occupancies down to about a tenth are covered by going ten times shorter. Longer trial periods are a
    question about modulation of the rate, which this test does not ask (see `clock_test`). The lower end is
    raised to four times the median timing error when `time_error` is given.
    """
    d = np.diff(events.time)
    d = d[d > 0]
    med = float(np.median(d)) if d.size else 1.0
    floor = 4.0 * float(np.median(events.time_error)) if events.time_error is not None else 0.0
    lo = max(med / 10.0, floor)
    return lo, max(2.0 * med, 2.0 * lo)


def clock_test(events: EventList, rng: np.random.Generator, *, min_period: float | None = None,
               max_period: float | None = None, n_surrogates: int = 200, alpha: float = 1e-2,
               jitter_factor: float = 5.0, null: str = "rotate", oversample: float = 4) -> ClockResult:
    """Do the events fall on the ticks of a regular clock, with some ticks possibly empty?

    Statistic: the largest Rayleigh power (see `rayleigh_spectrum`) over trial periods min_period..max_period
    (by default a tenth of the median interval to twice it; see `default_period_range`).

    Null. The surrogates use windows of width W = jitter_factor * max_period, by default 5 * max_period, which
    is 10 median intervals. Two requirements fix W. It must be several times the longest trial period, so that
    a clock anywhere in the searched range is erased from the surrogates. And it must be short compared with the
    rate changes a natural source has, so that they survive into the null: tying W to the median interval means
    every change of rate that involves more than about ten events is kept.

    The default null is `rotation_surrogates` ("rotate"), not plain `jitter_surrogates` ("jitter"), and this is
    a deliberate departure. Independent jitter also erases clustering on scales shorter than W (a burst of a
    few events), and clustering raises the Rayleigh power smoothly over the whole band of short periods: with
    jitter, 19% of clustered (Weibull, shape 0.6) event lists gave p <= 0.05, and the repeating-burster
    benchmark was called a clock. Rotating each window's contents as a block keeps that clustering, keeps the
    rate in every window exactly, and destroys only the phase relation between windows, which is what a clock
    has and clustering lacks. The price is a little power: a clock keeps in the surrogates the coherence of the
    roughly ten events inside one window (a Rayleigh power of about 20 R^2 for a phase concentration R, against
    about 2 N R^2 in the data).

    When `events.exposure` is given, the surrogates are drawn in live time (time counted only during
    exposure) and mapped back, so no surrogate event falls in a gap and the exposure edges are the same in data
    and null.

    Because the default search stops at twice the median interval, a periodic modulation of the rate with many
    events per cycle (for example activity windows) is not searched for, and would be kept by the null: it is
    rate, not timing. Pass max_period to ask about longer clocks; the null then keeps only rate changes slower
    than 5 * max_period, so a strong modulation of the rate at those periods will be reported as a clock.

    p-value: the fraction of surrogates whose maximum power reaches the observed maximum, by Besag and
    Clifford's sequential scheme (`sequential_p`). Unremarkable lists stop after ten exceedances and cost a few
    dozen surrogates; remarkable ones use all n_surrogates, and then p = 1 / (n_surrogates + 1) at best. So a
    detection at alpha = 0.01 needs n_surrogates >= 99.

    Fundamental versus harmonics: events on the ticks of period P are equally concentrated in phase at P / 2,
    P / 3 ...; the largest power can fall on any of these. The reported period is the longest multiple k P' of
    the best-power period P' (up to max_period) at which the events' circular scatter about the ticks, in time
    units, is no more than 25% worse; at 2 P, events on alternate ticks fall half a period apart and the scatter
    jumps. This cannot tell a clock of period P whose odd ticks are always empty from a clock of period 2 P; it
    reports 2 P, the shorter description. Nor can it see a clock whose ticks are all more than ten times as far
    apart as the typical event spacing (several events per tick) unless min_period and max_period are given.
    """
    times = np.sort(np.asarray(events.time, float))
    n = times.size
    lo_default, hi_default = default_period_range(events)
    lo = float(min_period) if min_period is not None else lo_default
    hi = float(max_period) if max_period is not None else hi_default
    width = jitter_factor * hi
    nan = float("nan")
    if n < _MIN_EVENTS or not hi > lo:
        return ClockResult(nan, nan, nan, nan, n, nan, nan, nan, width, False,
                           f"not tested: fewer than {_MIN_EVENTS} events or an empty period range", lo, hi)
    grid = _make_grid(times, lo, hi, oversample)
    power = _fft_power(times, grid)
    observed = float(power.max())
    best_f = float(grid.freq[int(np.argmax(power))])

    if null not in ("rotate", "jitter"):
        raise ValueError("null must be 'rotate' or 'jitter'")
    make = rotation_surrogates if null == "rotate" else jitter_surrogates
    live = _LiveTime(events.exposure) if events.exposure else None
    live_times = live.to_live(times) if live is not None else times

    def draw(m: int) -> np.ndarray:
        fake = make(live_times, width, rng, m)
        if live is not None:
            fake = np.sort(live.to_real(fake), axis=1)
        return _fft_power(fake, grid).max(axis=1)

    p = sequential_p(observed, draw, batch=min(25, n_surrogates), max_draws=n_surrogates)

    span = float(times[-1] - times[0])
    f = _refine(times, best_f, 1.0 / (grid.length * grid.dt))
    base_period = 1.0 / f
    _, _, base_scatter = _phase_stats(times, base_period)
    period, scatter, multiple = base_period, base_scatter, 1
    for k in range(2, int(np.floor(hi / base_period)) + 1):
        fk = _refine(times, f / k, 1.0 / (k * max(span, 1e-12) * 2.0))
        _, _, s_k = _phase_stats(times, 1.0 / fk)
        if s_k <= 1.25 * base_scatter + 1e-12:
            period, scatter, multiple = 1.0 / fk, s_k, k
    exact, _ = _exact_power(times, np.array([1.0 / period]))
    _, phase, _ = _phase_stats(times, period)
    residual = (times - phase) / period
    residual_rms = float(np.sqrt(np.mean((residual - np.rint(residual)) ** 2)))
    bits = to_bits(events, period, phase)
    known = bits.bits if bits.known is None else bits.bits[bits.known]
    occupancy = float(known.mean()) if known.size else nan
    if multiple == 1:
        note = ("the best-power period was kept: no multiple of it up to max_period keeps the events as tightly "
                "on the ticks")
    else:
        note = (f"the best power fell at {base_period:.6g}; the reported period is {multiple} times that, the "
                f"longest multiple up to max_period with the events' scatter about the ticks no more than 25% "
                f"worse (a clock is also concentrated at its sub-multiples)")
    detected = bool(p.p <= alpha)
    return ClockResult(float(period), phase, float(exact[0]), p, n, span / period, occupancy, residual_rms, width,
                       detected, note, lo, hi)


# ---------------------------------------------------------------------------------------------------------------
# Slot occupancy
# ---------------------------------------------------------------------------------------------------------------

def to_bits(events: EventList, period: float, phase: float, n_slots: int | None = None) -> BitStream:
    """Slot occupancy on a clock: bit k is 1 when at least one event's nearest tick is tick k, else 0.

    Ticks are at phase + k * period. Which ticks count: with `events.exposure`, every tick lying at least 0.45
    periods inside an exposure interval (so its slot, tick +- half a period, is observed to within 5% of a
    period), from the first such tick to the last; ticks in exposure gaps are marked unknown (`known` False).
    Without exposure, from the first event's tick to the last event's tick. n_slots, if given, fixes the number
    of slots counted from the first one. meta records the period, phase, the time of slot 0, the number of
    events further than a quarter period from their tick (`n_far`, a measure of how well the clock fits) and
    the number of slots holding more than one event (`n_multiple`).
    """
    times = np.asarray(events.time, float)
    period, phase = float(period), float(phase)
    slot = np.rint((times - phase) / period).astype(np.int64)
    offset = (times - phase) / period - slot
    n_far = int(np.sum(np.abs(offset) > 0.25))
    if events.exposure:
        lo_t = min(a for a, _ in events.exposure)
        hi_t = max(b for _, b in events.exposure)
        first = int(np.ceil((lo_t + 0.45 * period - phase) / period))
        last = int(np.floor((hi_t - 0.45 * period - phase) / period))
    elif times.size:
        first, last = int(slot.min()), int(slot.max())
    else:
        first, last = 0, -1
    if n_slots is not None:
        last = first + int(n_slots) - 1
    count = max(last - first + 1, 0)
    inside = (slot >= first) & (slot <= last)
    hits = np.bincount(slot[inside] - first, minlength=count)[:count]
    bits = (hits > 0).astype(np.uint8)
    known = None
    if events.exposure and count:
        ticks = phase + (first + np.arange(count)) * period
        covered = np.zeros(count, dtype=bool)
        for a, b in events.exposure:
            covered |= (ticks >= a + 0.45 * period) & (ticks <= b - 0.45 * period)
        known = None if covered.all() else covered
    meta = {"period": period, "phase": phase, "slot0_time": phase + first * period, "n_events": int(times.size),
            "n_far": n_far, "n_multiple": int(np.sum(hits > 1)), "n_outside": int(np.sum(~inside))}
    return BitStream(bits, known, meta)


# ---------------------------------------------------------------------------------------------------------------
# Spacing test
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class SpacingResult:
    """The outcome of `spacing_test`.

    quantum: the quantum found (same unit as the event times).
    comb: the coherence C(q) = |sum_i exp(2 pi i d_i / q)| / N of the intervals at that quantum: about
        1 / sqrt(N) for intervals unrelated to q, 1 for intervals that are exact multiples of it.
    p_value: probability that intervals drawn from a smooth (comb-free) interval distribution reach the observed
        maximum coherence over the trial quanta.
    n_levels, levels: the multiples of the quantum that occur (at least twice, within a quarter quantum).
    detected: p_value <= alpha.
    """

    quantum: float
    comb: float
    p_value: MonteCarloP | float
    n_levels: int
    levels: np.ndarray
    detected: bool
    min_quantum: float = float("nan")
    max_quantum: float = float("nan")
    bandwidth: float = float("nan")
    note: str = ""


def _comb_grid(d: np.ndarray, q_min: float, q_max: float, oversample: float) -> _Grid:
    dt = q_min / 8.0
    reach = float(np.quantile(d, 0.9))
    length = sfft.next_fast_len(int(np.ceil(oversample * reach / dt)) + 1, real=True)
    length = min(length, _MAX_FFT)
    df = 1.0 / (length * dt)
    k_lo = max(1, int(np.ceil(1.0 / (q_max * df))))
    k_hi = min(length // 2, int(np.floor(1.0 / (q_min * df))))
    return _Grid(0.0, dt, length, k_lo, max(k_hi, k_lo))


def _comb_fft(d: np.ndarray, grid: _Grid) -> np.ndarray:
    """C(q) on the grid for intervals (N,) or a batch (M, N): the Rayleigh power machinery on the intervals."""
    n = np.shape(d)[-1]
    return np.sqrt(_fft_power(d, grid) * n / 2.0) / n


def _coherence(d: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Exact C(q) = |sum exp(2 pi i d / q)| / N at each trial quantum."""
    q = np.atleast_1d(np.asarray(q, float))
    total = np.zeros(q.size, dtype=complex)
    for chunk in range(0, d.size, 4096):
        total += np.exp(2j * np.pi * np.outer(1.0 / q, d[chunk:chunk + 4096])).sum(axis=1)
    return np.abs(total) / d.size


def smoothed_interval_surrogates(d: np.ndarray, bandwidth: float, rng: np.random.Generator, n: int) -> np.ndarray:
    """n sets of len(d) intervals drawn from a Gaussian kernel density estimate of log-interval.

    Each surrogate interval is d_j * exp(bandwidth * z) for a randomly chosen observed interval d_j and a
    standard normal z: the shape of the interval distribution on scales coarser than the bandwidth (in log
    units) is kept, including its behaviour near zero, and any comb finer than about bandwidth * d is erased.
    """
    pick = rng.integers(0, d.size, size=(n, d.size))
    return d[pick] * np.exp(bandwidth * rng.standard_normal((n, d.size)))


def spacing_test(events: EventList, rng: np.random.Generator, *, n_surrogates: int = 200, alpha: float = 1e-2,
                 bandwidth: float = 0.3, min_quantum: float | None = None, max_quantum: float | None = None,
                 oversample: float = 4) -> SpacingResult:
    """Are the intervals between consecutive events whole multiples of one quantum (pulse-position coding)?

    Statistic: C(q) = |sum_i exp(2 pi i d_i / q)| / N over the N intervals d_i, maximised over trial quanta q
    from min_quantum to max_quantum (computed by binned FFT, then refined with the exact sum).

    Null: shuffling the intervals leaves C unchanged, so the null is a smoothed bootstrap (see
    `smoothed_interval_surrogates`) with a kernel of `bandwidth` in log-interval. A Gaussian kernel of width
    s in linear units multiplies the coherence at quantum q by exp(-2 pi^2 s^2 / q^2), which is below 1% once
    s >= q / 2. In log units the kernel width at interval d is about bandwidth * d, so the comb at q is erased
    for every interval longer than q / (2 * bandwidth). The default largest trial quantum is therefore
    2 * bandwidth * (25th percentile of the intervals): at least three quarters of the intervals have their comb
    erased in the surrogates, and the shortest quarter keeps a little of it, which makes the null if anything
    conservative. A bandwidth of 0.3 in log units is a little narrower than the log-spread of any natural
    interval distribution (1.28 for exponential intervals), so the envelope of the distribution survives.

    A comb with only one or two levels (intervals all near q, or near q and 2 q) is a narrow or bimodal interval
    distribution, which natural sources can have; it is detected only through its sub-multiples q / j, which lie
    inside the searched range when the timing scatter is small enough, and the quantum is then reported as the
    longest multiple of the best trial quantum whose coherence is no more than 10% lower.

    If `clock_test` finds a clock, the intervals are multiples of its period and this test is nearly redundant:
    it will usually find the same period (or a multiple) as its quantum. The timing resolution of the instrument
    is itself a quantum: pass min_quantum above it when times are quantised. p >= 1 / (n_surrogates + 1).
    """
    times = np.sort(np.asarray(events.time, float))
    d = np.diff(times)
    d = d[d > 0]
    nan = float("nan")
    if d.size < _MIN_EVENTS:
        return SpacingResult(nan, nan, nan, 0, np.zeros(0, int), False, note="not tested: too few intervals")
    q_max = float(max_quantum) if max_quantum is not None else 2.0 * bandwidth * float(np.quantile(d, 0.25))
    q_min = float(min_quantum) if min_quantum is not None else q_max / 32.0
    if events.time_error is not None and min_quantum is None:
        q_min = max(q_min, 4.0 * float(np.median(events.time_error)))
    if not q_max > q_min:
        return SpacingResult(nan, nan, nan, 0, np.zeros(0, int), False, q_min, q_max, bandwidth,
                             "not tested: the trial quanta are below the timing resolution")
    grid = _comb_grid(d, q_min, q_max, oversample)
    comb = _comb_fft(d, grid)
    observed = float(comb.max())
    q_best = float(1.0 / grid.freq[int(np.argmax(comb))])

    def draw(m: int) -> np.ndarray:
        return _comb_fft(smoothed_interval_surrogates(d, bandwidth, rng, m), grid).max(axis=1)

    p = sequential_p(observed, draw, batch=min(25, n_surrogates), max_draws=n_surrogates)

    reach = float(np.quantile(d, 0.9))
    q = _refine_quantum(d, q_best, q_best ** 2 / (oversample * reach))
    c_base = float(_coherence(d, q)[0])
    quantum, c_q = q, c_base
    for k in range(2, min(64, int(np.floor(np.median(d) / q))) + 1):
        qk = _refine_quantum(d, k * q, k * q * q / (2.0 * reach))
        ck = float(_coherence(d, qk)[0])
        if ck >= 0.9 * c_base:
            quantum, c_q = qk, ck
    multiple = np.rint(d / quantum).astype(np.int64)
    close = np.abs(d / quantum - multiple) < 0.25
    values, counts = np.unique(multiple[close & (multiple > 0)], return_counts=True)
    levels = values[counts >= 2]
    detected = bool(p.p <= alpha)
    note = "" if detected else "no quantum stands out from a smooth interval distribution"
    return SpacingResult(float(quantum), c_q, p, int(levels.size), levels, detected, q_min, q_max, bandwidth, note)


def _refine_quantum(d: np.ndarray, q0: float, half_width: float, steps: int = 2) -> float:
    q = q0
    for _ in range(steps):
        trial = q + np.linspace(-half_width, half_width, 41)
        trial = trial[trial > 0]
        q = float(trial[int(np.argmax(_coherence(d, trial)))])
        half_width /= 10.0
    return q
