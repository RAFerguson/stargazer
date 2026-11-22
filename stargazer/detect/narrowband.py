"""Narrowband drift search: a tone whose frequency drifts linearly with time in a dynamic spectrum.

A transmitter on a planet that rotates and orbits, seen from a telescope on a planet that rotates and orbits,
reaches us as a tone whose frequency changes almost linearly over a few minutes. The classic search
(Breakthrough Listen's turboSETI and its successors) adds up the power along every straight line through the
time x frequency plane, a "de-Doppler" search, and keeps the lines that stand far above the noise. This module
does that, with four things done carefully.

1. **What S/N means.** Each cell of a spectrogram is an average of `n_avg` raw power samples, so in pure noise
   it is gamma-distributed with shape `n_avg`, and a sum of K cells is gamma-distributed with shape K x
   `n_avg`. It is not Gaussian when `n_avg` is small. The bandpass shape is removed with a running median, the
   level of every integration is measured in coarse frequency blocks, and the S/N of a sum S of K normalised
   cells is (S - K) / sqrt(K / n_avg): the excess over the expected sum in units of its standard deviation
   under noise alone. Its local probability comes from the exact gamma tail (`Hit.p_local`). If the data do
   not say what `n_avg` is (`meta["n_avg"]`), it is estimated from the spread of the cells. Two corrections
   keep this honest for wide sums: the noise model is estimated again with the tracks of first-pass hits
   masked, because strong tones bias it; and the errors of the estimated bandpass and row levels, which every
   cell of one sum shares, are measured and added to the noise of a K-cell sum as a factor (1 + a K).
2. **No loss at high drift.** A tone drifting by more than one channel per integration spreads its power over
   several channels in every row. Summing one channel per row (as a plain Taylor tree does) then captures only
   a fraction of about 1 / |drift in channels per row| of the power, which is why Margot et al. (2023) found
   the standard fine-resolution pipeline recovered only about 6% of injected signals over +-8.9 Hz/s. Here each
   trial drift sums, in every row, the run of channels the tone actually crosses (a boxcar whose width grows
   with |drift|), and the noise scale is adjusted to the number of cells added.
3. **Calibration by decoys.** A drifting tone is coherent from row to row; noise is not. Rotating every row by
   its own random number of channels keeps the noise of every row and destroys any line, so the strongest hit
   of the real data is compared with the strongest hits of many rotated copies (`search`).
4. **Interference is a separate problem.** Almost every strong narrowband signal on Earth is human-made, and
   it is as coherent as an extraterrestrial one. Decoys do not calibrate against it; the ON/OFF cadence
   filter does (`cadence_filter`, `cadence_search`).

What the search cannot see: tones whose frequency is not close to linear over one scan (the wandering
interferer of `bench.spectra` is found weakly or in pieces); tones wider than a few channels (the boxcar is
built for a pure tone); drifts beyond `max_drift_hz_per_s`; signals much shorter than a scan (a pulse
switched on for one integration is diluted by the others); and anything below the S/N threshold. Two signals
whose tracks cross within one scan are reported as one, the stronger.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from functools import lru_cache
from typing import Any, Sequence

import numpy as np
from scipy import special, stats

from ..trust.decoys import MonteCarloP, sequential_p

SECONDS_PER_DAY = 86400.0


@dataclass(frozen=True)
class DriftConfig:
    """Everything that defines a narrowband drift search. Hash it and commit the hash before running."""

    max_drift_hz_per_s: float = 4.0
    snr_threshold: float = 10.0
    n_decoys: int = 200                  # row-rotation decoys for the global probability (the most drawn)
    boxcar: bool = True                  # sum every channel the tone crosses in each row (False: one per row)
    drift_step: float = 1.0              # trial drift spacing, as the change in total displacement over the scan,
    #                                      in channels (1.0: neighbouring trials end one channel apart)
    start_offsets: int = 2               # trial start positions per channel: 2 tries whole and half channels
    bandpass_window: int = 1025          # channels in the running median that removes the bandpass shape
    block_channels: int = 1024           # coarse blocks in which the level of every integration is measured
    normalisation: str = "gamma"         # "gamma": exact gamma tail; "empirical": robust per-drift rescaling
    fit_inflation: bool = True           # widen the noise of wide sums by the measured cost of the noise model
    inflation_probes: int = 48           # drifts at which that cost is measured before the main pass
    mask_signals: bool = True            # re-estimate the noise model with the tracks of first-pass hits masked
    max_masked_fraction: float = 0.3     # give up masking (and say so) if hits cover more than this
    hit_guard_channels: float = 2.0      # tracks closer than this (beyond their widths) belong to one signal
    max_candidates: int = 20000          # channels above threshold considered by the hit grouping, strongest first
    max_hits: int = 1000
    keep_plane_cells: int = 2 ** 24      # keep the full (drift x frequency) plane when it is no larger than this
    decoy_batch: int = 10
    stop_exceedances: int = 10
    off_snr_threshold: float = 5.0       # threshold for hits in OFF scans: interference can hide just below 10
    cadence_freq_tolerance_channels: float = 4.0
    cadence_drift_tolerance_steps: float = 4.0


@dataclass
class Spectrogram:
    """A dynamic spectrum: power in n_time integrations by n_freq channels, with what is needed to place it.

    Channel k is centred on f_start + k x f_step Hz (f_step is often negative in radio data). Row t covers the
    time from t x t_samp to (t + 1) x t_samp seconds after the start of the scan. `meta` may carry "n_avg", the
    number of independent raw power samples averaged into each cell (for Stokes I, both polarisations count),
    and "tstart", the start of the scan as a Modified Julian Date, which the cadence search uses.
    """

    power: np.ndarray
    t_samp: float
    f_start: float
    f_step: float
    source_id: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.power = np.asarray(self.power, dtype=float)
        if self.power.ndim != 2 or min(self.power.shape) < 1:
            raise ValueError("power must be a two-dimensional (n_time, n_freq) array")
        if not np.all(np.isfinite(self.power)) or np.any(self.power < 0):
            raise ValueError("power must be finite and non-negative")
        if not self.t_samp > 0 or self.f_step == 0:
            raise ValueError("t_samp must be positive and f_step non-zero")

    @property
    def n_time(self) -> int:
        return int(self.power.shape[0])

    @property
    def n_freq(self) -> int:
        return int(self.power.shape[1])

    @property
    def duration(self) -> float:
        """Length of the scan in seconds."""
        return self.n_time * self.t_samp

    def frequencies(self) -> np.ndarray:
        """Centre frequency of every channel, in Hz."""
        return self.f_start + self.f_step * np.arange(self.n_freq)

    def with_power(self, power: np.ndarray) -> "Spectrogram":
        """A copy with the same geometry and new power (used for injections and decoys)."""
        return Spectrogram(power, self.t_samp, self.f_start, self.f_step, self.source_id, dict(self.meta))


@dataclass
class Hit:
    """One narrowband signal found in one scan.

    frequency_hz is the frequency at the start of the scan; drift_hz_per_s is its rate of change. snr is the
    de-drifted S/N defined in the module docstring. p_local is the probability that this one trial (one start
    channel, one drift, one start offset) would reach this S/N in noise alone; it carries no correction for the
    millions of trials searched, which is what `SearchResult.decoy` is for.
    """

    frequency_hz: float
    drift_hz_per_s: float
    snr: float
    channel: float              # start position in channels (fractional when a half-channel offset fitted best)
    drift_index: int            # index into DedopplerResult.drift_hz_per_s
    width_channels: float       # mean number of channels summed per integration at this drift
    p_local: float
    n_cells: int = 0            # number of cells added along the track

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class NoiseModel:
    """How the spectrogram was put on a common scale before summing.

    Normalised power is power / (bandpass[f] x level[t, block of f]); in noise alone every normalised cell is
    gamma-distributed with shape n_avg and mean 1.
    """

    n_avg: float
    n_avg_source: str                  # "meta" or "estimated from the quartiles of the cells"
    bandpass: np.ndarray               # (n_freq,) running median of each channel's median over time
    level: np.ndarray                  # (n_time, n_blocks) residual level of each integration in each block
    block_edges: np.ndarray            # (n_blocks + 1,) channel boundaries of the blocks
    excess_scale: np.ndarray           # (n_blocks,) spread of normalised cells / that expected from n_avg

    def summary(self) -> dict:
        return {"n_avg": float(self.n_avg), "n_avg_source": self.n_avg_source,
                "n_blocks": int(self.level.shape[1]),
                "max_excess_scale": float(np.nanmax(self.excess_scale)) if np.isfinite(self.excess_scale).any()
                else float("nan"),
                "row_level_range": [float(self.level.min()), float(self.level.max())]}


@dataclass
class DedopplerResult:
    """The de-Doppler plane and what is needed to read it.

    plane[j, c] is the S/N of the best of the start offsets for a track starting at channel c with drift
    drift_hz_per_s[j]; it is kept only when no larger than `DriftConfig.keep_plane_cells` (None otherwise).
    snr_max and best_drift always hold, for every channel, the best S/N over all drifts and its drift index.
    n_cells[j] is the number of cells summed by a full-length track at drift j: the noise scale of the raw sum
    is sqrt(n_cells / n_avg) in normalised units. snr_scale[j] is the robust spread of the S/N over channels
    at drift j (1 if the noise model is right). At high drift neighbouring channels share most of their cells,
    so with a few thousand channels it scatters by 5 to 10% about 1 from drift to drift even when all is well.
    """

    drift_hz_per_s: np.ndarray
    drift_channels: np.ndarray         # the same drifts in channels per integration (sign of channel index)
    n_cells: np.ndarray                # (n_drift,) cells of the whole-channel start offset
    cells_by_offset: np.ndarray        # (n_drift, start_offsets) cells of every trial
    snr_max: np.ndarray
    best_drift: np.ndarray
    best_offset: np.ndarray            # start offset (fraction of a channel) of the best trial per channel
    noise: NoiseModel
    plane: np.ndarray | None
    snr_scale: np.ndarray
    n_trials: int                      # start channels x drifts x start offsets
    n_time: int
    t_samp: float
    f_start: float
    f_step: float
    normalisation: str
    inflation: float = 0.0             # a in: variance of a K-cell sum = (K / n_avg)(1 + a K); see fit_inflation
    masked_fraction: float = 0.0       # share of cells masked as signal before the noise model was re-estimated

    @property
    def drift_step_hz_per_s(self) -> float:
        return float(abs(self.drift_hz_per_s[1] - self.drift_hz_per_s[0])) if self.drift_hz_per_s.size > 1 \
            else abs(self.f_step) / (self.n_time * self.t_samp)

    def p_local(self, snr: float | np.ndarray, drift_index: int | np.ndarray,
                offset_index: int | np.ndarray = 0) -> np.ndarray:
        """Probability of reaching `snr` or more in one trial (drift, start offset), in noise alone.

        Read from the exact gamma tail of a K-cell sum. Where an inflation a > 0 was fitted, the S/N has already
        been divided by sqrt(1 + a K) and the same tail is used, which treats the extra variance as if it were
        gamma-shaped: an approximation, good when a K is small (it is under 0.1 in the tests here).
        """
        snr = np.asarray(snr, float)
        if self.normalisation == "empirical":
            return stats.norm.sf(snr)
        shape = self.cells_by_offset[np.asarray(drift_index), np.asarray(offset_index)] * self.noise.n_avg
        return special.gammaincc(shape, shape + snr * np.sqrt(shape))


@dataclass
class SearchResult:
    """Hits in one scan, and the probability of the strongest one against noise."""

    source_id: str
    hits: list[Hit]
    statistic: float                   # S/N of the strongest trial anywhere in the plane
    decoy: MonteCarloP | None          # fraction of row-rotation decoys whose strongest trial matches it
    n_trials: int
    noise: dict
    flags: list[str] = field(default_factory=list)

    @property
    def strongest(self) -> Hit | None:
        return max(self.hits, key=lambda h: h.snr) if self.hits else None

    @property
    def p(self) -> float:
        return self.decoy.p if self.decoy else float("nan")

    def to_dict(self) -> dict:
        return {"source_id": self.source_id, "hits": [h.to_dict() for h in self.hits], "statistic": self.statistic,
                "decoy": self.decoy.to_dict() if self.decoy else None, "n_trials": self.n_trials,
                "noise": self.noise, "flags": list(self.flags)}


# ----------------------------------------------------------------------------------------------- noise model

@lru_cache(maxsize=256)
def _median_factor(n_avg: float) -> float:
    """Median of a gamma variable of shape n_avg and mean 1."""
    return float(stats.gamma.median(n_avg) / n_avg)


@lru_cache(maxsize=256)
def _iqr_factor(n_avg: float) -> float:
    """Interquartile range of a gamma variable of shape n_avg and mean 1."""
    return float((stats.gamma.ppf(0.75, n_avg) - stats.gamma.ppf(0.25, n_avg)) / n_avg)


@lru_cache(maxsize=1)
def _quartile_table() -> tuple[np.ndarray, np.ndarray]:
    shape = np.logspace(-1, 6, 600)
    ratio = (stats.gamma.ppf(0.75, shape) - stats.gamma.ppf(0.25, shape)) / stats.gamma.ppf(0.5, shape)
    return ratio[::-1], shape[::-1]           # increasing ratio for np.interp


def shape_from_cells(z: np.ndarray) -> float:
    """The gamma shape whose interquartile range over median matches that of the cells z.

    The ratio does not depend on the scale and is barely moved by a few bright cells, so it estimates the
    number of raw samples averaged per cell even when tones and interference are present. It assumes the
    cells are gamma-distributed; quantised or clipped data break that assumption.
    """
    q25, q50, q75 = np.quantile(z, [0.25, 0.5, 0.75])
    ratio, shape = _quartile_table()
    return float(np.interp((q75 - q25) / q50, ratio, shape))


def _block_edges(n_freq: int, block: int) -> np.ndarray:
    n_blocks = max(1, int(round(n_freq / block)))
    return np.linspace(0, n_freq, n_blocks + 1).round().astype(int)


def _running_median(values: np.ndarray, window: int) -> np.ndarray:
    """Running median over `window` channels, evaluated every window / 8 channels and interpolated.

    The ends are padded by odd reflection about a robust straight-line fit to the last half-window of channels
    (a linear continuation), so a sloping bandpass is not bent at the band edges. Reflecting about the single
    edge channel instead would put that one channel's noise into half of every window near the edge.
    """
    n = values.size
    window = int(min(window, n if n % 2 else n - 1))
    if window < 3:
        return values.copy()
    half = window // 2

    def pivot(end: np.ndarray) -> float:
        k = np.arange(end.size, dtype=float)
        keep = np.ones(end.size, dtype=bool)
        for _ in range(2):
            slope, intercept = np.polyfit(k[keep], end[keep], 1)
            resid = end - (intercept + slope * k)
            keep = np.abs(resid) <= 4.0 * 1.4826 * np.median(np.abs(resid)) + 1e-300
        return float(intercept)

    left, right = pivot(values[:half + 1]), pivot(values[::-1][:half + 1])
    padded = np.concatenate([2 * left - values[1:half + 1][::-1], values, 2 * right - values[-half - 1:-1][::-1]])
    step = max(1, window // 8)
    centres = np.unique(np.concatenate([np.arange(0, n, step), [n - 1]]))
    windows = np.lib.stride_tricks.sliding_window_view(padded, window)[centres]
    return np.interp(np.arange(n), centres, np.median(windows, axis=1))


def _fill_missing(values: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    """Replace non-finite entries by linear interpolation between finite neighbours (or by the fallback)."""
    bad = ~np.isfinite(values)
    if not bad.any():
        return values
    if bad.all():
        return fallback
    good = np.flatnonzero(~bad)
    out = values.copy()
    out[bad] = np.interp(np.flatnonzero(bad), good, values[good])
    return out


def normalise(spec: Spectrogram, config: DriftConfig = DriftConfig(), mask: np.ndarray | None = None
              ) -> tuple[np.ndarray, NoiseModel]:
    """Normalised power (mean 1 per cell in noise) and the noise model used.

    The bandpass is the running median (over `bandpass_window` channels) of the median over time of each
    channel. The median over time matters: a tone drifting fast holds any one channel for only an integration
    or two, so it does not raise that channel's median, whereas it would raise the time average over hundreds
    of channels and the bandpass would partly absorb it. Each integration's level in each block of
    `block_channels` channels is then measured by its median, converted to a mean with the gamma median-to-mean
    ratio, so a broadband pulse or a gain change in one integration is divided out rather than summed. n_avg is
    taken from `spec.meta["n_avg"]` when present and otherwise estimated from the quartiles of the cells
    (`shape_from_cells`). Cells where `mask` is True (the tracks of signals already found) are left out of every
    one of these estimates.

    The estimates are themselves noisy and the errors are shared by all the cells of a wide sum: the bandpass by
    about 1.6 / sqrt(n_time x n_avg x bandpass_window) per channel, the same in every row, and each row's level
    by about 1.25 / sqrt(n_avg x block_channels). They add a fraction of about K x n_avg x error^2 to the
    variance of a K-cell sum, which matters at high drift; `fit_inflation` measures and corrects it.
    """
    power = spec.power
    if mask is not None and not mask.any():
        mask = None
    held = power if mask is None else np.where(mask, np.nan, power)
    with np.errstate(all="ignore"), _quiet():
        per_channel = np.median(power, axis=0) if mask is None else np.nanmedian(held, axis=0)
    per_channel = _fill_missing(per_channel, np.median(power, axis=0))
    bandpass = _running_median(per_channel, config.bandpass_window)
    positive = bandpass[bandpass > 0]
    floor = float(np.median(positive)) * 1e-6 if positive.size else 1.0
    bandpass = np.where(bandpass > floor, bandpass, floor if positive.size else 1.0)
    y = power / bandpass
    y_held = y if mask is None else np.where(mask, np.nan, y)
    edges = _block_edges(spec.n_freq, config.block_channels)
    medians = np.empty((spec.n_time, edges.size - 1))
    with np.errstate(all="ignore"), _quiet():
        for b in range(edges.size - 1):
            part = y_held[:, edges[b]:edges[b + 1]]
            medians[:, b] = np.median(part, axis=1) if mask is None else np.nanmedian(part, axis=1)
    medians = np.where(np.isfinite(medians) & (medians > 0), medians, 1.0)
    block_of = np.repeat(np.arange(edges.size - 1), np.diff(edges))
    z = y_held / medians[:, block_of]
    z = z[np.isfinite(z)] if mask is not None else z
    if spec.meta.get("n_avg"):
        n_avg, source = float(spec.meta["n_avg"]), "meta"
    else:
        n_avg, source = shape_from_cells(z), "estimated from the quartiles of the cells"
    level = medians / _median_factor(n_avg)
    x = y / level[:, block_of]
    x_held = x if mask is None else np.where(mask, np.nan, x)
    expected_iqr = _iqr_factor(n_avg)
    excess = np.empty(edges.size - 1)
    for b in range(edges.size - 1):
        part = x_held[:, edges[b]:edges[b + 1]]
        part = part[np.isfinite(part)]
        excess[b] = np.subtract(*np.quantile(part, [0.75, 0.25])) / expected_iqr if part.size > 10 else np.nan
    return x, NoiseModel(n_avg, source, bandpass, level, edges, excess)


def signal_mask(hits: Sequence["Hit"], result: "DedopplerResult", n_freq: int, guard: int = 2) -> np.ndarray:
    """Cells (n_time, n_freq) on the tracks of the given hits, widened by `guard` channels on each side."""
    mask = np.zeros((result.n_time, n_freq), dtype=bool)
    for h in hits:
        c = int(np.floor(h.channel + 1e-9))
        lo, hi = track_ranges(float(result.drift_channels[h.drift_index]), h.channel - c, result.n_time, True)
        for t in range(result.n_time):
            if lo[t] <= hi[t]:
                mask[t, max(0, c + lo[t] - guard):max(0, min(n_freq, c + hi[t] + guard + 1))] = True
    return mask


class _quiet:
    """Silence numpy's warnings about all-NaN slices (handled explicitly)."""

    def __enter__(self) -> None:
        import warnings

        self._ctx = warnings.catch_warnings()
        self._ctx.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def __exit__(self, *exc: object) -> None:
        self._ctx.__exit__(*exc)


# ---------------------------------------------------------------------------------------------- trial tracks

def drift_grid(spec: Spectrogram, config: DriftConfig = DriftConfig()) -> tuple[np.ndarray, np.ndarray]:
    """Trial drifts, as (Hz per second, channels per integration), in increasing Hz per second.

    The spacing is `drift_step` / n_time channels per integration: neighbouring trials end the scan
    `drift_step` channels apart. Zero drift is always a trial.
    """
    step = config.drift_step / spec.n_time
    r_max = config.max_drift_hz_per_s * spec.t_samp / abs(spec.f_step)
    j_max = int(np.floor(r_max / step + 1e-9))
    channels = np.arange(-j_max, j_max + 1) * step
    hz = channels * spec.f_step / spec.t_samp
    order = np.argsort(hz, kind="stable")
    return hz[order], channels[order]


def _overlaps(start: float, rate: float, n_time: int) -> tuple[int, np.ndarray]:
    """Fraction of each row's tone power in each channel, for a tone at `start` + `rate` x t (t in rows).

    Channel k covers positions [k - 0.5, k + 0.5). Returns (first channel, array (n_time, n_channels)).
    """
    a = start + rate * np.arange(n_time)
    b = a + rate
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    k0 = int(np.floor(lo.min() + 0.5))
    k1 = int(np.floor(hi.max() + 0.5))
    k = np.arange(k0, k1 + 1)[None, :]
    if abs(rate) < 1e-12:
        out = (np.floor(lo[:, None] + 0.5) == k).astype(float)
    else:
        cover = np.minimum(hi[:, None], k + 0.5) - np.maximum(lo[:, None], k - 0.5)
        out = np.clip(cover, 0.0, None) / (hi - lo)[:, None]
    return k0, out


def track_ranges(rate: float, offset: float, n_time: int, boxcar: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """The channels summed in each row for a trial track, relative to its start channel.

    Returns (lo, hi), each (n_time,): row t sums channels c + lo[t] ... c + hi[t] (empty when lo > hi).

    With the boxcar, the set of cells is the one that maximises signal / sqrt(number of cells) for a pure tone
    starting `offset` channels after the centre of channel c: cells are taken in decreasing order of the share
    of tone power they hold, for as long as adding one raises the matched S/N. In each row this is a run of
    consecutive channels whose length grows with |rate|. Without it (the naive sum), each row contributes the
    one channel that holds the tone at the middle of the integration.
    """
    if not boxcar:
        centre = np.floor(offset + rate * (np.arange(n_time) + 0.5) + 0.5).astype(int)
        return centre, centre.copy()
    k0, share = _overlaps(offset, rate, n_time)
    flat = np.sort(share.ravel())[::-1]
    flat = flat[flat > 0]
    gain = np.cumsum(flat) / np.sqrt(np.arange(1, flat.size + 1))
    cut = flat[int(np.argmax(gain))]
    chosen = share >= cut - 1e-12
    lo = np.full(n_time, 1, dtype=int)
    hi = np.zeros(n_time, dtype=int)
    for t in range(n_time):
        idx = np.flatnonzero(chosen[t])
        if idx.size:
            lo[t], hi[t] = k0 + idx[0], k0 + idx[-1]
    return lo, hi


@lru_cache(maxsize=32)
def _trial_table(rates: tuple[float, ...], n_offsets: int, n_time: int, boxcar: bool
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """lo, hi (n_drift, n_offsets, n_time), cells (n_drift, n_offsets) and offsets for every trial."""
    offsets = np.arange(n_offsets) / n_offsets
    lo = np.empty((len(rates), n_offsets, n_time), dtype=int)
    hi = np.empty_like(lo)
    for j, rate in enumerate(rates):
        for o, off in enumerate(offsets):
            lo[j, o], hi[j, o] = track_ranges(rate, off, n_time, boxcar)
    cells = np.clip(hi - lo + 1, 0, None).sum(axis=2)
    return lo, hi, cells, offsets


# -------------------------------------------------------------------------------------------- the de-Doppler

def _edge_cells(lo: np.ndarray, hi: np.ndarray, channels: np.ndarray, n_freq: int) -> np.ndarray:
    """Number of in-band cells summed by tracks starting at `channels`."""
    a = np.maximum(channels[:, None] + lo[None, :], 0)
    b = np.minimum(channels[:, None] + hi[None, :], n_freq - 1)
    return np.clip(b - a + 1, 0, None).sum(axis=1)


@lru_cache(maxsize=4096)
def _lower_spread(shape: float) -> float:
    """(median - 15.87th percentile) of a standardised gamma variable of this shape (Gaussian: 1)."""
    if not np.isfinite(shape) or shape > 1e7:
        return 1.0
    q16, q50 = stats.gamma.ppf([stats.norm.cdf(-1.0), 0.5], shape)
    return float((q50 - q16) / np.sqrt(shape))


def _robust_scale(values: np.ndarray, shape: float = np.inf) -> float:
    """The spread of S/N values, measured from below the median only.

    Signals and their butterflies only ever raise the S/N, so the lower half of the distribution is clean even
    when strong tones touch many channels. The distance from the median down to the 15.87th percentile is
    compared with the same distance for a standardised gamma variable of the given shape (K x n_avg).
    """
    values = values[np.isfinite(values)]
    if values.size < 10:
        return float("nan")
    q16, q50 = np.quantile(values, [stats.norm.cdf(-1.0), 0.5])
    return float((q50 - q16) / _lower_spread(float(shape)))


class _Summer:
    """Sums along trial tracks of normalised power x, by cumulative sums along frequency."""

    def __init__(self, x: np.ndarray, lo: np.ndarray, hi: np.ndarray, cells: np.ndarray, n_avg: float) -> None:
        self.n_time, self.n_freq = x.shape
        self.lo, self.hi, self.cells, self.n_avg = lo, hi, cells, n_avg
        nonempty = hi >= lo
        self.nonempty = nonempty
        self.pad = int(max(1, np.max(np.abs(np.where(nonempty, lo, 0))),
                           np.max(np.abs(np.where(nonempty, hi, 0))))) + 1
        self.cumulative = np.zeros((self.n_time, self.n_freq + 2 * self.pad + 1))
        np.cumsum(np.pad(x, ((0, 0), (self.pad, self.pad))), axis=1, out=self.cumulative[:, 1:])
        self.channels = np.arange(self.n_freq)

    def total(self, j: int, o: int) -> np.ndarray:
        """The plain sum along trial (j, o) for every start channel."""
        n_freq, pad, cumulative = self.n_freq, self.pad, self.cumulative
        total = np.zeros(n_freq)
        for t in range(self.n_time):
            a, b = self.lo[j, o, t], self.hi[j, o, t]
            if a <= b:
                total += cumulative[t, pad + b + 1:pad + b + 1 + n_freq]
                total -= cumulative[t, pad + a:pad + a + n_freq]
        return total

    def snr(self, j: int, o: int) -> np.ndarray:
        """(sum - K) / sqrt(K / n_avg) for every start channel; tracks mostly out of band get -inf."""
        n_freq = self.n_freq
        total = self.total(j, o)
        k = float(self.cells[j, o])
        snr = (total - k) / np.sqrt(k / self.n_avg)
        row_lo, row_hi = self.lo[j, o][self.nonempty[j, o]], self.hi[j, o][self.nonempty[j, o]]
        left, right = max(0, -int(row_lo.min())), max(0, int(row_hi.max()))
        edge = np.unique(np.concatenate([self.channels[:min(left, n_freq)],
                                         self.channels[max(n_freq - right, 0):]]))
        if edge.size:
            k_edge = _edge_cells(row_lo, row_hi, edge, n_freq).astype(float)
            ok = k_edge >= 0.5 * k
            snr[edge] = np.where(ok, (total[edge] - k_edge) / np.sqrt(np.maximum(k_edge, 1) / self.n_avg), -np.inf)
        return snr


def _duplicates(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """True for trials (drift, offset) that repeat offset 0's cells shifted by whole channels (zero drift)."""
    out = np.zeros(lo.shape[:2], dtype=bool)
    for o in range(1, lo.shape[1]):
        shift = lo[:, o, :1] - lo[:, 0, :1]
        out[:, o] = np.all(lo[:, o] - lo[:, 0] == shift, axis=1) & np.all(hi[:, o] - hi[:, 0] == shift, axis=1)
    return out


def _sample(n_freq: int, cap: int = 16384) -> slice:
    return slice(0, n_freq, max(1, n_freq // cap))


def fit_inflation(scales: np.ndarray, cells: np.ndarray) -> float:
    """The coefficient a in variance = (K / n_avg) x (1 + a x K), fitted to measured S/N spreads.

    Errors in the estimated bandpass and row levels are shared by all the cells of one sum, so they add
    variance in proportion to K squared, while noise adds it in proportion to K. Given the robust spread of the
    S/N over channels at several drifts (scales) and their cell counts, this is the least-squares a >= 0.
    """
    ok = np.isfinite(scales)
    if not np.any(ok):
        return 0.0
    k = cells[ok].astype(float)
    return float(max(0.0, np.sum(k * (scales[ok] ** 2 - 1.0)) / np.sum(k * k)))


def _dedoppler_normalised(x: np.ndarray, noise: NoiseModel, hz: np.ndarray, rates: np.ndarray,
                          config: DriftConfig, keep_plane: bool, meta: dict, mask: np.ndarray | None = None
                          ) -> DedopplerResult:
    """The de-Doppler search of already-normalised power x (n_time, n_freq).

    Where `mask` marks the cells of signals already found, trials that touch them are left out of the spread
    measurements (the inflation fit and snr_scale): a strong tone raises the S/N of every track that crosses it,
    which can be a large share of all channels at a given drift.
    """
    n_time, n_freq = x.shape
    lo, hi, cells, offsets = _trial_table(tuple(float(r) for r in rates), max(1, config.start_offsets), n_time,
                                          config.boxcar)
    n_drift, n_off = cells.shape
    summer = _Summer(x, lo, hi, cells, noise.n_avg)
    skip = _duplicates(lo, hi)
    sample = _sample(n_freq)
    touched = None if mask is None or not mask.any() else _Summer(mask.astype(float), lo, hi, cells, 1.0)

    def clean(values: np.ndarray, j: int, o: int) -> np.ndarray:
        if touched is None:
            return values[sample]
        return values[sample][touched.total(j, o)[sample] < 0.5]

    empirical = config.normalisation == "empirical"
    inflation = 0.0
    if not empirical and config.fit_inflation and n_drift > 1:
        probe = np.unique(np.linspace(0, n_drift - 1, min(n_drift, config.inflation_probes)).round().astype(int))
        scales = np.array([_robust_scale(clean(summer.snr(j, 0), j, 0), cells[j, 0] * noise.n_avg)
                           for j in probe])
        inflation = fit_inflation(scales, cells[probe, 0])
    snr_max = np.full(n_freq, -np.inf)
    best_drift = np.zeros(n_freq, dtype=int)
    best_offset = np.zeros(n_freq)
    plane = np.empty((n_drift, n_freq), dtype=np.float32) if keep_plane else None
    snr_scale = np.empty(n_drift)
    for j in range(n_drift):
        row_best = np.full(n_freq, -np.inf)
        row_off = np.zeros(n_freq)
        for o in range(n_off):
            if skip[j, o]:
                continue
            snr = summer.snr(j, o)
            if empirical:
                finite = snr[np.isfinite(snr)]
                scale = _robust_scale(clean(snr, j, o))
                snr = (snr - np.median(finite)) / (scale if scale > 0 else 1.0)
            elif inflation > 0:
                snr = snr / np.sqrt(1.0 + inflation * cells[j, o])
            if o == 0:
                snr_scale[j] = _robust_scale(clean(snr, j, 0), np.inf if empirical else cells[j, 0] * noise.n_avg)
            better = snr > row_best
            row_best = np.where(better, snr, row_best)
            row_off = np.where(better, offsets[o], row_off)
        if plane is not None:
            plane[j] = row_best
        better = row_best > snr_max
        snr_max = np.where(better, row_best, snr_max)
        best_drift = np.where(better, j, best_drift)
        best_offset = np.where(better, row_off, best_offset)
    return DedopplerResult(hz, rates, cells[:, 0], cells, snr_max, best_drift, best_offset, noise, plane,
                           snr_scale, int(n_freq * np.sum(~skip)), n_time, float(meta["t_samp"]),
                           float(meta["f_start"]), float(meta["f_step"]), config.normalisation, inflation)


def _two_pass(spec: Spectrogram, config: DriftConfig, keep_plane: bool) -> tuple[np.ndarray, DedopplerResult]:
    """Normalise and search; if hits are found, mask their tracks, normalise again without them, search again."""
    hz, rates = drift_grid(spec, config)
    meta = {"t_samp": spec.t_samp, "f_start": spec.f_start, "f_step": spec.f_step}
    x, noise = normalise(spec, config)
    result = _dedoppler_normalised(x, noise, hz, rates, config, keep_plane, meta)
    if config.mask_signals:
        hits, _ = _hits_from(result, config.snr_threshold, config)
        if hits:
            mask = signal_mask(hits, result, spec.n_freq)
            if mask.mean() <= config.max_masked_fraction:
                x, noise = normalise(spec, config, mask)
                result = _dedoppler_normalised(x, noise, hz, rates, config, keep_plane, meta, mask)
                result.masked_fraction = float(mask.mean())
            else:
                result.masked_fraction = float("nan")
    return x, result


def dedoppler(spec: Spectrogram, config: DriftConfig = DriftConfig()) -> DedopplerResult:
    """The de-Doppler plane of one scan: the S/N of the sum along every trial track.

    A trial is a start channel, a drift from `drift_grid` and a start offset (whole or half channel). Its sum
    adds, in every row, the normalised power of the channels given by `track_ranges`; tracks that leave the
    band are scored on the cells that stay in it, and dropped when fewer than half do. The S/N is
    (sum - K) / sqrt(K / n_avg) for K cells: exact standardisation under the gamma noise model, divided by
    sqrt(1 + a K) where `fit_inflation` finds a > 0 (with `normalisation="empirical"`, each drift's S/N is
    instead re-centred and re-scaled by its own median and lower-side spread over channels, and probabilities
    are read from a Gaussian, which understates the tail when n_avg is small).

    Strong signals bias the noise model they are measured against (a tone at 3 Hz/s puts power into a few
    per cent of every row's cells). So when the first pass finds hits, their tracks are masked, the noise model
    is estimated again without them, and the search is repeated (`mask_signals`).

    Cost is about n_drift x start_offsets x n_time x 2 array additions of length n_freq, twice when there are
    hits: brute force, with cumulative sums so that a boxcar of any width costs the same as one channel.
    """
    keep = drift_grid(spec, config)[0].size * spec.n_freq <= config.keep_plane_cells
    return _two_pass(spec, config, keep)[1]


# ---------------------------------------------------------------------------------------------------- hits

def _hits_from(result: DedopplerResult, threshold: float, config: DriftConfig) -> tuple[list[Hit], list[str]]:
    flags: list[str] = []
    candidates = np.flatnonzero(result.snr_max >= threshold)
    if candidates.size == 0:
        return [], flags
    order = candidates[np.argsort(result.snr_max[candidates], kind="stable")[::-1]]
    if order.size > config.max_candidates:
        flags.append(f"{order.size} channels above threshold; only the strongest {config.max_candidates} grouped")
        order = order[:config.max_candidates]
    n_time = result.n_time
    starts, rates, widths, hits = [], [], [], []
    for c in order:
        j = int(result.best_drift[c])
        rate = float(result.drift_channels[j])
        start = float(c + result.best_offset[c])
        o = int(round(float(result.best_offset[c]) * result.cells_by_offset.shape[1]))
        cells = int(result.cells_by_offset[j, o])
        width = cells / n_time
        if starts:
            s0, r0, w0 = np.array(starts), np.array(rates), np.array(widths)
            d_start = start - s0
            d_end = d_start + (rate - r0) * n_time
            crosses = d_start * d_end <= 0
            near = np.minimum(np.abs(d_start), np.abs(d_end)) <= 0.5 * (w0 + width) + config.hit_guard_channels
            if np.any(crosses | near):
                continue
        starts.append(start)
        rates.append(rate)
        widths.append(width)
        snr = float(result.snr_max[c])
        hits.append(Hit(result.f_start + start * result.f_step, float(result.drift_hz_per_s[j]), snr, start, j,
                        width, float(result.p_local(snr, j, o)), cells))
        if len(hits) >= config.max_hits:
            flags.append(f"stopped at {config.max_hits} hits")
            break
    return hits, flags


def find_hits(spec: Spectrogram | DedopplerResult, config: DriftConfig = DriftConfig(),
              threshold: float | None = None) -> list[Hit]:
    """Signals above the S/N threshold, one per signal, strongest first.

    A strong tone raises the S/N of every trial track that shares some of its cells: neighbouring channels at
    the same drift and tracks at other drifts that cross it (the "butterfly"). Channels above threshold are
    taken strongest first, each with its best drift, and a channel is dropped if its track crosses, or comes
    within the two tracks' half-widths plus `hit_guard_channels` of, the track of a stronger hit already kept.
    The price is that two real signals whose tracks cross within the scan are reported as one.

    Accepts a spectrogram (searched with `dedoppler`) or a DedopplerResult already computed, so that the same
    plane can be read at two thresholds (as the cadence search does for OFF scans).
    """
    result = spec if isinstance(spec, DedopplerResult) else dedoppler(spec, config)
    return _hits_from(result, config.snr_threshold if threshold is None else threshold, config)[0]


def _rotated(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n_time, n_freq = x.shape
    shifts = rng.integers(0, n_freq, size=n_time)
    rows = np.arange(n_freq)[None, :] - shifts[:, None]
    return np.take_along_axis(x, rows % n_freq, axis=1)


def search(spec: Spectrogram, rng: np.random.Generator, config: DriftConfig = DriftConfig()) -> SearchResult:
    """Hits in one scan and the probability, against noise, of the strongest trial.

    The statistic is the largest S/N anywhere in the de-Doppler plane. Decoys rotate every row of the
    normalised spectrogram by its own random number of channels (a circular shift) and run the whole search
    again, normalisation included; the probability is the fraction of decoys whose largest S/N matches or beats
    the real one (`sequential_p`, at most `n_decoys`). It includes the look-elsewhere cost of every channel,
    drift and start offset searched, and it does not depend on the gamma noise model being right.

    The null is "no structure coherent from row to row". Rotation destroys a drifting tone, but it equally
    destroys zero-drift interference that is constant in frequency, a comb of instrumental lines, or any
    persistent bandpass feature. A tiny probability here therefore says "not noise", never "not
    interference". Telling a transmitter from interference is the job of the cadence filter. The errors of the
    estimated bandpass are coherent across rows and the decoys scramble them, so the real data's wide sums are
    noisier than the decoys'; the variance inflation fitted to each (`fit_inflation`) puts both back on a unit
    scale, and `DedopplerResult.snr_scale` shows how well.
    """
    keep = drift_grid(spec, config)[0].size * spec.n_freq <= config.keep_plane_cells
    x, result = _two_pass(spec, config, keep)
    noise = result.noise
    hits, flags = _hits_from(result, config.snr_threshold, config)
    finite = result.snr_max[np.isfinite(result.snr_max)]
    observed = float(finite.max()) if finite.size else float("nan")
    decoy_meta = {"n_avg": spec.meta["n_avg"]} if spec.meta.get("n_avg") else {}

    def draw(n: int) -> np.ndarray:
        out = np.empty(n)
        for i in range(n):
            fake = Spectrogram(_rotated(x, rng), spec.t_samp, spec.f_start, spec.f_step, spec.source_id,
                               decoy_meta)
            out[i] = np.max(_two_pass(fake, config, False)[1].snr_max)
        return out

    decoy = None
    if np.isfinite(observed) and config.n_decoys > 0:
        decoy = sequential_p(observed, draw, batch=config.decoy_batch, max_draws=config.n_decoys,
                             stop_exceedances=config.stop_exceedances)
    scale = result.snr_scale[np.isfinite(result.snr_scale)]
    if scale.size and (abs(np.median(scale) - 1.0) > 0.1 or scale.max() > 1.3):
        flags.append(f"S/N spread over channels departs from 1 (median {np.median(scale):.2f}, range "
                     f"{scale.min():.2f} to {scale.max():.2f}): the noise model is not exact here")
    if np.isfinite(noise.excess_scale).any() and np.nanmax(noise.excess_scale) > 1.15:
        flags.append("some blocks are noisier than the gamma model with this n_avg predicts")
    if noise.n_avg_source != "meta":
        flags.append(f"n_avg not given; estimated as {noise.n_avg:.1f} from the cells' quartiles")
    if not np.isfinite(result.masked_fraction):
        flags.append("hits cover too much of the band to mask; the noise model includes them")
    summary = noise.summary()
    summary["inflation"] = result.inflation
    summary["masked_fraction"] = result.masked_fraction
    summary["snr_scale_range"] = [float(scale.min()), float(scale.max())] if scale.size else []
    return SearchResult(spec.source_id, hits, observed, decoy, result.n_trials, summary, flags)


# -------------------------------------------------------------------------------------------------- injection

def track_cells(spec: Spectrogram, frequency_of_time: Any, width_hz: float = 0.0) -> np.ndarray:
    """Share of a tone's total power that falls in each cell, for a tone following frequency_of_time(t).

    frequency_of_time takes seconds since the start of the scan (an array) and returns Hz. The tone is
    followed at many instants within every integration; a tone of non-zero width spreads its power uniformly
    over width_hz. Power that falls outside the band is lost. Returns an (n_time, n_freq) array summing to at
    most 1. No spectral leakage between channels is modelled (a real filterbank's channel response is not a
    perfect boxcar, which smears a tone a little further).
    """
    n_time, n_freq = spec.power.shape
    probe = np.linspace(0.0, spec.duration, 4 * n_time + 1)
    pos = (np.asarray(frequency_of_time(probe), float) - spec.f_start) / spec.f_step
    width_ch = abs(width_hz / spec.f_step)
    sweep = np.max(np.abs(np.diff(pos))) * 4 if pos.size > 1 else 0.0
    n_sub = int(32 * max(1, np.ceil(sweep + width_ch)))
    t = (np.arange(n_time)[:, None] + (np.arange(n_sub)[None, :] + 0.5) / n_sub) * spec.t_samp
    p = (np.asarray(frequency_of_time(t.ravel()), float) - spec.f_start) / spec.f_step
    rows = np.repeat(np.arange(n_time), n_sub)
    out = np.zeros((n_time, n_freq))
    weight = 1.0 / (n_time * n_sub)
    if width_ch <= 0:
        k = np.floor(p + 0.5).astype(int)
        ok = (k >= 0) & (k < n_freq)
        np.add.at(out, (rows[ok], k[ok]), weight)
        return out
    a, b = p - width_ch / 2, p + width_ch / 2
    first = np.floor(a + 0.5).astype(int)
    for step in range(int(np.ceil(width_ch)) + 2):
        k = first + step
        cover = np.clip(np.minimum(b, k + 0.5) - np.maximum(a, k - 0.5), 0.0, None) / width_ch
        ok = (k >= 0) & (k < n_freq) & (cover > 0)
        np.add.at(out, (rows[ok], k[ok]), weight * cover[ok])
    return out


def inject_track(spec: Spectrogram, frequency_of_time: Any, snr: float, rng: np.random.Generator | None = None,
                 width_hz: float = 0.0, config: DriftConfig = DriftConfig()) -> Spectrogram:
    """A copy of `spec` with a tone following frequency_of_time(t) added at the given ideal S/N.

    snr is the expected S/N of the ideal matched sum: the sum over every cell the tone touches, each weighted
    by its share of the tone's power, divided by that sum's noise standard deviation. No search with
    unweighted sums can exceed it, so a recovered S/N divided by it is the efficiency of the search. The noise
    level at each cell comes from the same noise model the search uses (`normalise`).

    Without rng the expected tone power is added. With rng the extra scatter that a tone adds to the power it
    shares a cell with (variance 2 x tone / n_avg in normalised units) is drawn too.
    """
    _, noise = normalise(spec, config)
    share = track_cells(spec, frequency_of_time, width_hz)
    norm = float(np.sqrt(np.sum(share ** 2)))
    if norm == 0:
        return spec.with_power(spec.power.copy())
    amplitude = snr / (np.sqrt(noise.n_avg) * norm)
    tone = amplitude * share                                   # normalised units: noise has mean 1
    if rng is not None:
        tone = tone + np.sqrt(2.0 * tone / noise.n_avg) * rng.standard_normal(tone.shape)
    block_of = np.repeat(np.arange(noise.level.shape[1]), np.diff(noise.block_edges))
    scale = noise.bandpass[None, :] * noise.level[:, block_of]
    return spec.with_power(np.clip(spec.power + tone * scale, 0.0, None))


def inject_tone(spec: Spectrogram, f0_hz: float, drift_hz_per_s: float, snr: float,
                rng: np.random.Generator | None = None, width_hz: float = 0.0,
                config: DriftConfig = DriftConfig()) -> Spectrogram:
    """A copy of `spec` with a linearly drifting tone added: f(t) = f0_hz + drift_hz_per_s x t.

    f0_hz is the frequency at the start of the scan; snr is the ideal matched S/N defined in `inject_track`.
    """
    return inject_track(spec, lambda t: f0_hz + drift_hz_per_s * np.asarray(t, float), snr, rng, width_hz, config)


# ---------------------------------------------------------------------------------------------- the cadence

@dataclass
class CadenceCandidate:
    """A group of hits followed across the scans of a cadence, accepted or rejected, and why.

    The group starts from a hit in the first ON scan and follows the line f(T) = frequency_hz + drift x (T -
    T0) through the later scans, T0 being the start of the first ON scan. on_hits holds the hit matched in each
    ON scan (None where none was), off_hits the OFF-scan hits that lie on the line and so veto it.
    """

    accepted: bool
    reason: str
    frequency_hz: float
    drift_hz_per_s: float
    on_hits: list[Hit | None]
    off_hits: list[tuple[int, Hit]] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)

    @property
    def snr(self) -> float:
        """The weakest ON-scan S/N (nan if an ON scan has no match)."""
        values = [h.snr for h in self.on_hits if h is not None]
        return float(min(values)) if values and len(values) == len(self.on_hits) else float("nan")

    def to_dict(self) -> dict:
        return {"accepted": self.accepted, "reason": self.reason, "frequency_hz": self.frequency_hz,
                "drift_hz_per_s": self.drift_hz_per_s,
                "on_hits": [h.to_dict() if h else None for h in self.on_hits],
                "off_hits": [(k, h.to_dict()) for k, h in self.off_hits], "flags": list(self.flags)}


def _line_fit(times: list[float], freqs: list[float], drifts: list[float]) -> tuple[float, float]:
    """Frequency at T = 0 and drift of the line through matched hits (start frequencies at scan starts).

    With one point the hit's own drift is used; with more, the slope through the points, which is much more
    precise because the scans are minutes apart.
    """
    if len(times) == 1:
        return freqs[0] - drifts[0] * times[0], drifts[0]
    slope, intercept = np.polyfit(np.asarray(times), np.asarray(freqs), 1)
    return float(intercept), float(slope)


def cadence_filter(on: Sequence[Sequence[Hit]], off: Sequence[Sequence[Hit]], on_starts_s: Sequence[float],
                   off_starts_s: Sequence[float], scan_length_s: float, *, freq_tolerance_hz: float,
                   drift_tolerance_hz_per_s: float, off_snr_threshold: float = 0.0) -> list[CadenceCandidate]:
    """Follow hits across an ON/OFF cadence and keep those seen in every ON scan and in no OFF scan.

    on and off hold the hits of each ON and OFF scan; on_starts_s and off_starts_s give each scan's start in
    seconds on a common clock. For every hit of the first ON scan, strongest first, the line through it is
    extrapolated to each later ON scan; a hit there matches if its drift is within drift_tolerance_hz_per_s of
    the line's and its start frequency is within freq_tolerance_hz plus drift_tolerance_hz_per_s x (time
    elapsed) of the line's. After each match the line is refitted through all matched start frequencies. A
    group is rejected if any ON scan has no match, or if any OFF hit (at or above off_snr_threshold, whatever
    its drift) has a track that crosses or comes within the same frequency window of the line during its OFF
    scan. The OFF hits should come from a lower threshold than the ON hits: BLC1 and its look-alikes showed that
    interference can sit just below threshold in the OFF scans.

    The filter rejects interference that is present when the telescope looks away. It cannot reject
    interference whose on/off pattern happens to follow the cadence (BLC1), and it rejects a real transmitter
    that is not on throughout the ON scans. Its false-alarm rate is not computed here; `cadence_search` reports
    the count with ON and OFF labels swapped, which measures chance coincidence directly. Candidates at zero
    drift are flagged, since a signal with no drift at all is more often local than not.
    """
    if not on:
        return []
    t0 = float(on_starts_s[0])
    on_t = [float(t) - t0 for t in on_starts_s]
    off_t = [float(t) - t0 for t in off_starts_s]
    claimed: set[tuple[int, int]] = set()
    out: list[CadenceCandidate] = []
    for seed in sorted(on[0], key=lambda h: -h.snr):
        times, freqs, drifts = [0.0], [seed.frequency_hz], [seed.drift_hz_per_s]
        matched: list[Hit | None] = [seed]
        reason = ""
        for j in range(1, len(on)):
            f_ref, drift = _line_fit(times, freqs, drifts)
            predicted = f_ref + drift * on_t[j]
            window = freq_tolerance_hz + drift_tolerance_hz_per_s * abs(on_t[j] - times[-1])
            best, best_gap = None, np.inf
            for i, h in enumerate(on[j]):
                gap = abs(h.frequency_hz - predicted)
                if gap <= window and abs(h.drift_hz_per_s - drift) <= drift_tolerance_hz_per_s and gap < best_gap:
                    if (j, i) not in claimed:
                        best, best_gap = (i, h), gap
            if best is None:
                matched.append(None)
                reason = reason or (f"absent from ON scan {j + 1}: no unclaimed hit within {window:.1f} Hz of "
                                    f"{predicted:.1f} Hz at drift {drift:.4f} +- "
                                    f"{drift_tolerance_hz_per_s:.4f} Hz/s")
                continue
            claimed.add((j, best[0]))
            matched.append(best[1])
            times.append(on_t[j])
            freqs.append(best[1].frequency_hz)
            drifts.append(best[1].drift_hz_per_s)
        f_ref, drift = _line_fit(times, freqs, drifts)
        vetoes: list[tuple[int, Hit]] = []
        if not reason:
            for k, hits in enumerate(off):
                predicted = f_ref + drift * off_t[k]
                elapsed = min(abs(off_t[k] - t) for t in times) + scan_length_s
                window = freq_tolerance_hz + drift_tolerance_hz_per_s * elapsed
                for h in hits:
                    if h.snr < off_snr_threshold:
                        continue
                    d_start = h.frequency_hz - predicted
                    d_end = d_start + (h.drift_hz_per_s - drift) * scan_length_s
                    if d_start * d_end <= 0 or min(abs(d_start), abs(d_end)) <= window:
                        vetoes.append((k, h))
            if vetoes:
                k, h = max(vetoes, key=lambda kh: kh[1].snr)
                reason = (f"on the line in OFF scan {k + 1}: hit at {h.frequency_hz:.1f} Hz, drift "
                          f"{h.drift_hz_per_s:.4f} Hz/s, S/N {h.snr:.1f}"
                          + (f" (and {len(vetoes) - 1} more OFF hits)" if len(vetoes) > 1 else ""))
        accepted = not reason
        verdict = reason or f"in all {len(on)} ON scans and on the line in none of the {len(off)} OFF scans"
        candidate = CadenceCandidate(accepted, verdict, f_ref, drift, matched, vetoes)
        if accepted and abs(drift) < drift_tolerance_hz_per_s:
            candidate.flags.append("zero drift within tolerance: a transmitter that compensates for its "
                                   "acceleration, or, far more often, local interference")
        out.append(candidate)
    return out


@dataclass
class CadenceResult:
    """The outcome of a cadence search: every group, the hits it came from and the label-swap null count."""

    candidates: list[CadenceCandidate]
    on_hits: list[list[Hit]]
    off_hits: list[list[Hit]]
    swapped_accepted: int                # candidates accepted when ON and OFF labels are swapped
    tolerances: dict

    @property
    def accepted(self) -> list[CadenceCandidate]:
        return [c for c in self.candidates if c.accepted]


def scan_start_seconds(spec: Spectrogram) -> float:
    """Start of the scan in seconds, from meta["tstart"] (MJD) or meta["tstart_s"]."""
    if "tstart_s" in spec.meta:
        return float(spec.meta["tstart_s"])
    if "tstart" in spec.meta:
        return float(spec.meta["tstart"]) * SECONDS_PER_DAY
    raise ValueError(f"scan {spec.source_id!r} has no start time (meta['tstart'] in MJD or meta['tstart_s'])")


def cadence_search(scans: Sequence[Spectrogram], is_on: Sequence[bool], config: DriftConfig = DriftConfig()
                   ) -> CadenceResult:
    """Search every scan of a cadence and run `cadence_filter`.

    Each scan is de-Dopplered once. ON scans give hits at `snr_threshold`, OFF scans at the lower
    `off_snr_threshold`. Tolerances are `cadence_freq_tolerance_channels` channels and
    `cadence_drift_tolerance_steps` drift steps of the first scan. The same filter is then run with the labels
    swapped (OFF scans as ON, at the full threshold, and ON scans as OFF, at the low one); with three of each,
    what passes there passes by chance or by a cadence-shaped coincidence, and is reported as a count.
    """
    if len(scans) != len(is_on):
        raise ValueError("is_on must give one label per scan")
    results = [dedoppler(s, config) for s in scans]
    high = [_hits_from(r, config.snr_threshold, config)[0] for r in results]
    low = [_hits_from(r, config.off_snr_threshold, config)[0] for r in results]
    starts = [scan_start_seconds(s) for s in scans]
    first = scans[0]
    freq_tol = config.cadence_freq_tolerance_channels * abs(first.f_step)
    drift_tol = config.cadence_drift_tolerance_steps * results[0].drift_step_hz_per_s
    length = first.duration
    on_idx = [i for i, flag in enumerate(is_on) if flag]
    off_idx = [i for i, flag in enumerate(is_on) if not flag]
    kwargs = {"freq_tolerance_hz": freq_tol, "drift_tolerance_hz_per_s": drift_tol,
              "off_snr_threshold": config.off_snr_threshold}
    candidates = cadence_filter([high[i] for i in on_idx], [low[i] for i in off_idx], [starts[i] for i in on_idx],
                                [starts[i] for i in off_idx], length, **kwargs)
    swapped = cadence_filter([high[i] for i in off_idx], [low[i] for i in on_idx], [starts[i] for i in off_idx],
                             [starts[i] for i in on_idx], length, **kwargs)
    return CadenceResult(candidates, [high[i] for i in on_idx], [low[i] for i in off_idx],
                         sum(c.accepted for c in swapped),
                         {"freq_tolerance_hz": freq_tol, "drift_tolerance_hz_per_s": drift_tol,
                          "scan_length_s": length})
