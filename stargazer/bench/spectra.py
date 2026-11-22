"""Synthetic dynamic spectra with known answers, for testing the narrowband drift search.

Each generator returns (Spectrogram, truth) or, for a cadence, (list of Spectrograms, truth). truth["signals"]
lists every tone put in, with its frequency at the start of the scan (f0_hz), drift, ideal matched S/N (as
defined by `detect.narrowband.inject_track`) and origin ("planted" for a stand-in transmitter, "interference"
for anything a search should not report as a candidate).

Noise is what an ideal spectrometer produces: every cell is an average of n_avg exponentially distributed raw
power samples (gamma-distributed with shape n_avg), times a smooth bandpass. Real data add things these do not
model: channel leakage, scalloping, a DC spike, quantisation, non-stationary interference. The default geometry
is that of Breakthrough Listen's high-spectral-resolution product (2.79 Hz channels, 18.25 s integrations).

These are self-contained and deliberately not registered in the bit-stream zoo.
"""
from __future__ import annotations

from typing import Any, Sequence

import numpy as np

from ..detect.narrowband import DriftConfig, Spectrogram, inject_tone, inject_track

T_SAMP = 18.253611
F_STEP = -2.793968
F_START = 1.42e9
N_AVG = 51
MJD0 = 60000.0


def _bandpass(n_freq: int, rng: np.random.Generator, period_channels: float = 16384.0) -> np.ndarray:
    """A gentle bandpass: a slope of 10% per 16,384 channels and an 8% ripple with that period.

    Its shape is fixed in channels (as a real bandpass is fixed in hertz), so a narrow band sees a small piece
    of it; the running median of the search follows it to about 0.1%.
    """
    k = np.arange(n_freq) / period_channels
    phase = rng.uniform(0, 2 * np.pi)
    return 1.0 + 0.1 * (k - k.mean()) + 0.08 * np.sin(2 * np.pi * k + phase)


def noise(rng: np.random.Generator, n_time: int = 16, n_freq: int = 4096, *, t_samp: float = T_SAMP,
          f_start: float = F_START, f_step: float = F_STEP, n_avg: int = N_AVG, bandpass: bool = True,
          level: float = 1000.0, tstart_mjd: float = MJD0, source_id: str = "noise", record_n_avg: bool = True
          ) -> tuple[Spectrogram, dict]:
    """Pure noise: gamma(n_avg) cells times a bandpass. With record_n_avg=False, meta does not say what n_avg is,
    so the search has to estimate it."""
    shape = _bandpass(n_freq, rng) if bandpass else np.ones(n_freq)
    power = level * shape[None, :] * rng.gamma(n_avg, 1.0 / n_avg, size=(n_time, n_freq))
    meta: dict[str, Any] = {"tstart": tstart_mjd}
    if record_n_avg:
        meta["n_avg"] = n_avg
    spec = Spectrogram(power, t_samp, f_start, f_step, source_id, meta)
    return spec, {"kind": "noise", "signals": [], "n_avg": n_avg}


def _channel_hz(spec: Spectrogram, channel: float) -> float:
    return spec.f_start + channel * spec.f_step


def drifting_tones(rng: np.random.Generator, tones: Sequence[dict] | None = None, n_tones: int = 3,
                   snr: float = 20.0, max_drift_hz_per_s: float = 4.0, **noise_kw: Any
                   ) -> tuple[Spectrogram, dict]:
    """Noise plus linearly drifting tones.

    `tones` is a list of dicts with f0_hz (or channel), drift_hz_per_s, snr and optionally width_hz. Without it,
    n_tones tones are drawn at random start channels, spaced apart, with drifts uniform in
    +-max_drift_hz_per_s, each kept inside the band for the whole scan.
    """
    spec, truth = noise(rng, **noise_kw)
    truth["kind"] = "drifting tones"
    if tones is None:
        sweep = max_drift_hz_per_s * spec.duration / abs(spec.f_step)
        margin = int(np.ceil(sweep)) + 8
        slots = np.linspace(margin, spec.n_freq - margin, n_tones + 2)[1:-1]
        tones = []
        for centre in slots:
            drift = float(rng.uniform(-max_drift_hz_per_s, max_drift_hz_per_s))
            channel = float(centre + rng.uniform(-0.5, 0.5)) - 0.5 * drift * spec.duration / spec.f_step
            tones.append({"channel": channel, "drift_hz_per_s": drift, "snr": snr})
    for tone in tones:
        f0 = float(tone["f0_hz"]) if "f0_hz" in tone else _channel_hz(spec, float(tone["channel"]))
        spec = inject_tone(spec, f0, float(tone["drift_hz_per_s"]), float(tone["snr"]),
                           width_hz=float(tone.get("width_hz", 0.0)))
        truth["signals"].append({"f0_hz": f0, "channel": (f0 - spec.f_start) / spec.f_step,
                                 "drift_hz_per_s": float(tone["drift_hz_per_s"]), "snr": float(tone["snr"]),
                                 "width_hz": float(tone.get("width_hz", 0.0)), "origin": "planted"})
    return spec, truth


def constant_interference(rng: np.random.Generator, snr: float = 50.0, channel: float | None = None,
                          **noise_kw: Any) -> tuple[Spectrogram, dict]:
    """Noise plus a zero-drift tone: the signature of local interference (or of a transmitter that corrects
    for its own acceleration, which is why it cannot simply be thrown away)."""
    spec, truth = noise(rng, **noise_kw)
    channel = float(rng.integers(100, spec.n_freq - 100)) if channel is None else float(channel)
    f0 = _channel_hz(spec, channel)
    spec = inject_tone(spec, f0, 0.0, snr)
    truth["kind"] = "constant interference"
    truth["signals"].append({"f0_hz": f0, "channel": channel, "drift_hz_per_s": 0.0, "snr": snr,
                             "origin": "interference"})
    return spec, truth


def broadband_pulse(rng: np.random.Generator, row: int | None = None, strength: float = 50.0,
                    **noise_kw: Any) -> tuple[Spectrogram, dict]:
    """Noise with one integration made brighter across the whole band.

    strength is the excess in every cell of that row, in units of the per-cell noise standard deviation. A sum
    over w channels in each of n_time rows would gain strength x sqrt(w / n_time) in S/N in every channel at
    once (12.5 at zero drift and about 65 for a 27-channel boxcar, with the defaults), unless the search
    divides out each row's level, which it does.
    """
    spec, truth = noise(rng, **noise_kw)
    row = int(rng.integers(0, spec.n_time)) if row is None else int(row)
    power = spec.power.copy()
    power[row] *= 1.0 + strength / np.sqrt(truth["n_avg"])
    truth["kind"] = "broadband pulse"
    truth["pulse_row"] = row
    truth["strength"] = strength
    return spec.with_power(power), truth


def wandering_interferer(rng: np.random.Generator, snr: float = 40.0, max_drift_hz_per_s: float = 1.0,
                         channel: float | None = None, **noise_kw: Any) -> tuple[Spectrogram, dict]:
    """Noise plus a tone whose drift changes sign halfway: df/dt runs linearly from -max to +max.

    No straight line fits it, so a linear de-Doppler search either misses it or finds it in pieces, with less
    than its ideal S/N. truth["signals"][0]["snr"] is the ideal matched S/N along the true (curved) track.
    """
    spec, truth = noise(rng, **noise_kw)
    channel = spec.n_freq / 2 if channel is None else float(channel)
    f0 = _channel_hz(spec, channel)
    span = spec.duration

    def freq(t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, float)
        return f0 + max_drift_hz_per_s * (t * t / span - t)

    spec = inject_track(spec, freq, snr)
    truth["kind"] = "wandering interferer"
    truth["signals"].append({"f0_hz": f0, "channel": channel, "drift_hz_per_s": float("nan"), "snr": snr,
                             "max_drift_hz_per_s": max_drift_hz_per_s, "origin": "interference",
                             "track": "f0 + max_drift x (t^2 / T - t)"})
    return spec, truth


def cadence(rng: np.random.Generator, n_time: int = 16, n_freq: int = 4096, *, slew_s: float = 20.0,
            candidate: dict | None = None, constant: dict | None = None, hidden: dict | None = None,
            config: DriftConfig = DriftConfig(), **noise_kw: Any) -> tuple[list[Spectrogram], dict]:
    """A six-scan ON, OFF, ON, OFF, ON, OFF cadence with three signals, each following one line through time.

    * candidate (planted): a drifting tone in the ON scans only, as a distant transmitter would be.
    * constant (interference): a zero-drift tone in all six scans.
    * hidden (interference): a drifting tone in all six scans, strong in the ONs but just below the usual
      threshold of 10 in the OFFs, as BLC1's look-alikes were. A filter that searches the OFF scans at the ON
      threshold lets it through.

    Each dict may set "channel" (at the start of the first scan), "drift_hz_per_s", "snr" and, for hidden,
    "snr_off". Scans are n_time integrations long and separated by slew_s seconds; meta["tstart"] is in MJD.
    truth["is_on"] gives the labels.
    """
    candidate = {"channel": 1300.5, "drift_hz_per_s": 0.31, "snr": 25.0, **(candidate or {})}
    constant = {"channel": 2600.0, "drift_hz_per_s": 0.0, "snr": 40.0, **(constant or {})}
    hidden = {"channel": 3100.25, "drift_hz_per_s": -0.17, "snr": 20.0, "snr_off": 7.0, **(hidden or {})}
    is_on = [True, False, True, False, True, False]
    scans: list[Spectrogram] = []
    length = n_time * noise_kw.get("t_samp", T_SAMP)
    tstart_mjd = noise_kw.pop("tstart_mjd", MJD0)
    starts = [k * (length + slew_s) for k in range(6)]
    f_of = None
    for k, on in enumerate(is_on):
        spec, _ = noise(rng, n_time, n_freq, tstart_mjd=tstart_mjd + starts[k] / 86400.0,
                        source_id=f"scan{k + 1}-{'ON' if on else 'OFF'}", **noise_kw)
        f_of = f_of or (lambda ch, s=spec: s.f_start + ch * s.f_step)
        for name, sig, present, snr in (("candidate", candidate, on, candidate["snr"]),
                                        ("constant", constant, True, constant["snr"]),
                                        ("hidden", hidden, True, hidden["snr"] if on else hidden["snr_off"])):
            if present and snr > 0:
                f0 = f_of(sig["channel"]) + sig["drift_hz_per_s"] * starts[k]
                spec = inject_tone(spec, f0, sig["drift_hz_per_s"], snr, config=config)
        scans.append(spec)
    truth = {"kind": "cadence", "is_on": is_on, "starts_s": starts, "scan_length_s": length, "signals": []}
    for name, sig, origin in (("candidate", candidate, "planted"), ("constant", constant, "interference"),
                              ("hidden", hidden, "interference")):
        truth["signals"].append({"name": name, "f0_hz": f_of(sig["channel"]), "channel": sig["channel"],
                                 "drift_hz_per_s": sig["drift_hz_per_s"], "snr": sig["snr"], "origin": origin,
                                 **({"snr_off": sig["snr_off"]} if "snr_off" in sig else {})})
    return scans, truth
