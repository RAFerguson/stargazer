"""Injection and recovery: adding a known signal to real data to measure what a search can and cannot see.

A search that finds nothing has only said something if it could have found something. Every run therefore
injects synthetic signals into the very series it searched and reports the fraction recovered as a function
of the signal's strength and shape. That table is what turns "nothing found" into an upper limit.

All injections return a new Series and leave the input untouched. Amplitudes are in the unit of the series.
For magnitudes a positive `depth` in `on_off_keying` means fainter.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from ..core import Series


def sinusoid(series: Series, period_days: float, amplitude: float, phase_cycles: float = 0.0) -> Series:
    return series.with_values(series.value + amplitude * np.sin(2 * np.pi * (series.time / period_days + phase_cycles)))


def sawtooth(series: Series, period_days: float, amplitude: float, phase_cycles: float = 0.0) -> Series:
    """A slow rise and sudden fall each cycle: the shape simulations predict for accreting black-hole pairs."""
    frac = (series.time / period_days + phase_cycles) % 1.0
    return series.with_values(series.value + amplitude * (2.0 * frac - 1.0))


def periodic_flares(series: Series, period_days: float, amplitude: float, duty: float = 0.1,
                    phase_cycles: float = 0.0) -> Series:
    """A Gaussian flare once per cycle, lasting the fraction `duty` of the period (full width at half maximum)."""
    frac = (series.time / period_days + phase_cycles + 0.5) % 1.0 - 0.5
    width = duty / 2.355
    return series.with_values(series.value + amplitude * np.exp(-0.5 * (frac / width) ** 2))


def on_off_keying(series: Series, bits: np.ndarray, slot_days: float, depth: float,
                  start: float | None = None, repeat: bool = True) -> Series:
    """A coded modulation: during each slot the source is shifted by `depth` if the slot's bit is one.

    The bit pattern starts at `start` (default: the first observation) and repeats if it is shorter than
    the data. This is the simplest deliberate signal a slowly sampled light curve could carry.
    """
    bits = np.asarray(bits).astype(int).ravel()
    origin = series.time[0] if start is None else start
    slot = np.floor((series.time - origin) / slot_days).astype(int)
    inside = slot >= 0
    if repeat:
        value = bits[slot % bits.size]
    else:
        inside &= slot < bits.size
        value = bits[np.clip(slot, 0, bits.size - 1)]
    return series.with_values(series.value + depth * value * inside)


@dataclass
class RecoveryCell:
    """One cell of a completeness table."""

    label: str
    parameters: dict
    n_trials: int
    n_recovered: int

    @property
    def completeness(self) -> float:
        return self.n_recovered / self.n_trials if self.n_trials else float("nan")

    def interval(self) -> tuple[float, float]:
        """Wilson 95% interval on the completeness."""
        n, k, z = self.n_trials, self.n_recovered, 1.96
        if n == 0:
            return float("nan"), float("nan")
        centre = (k + z * z / 2) / (n + z * z)
        half = z * np.sqrt(k * (n - k) / n + z * z / 4) / (n + z * z)
        return float(centre - half), float(centre + half)

    def to_dict(self) -> dict:
        low, high = self.interval()
        return {"label": self.label, **self.parameters, "trials": self.n_trials, "recovered": self.n_recovered,
                "completeness": self.completeness, "low95": low, "high95": high}


def recovery_table(hosts: list[Series], inject: Callable[..., Series], grid: list[dict],
                   recovered: Callable[[Series, dict], bool], *, label: str = "") -> list[RecoveryCell]:
    """Inject every signal in `grid` into every host series and count how often `recovered` says yes.

    `inject(series, **parameters)` builds the injected series; `recovered(series, parameters)` runs the
    search and decides, by the same rule the real search uses, whether the injected signal was found.
    Hosts should be real series that the search itself passed over, so that the noise is the real noise.
    """
    cells = []
    for parameters in grid:
        hits = sum(bool(recovered(inject(host, **parameters), parameters)) for host in hosts)
        cells.append(RecoveryCell(label, dict(parameters), len(hosts), hits))
    return cells
