"""The measurement containers every module exchanges.

All are plain dataclasses over numpy arrays. They validate their shapes on construction, sort by time
where time exists, and convert to and from plain dictionaries so they can be cached and hashed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Provenance:
    """Where a piece of data came from: enough to fetch it again and to cite it."""

    source: str                      # registry id, for example "ztf_dr24_lc"
    query: dict[str, Any]
    retrieved_utc: str
    content_sha256: str
    citation: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {"source": self.source, "query": self.query, "retrieved_utc": self.retrieved_utc,
                "content_sha256": self.content_sha256, "citation": self.citation}

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Provenance":
        return cls(d["source"], d.get("query", {}), d.get("retrieved_utc", ""), d.get("content_sha256", ""),
                   d.get("citation", ""))


@dataclass
class Series:
    """A time-ordered measurement of one source in one band, with an error on every point.

    time is in days (any Julian-date flavour; `time_standard` says which). value and error are in the same
    unit, given by `kind`: "mag" (astronomical magnitudes, smaller is brighter), "flux" or "rate".
    """

    time: np.ndarray
    value: np.ndarray
    error: np.ndarray
    band: str = ""
    source_id: str = ""
    kind: str = "mag"
    time_standard: str = ""
    meta: dict[str, Any] = field(default_factory=dict)
    provenance: Provenance | None = None

    def __post_init__(self) -> None:
        self.time = np.asarray(self.time, dtype=float)
        self.value = np.asarray(self.value, dtype=float)
        self.error = np.asarray(self.error, dtype=float)
        if not (self.time.shape == self.value.shape == self.error.shape) or self.time.ndim != 1:
            raise ValueError("time, value and error must be one-dimensional arrays of equal length")
        if not (np.all(np.isfinite(self.time)) and np.all(np.isfinite(self.value)) and np.all(np.isfinite(self.error))):
            raise ValueError("time, value and error must be finite; drop bad points before building a Series")
        if np.any(self.error <= 0):
            raise ValueError("every error must be positive")
        order = np.argsort(self.time, kind="stable")
        self.time, self.value, self.error = self.time[order], self.value[order], self.error[order]

    def __len__(self) -> int:
        return int(self.time.size)

    @property
    def span(self) -> float:
        """Time between the first and last point, in days."""
        return float(self.time[-1] - self.time[0]) if len(self) else 0.0

    def select(self, mask: np.ndarray) -> "Series":
        """A new Series holding the points where mask is true."""
        mask = np.asarray(mask, dtype=bool)
        return Series(self.time[mask], self.value[mask], self.error[mask], self.band, self.source_id, self.kind,
                      self.time_standard, dict(self.meta), self.provenance)

    def with_values(self, value: np.ndarray) -> "Series":
        """A copy with the same sampling and errors and new values (used for simulations and injections)."""
        return Series(self.time.copy(), np.asarray(value, dtype=float), self.error.copy(), self.band, self.source_id,
                      self.kind, self.time_standard, dict(self.meta), self.provenance)

    def binned(self, width_days: float = 1.0) -> "Series":
        """Inverse-variance weighted means in bins of the given width (default: one per night)."""
        if not len(self):
            return self
        index = np.floor((self.time - np.floor(self.time[0])) / width_days).astype(np.int64)
        _, start = np.unique(index, return_index=True)
        weight = 1.0 / self.error ** 2
        w_sum = np.add.reduceat(weight, start)
        t = np.add.reduceat(weight * self.time, start) / w_sum
        v = np.add.reduceat(weight * self.value, start) / w_sum
        meta = dict(self.meta)
        meta["binned_days"] = width_days
        return Series(t, v, 1.0 / np.sqrt(w_sum), self.band, self.source_id, self.kind, self.time_standard, meta,
                      self.provenance)

    def seasons(self, gap_days: float = 60.0) -> np.ndarray:
        """An integer season label for every point; a new season starts after any gap longer than gap_days."""
        if not len(self):
            return np.zeros(0, dtype=int)
        return np.concatenate([[0], np.cumsum(np.diff(self.time) > gap_days)]).astype(int)

    def to_dict(self) -> dict[str, Any]:
        return {"time": self.time.tolist(), "value": self.value.tolist(), "error": self.error.tolist(),
                "band": self.band, "source_id": self.source_id, "kind": self.kind,
                "time_standard": self.time_standard, "meta": self.meta,
                "provenance": self.provenance.to_dict() if self.provenance else None}

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Series":
        prov = Provenance.from_dict(d["provenance"]) if d.get("provenance") else None
        return cls(d["time"], d["value"], d["error"], d.get("band", ""), d.get("source_id", ""), d.get("kind", "mag"),
                   d.get("time_standard", ""), d.get("meta", {}), prov)


@dataclass
class EventList:
    """Time-tagged events: photon arrival times, pulse arrival times, burst times.

    amplitude and time_error are optional per-event arrays. exposure, when known, is a list of
    (start, stop) intervals during which the instrument was looking; tests need it to tell silence from absence.
    """

    time: np.ndarray
    amplitude: np.ndarray | None = None
    time_error: np.ndarray | None = None
    exposure: list[tuple[float, float]] | None = None
    source_id: str = ""
    time_unit: str = "s"
    meta: dict[str, Any] = field(default_factory=dict)
    provenance: Provenance | None = None

    def __post_init__(self) -> None:
        self.time = np.asarray(self.time, dtype=float)
        if self.time.ndim != 1 or not np.all(np.isfinite(self.time)):
            raise ValueError("time must be a finite one-dimensional array")
        order = np.argsort(self.time, kind="stable")
        self.time = self.time[order]
        for name in ("amplitude", "time_error"):
            arr = getattr(self, name)
            if arr is not None:
                arr = np.asarray(arr, dtype=float)
                if arr.shape != self.time.shape:
                    raise ValueError(f"{name} must match time in length")
                setattr(self, name, arr[order])

    def __len__(self) -> int:
        return int(self.time.size)

    @property
    def intervals(self) -> np.ndarray:
        """Waiting times between consecutive events."""
        return np.diff(self.time)


@dataclass
class BitStream:
    """A sequence of bits of unknown format. `known` marks bits that were actually received (False = erased)."""

    bits: np.ndarray
    known: np.ndarray | None = None
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.bits = np.asarray(self.bits).astype(np.uint8).ravel()
        if self.bits.size and self.bits.max() > 1:
            raise ValueError("bits must be 0 or 1")
        if self.known is not None:
            self.known = np.asarray(self.known, dtype=bool).ravel()
            if self.known.shape != self.bits.shape:
                raise ValueError("known must match bits in length")

    def __len__(self) -> int:
        return int(self.bits.size)

    @classmethod
    def from_text(cls, text: str) -> "BitStream":
        """Build from a string of '0' and '1' characters; anything else is ignored, '?' marks an erased bit."""
        chars = [c for c in text if c in "01?"]
        bits = np.array([1 if c == "1" else 0 for c in chars], dtype=np.uint8)
        known = np.array([c != "?" for c in chars], dtype=bool)
        return cls(bits, None if known.all() else known)

    def to_text(self) -> str:
        out = self.bits.astype(str)
        if self.known is not None:
            out = np.where(self.known, out, "?")
        return "".join(out.tolist())


@dataclass
class SymbolSequence:
    """An ordered sequence of discrete symbols, optionally uncertain.

    symbols holds the most probable symbol at each position (integers in [0, alphabet_size)).
    probabilities, when present, has shape (length, alphabet_size): the probability of each symbol at each
    position. A measured signal turned into symbols is never certain, and statistics computed on the sequence
    should be averaged over `sample()` draws so that uncertainty reaches the final number.
    """

    symbols: np.ndarray
    alphabet_size: int
    probabilities: np.ndarray | None = None
    time: np.ndarray | None = None
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.symbols = np.asarray(self.symbols, dtype=np.int64).ravel()
        if self.symbols.size and (self.symbols.min() < 0 or self.symbols.max() >= self.alphabet_size):
            raise ValueError("symbols must lie in [0, alphabet_size)")
        if self.probabilities is not None:
            self.probabilities = np.asarray(self.probabilities, dtype=float)
            if self.probabilities.shape != (self.symbols.size, self.alphabet_size):
                raise ValueError("probabilities must have shape (length, alphabet_size)")
        if self.time is not None:
            self.time = np.asarray(self.time, dtype=float)
            if self.time.shape != self.symbols.shape:
                raise ValueError("time must match symbols in length")

    def __len__(self) -> int:
        return int(self.symbols.size)

    def sample(self, rng: np.random.Generator) -> np.ndarray:
        """One draw of a certain sequence from the per-position probabilities (the symbols themselves if certain)."""
        if self.probabilities is None:
            return self.symbols.copy()
        cumulative = np.cumsum(self.probabilities, axis=1)
        cumulative /= cumulative[:, -1:]
        return (rng.random((len(self), 1)) > cumulative).sum(axis=1).clip(0, self.alphabet_size - 1)
