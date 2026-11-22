"""Channel impairments for benchmark payloads, and `impaired` to apply one to a ZooItem.

Every channel returns a new payload and never modifies its input. Each also appends a record of what it did to
the new payload's `meta["channels"]` (for example the positions it flipped or deleted), which `impaired` copies
into the item's truth so a scorer can account for it.
"""
from __future__ import annotations

import copy
from typing import Any, Callable

import numpy as np

from stargazer.core import BitStream, EventList, Series

from .zoo import ZooItem


def _log(meta: dict[str, Any], record: dict[str, Any]) -> dict[str, Any]:
    out = dict(meta)
    out["channels"] = list(meta.get("channels", [])) + [record]
    return out


def _bits_out(stream: BitStream, bits: np.ndarray, known: np.ndarray | None, record: dict[str, Any]) -> BitStream:
    return BitStream(np.array(bits, dtype=np.uint8, copy=True), None if known is None else np.array(known, copy=True),
                     _log(stream.meta, record))


def flip_bits(stream: BitStream, rate: float, rng: np.random.Generator) -> BitStream:
    """Independent bit flips (a binary symmetric channel) with probability `rate` per bit."""
    flips = rng.random(len(stream)) < rate
    return _bits_out(stream, stream.bits ^ flips.astype(np.uint8), stream.known,
                     {"channel": "flip_bits", "rate": rate, "flipped": np.flatnonzero(flips)})


def erase_bits(stream: BitStream, fraction: float, rng: np.random.Generator) -> BitStream:
    """Erase a random `fraction` of the bits (exactly round(fraction * n)): marked unknown in `known`, set to 0."""
    n = len(stream)
    erased = np.zeros(n, dtype=bool)
    erased[rng.choice(n, size=int(round(fraction * n)), replace=False)] = True
    known = (np.ones(n, dtype=bool) if stream.known is None else stream.known.copy()) & ~erased
    bits = np.where(erased, 0, stream.bits)
    return _bits_out(stream, bits, known, {"channel": "erase_bits", "fraction": fraction,
                                           "erased": np.flatnonzero(erased)})


def prepend_random(stream: BitStream, n: int, rng: np.random.Generator) -> BitStream:
    """Put n random fair bits in front of the stream, so that its alignment is unknown."""
    head = (rng.random(n) < 0.5).astype(np.uint8)
    known = None if stream.known is None else np.concatenate([np.ones(n, dtype=bool), stream.known])
    return _bits_out(stream, np.concatenate([head, stream.bits]), known,
                     {"channel": "prepend_random", "n": n, "offset": n})


def truncate(stream: BitStream, n: int) -> BitStream:
    """Keep only the first n bits."""
    known = None if stream.known is None else stream.known[:n]
    return _bits_out(stream, stream.bits[:n], known, {"channel": "truncate", "n": n})


def drop_bits(stream: BitStream, fraction: float, rng: np.random.Generator) -> BitStream:
    """Delete a random `fraction` of the bits (exactly round(fraction * n)); everything after a deletion shifts."""
    n = len(stream)
    deleted = np.zeros(n, dtype=bool)
    deleted[rng.choice(n, size=int(round(fraction * n)), replace=False)] = True
    known = None if stream.known is None else stream.known[~deleted]
    return _bits_out(stream, stream.bits[~deleted], known, {"channel": "drop_bits", "fraction": fraction,
                                                             "deleted": np.flatnonzero(deleted)})


def add_noise(series: Series, sigma: float, rng: np.random.Generator) -> Series:
    """Add white Gaussian noise of standard deviation sigma; errors grow in quadrature."""
    value = series.value + sigma * rng.standard_normal(len(series))
    return Series(series.time.copy(), value, np.sqrt(series.error ** 2 + sigma ** 2), series.band, series.source_id,
                  series.kind, series.time_standard, _log(series.meta, {"channel": "add_noise", "sigma": sigma}),
                  series.provenance)


def drop_points(series: Series, fraction: float, rng: np.random.Generator) -> Series:
    """Remove a random `fraction` of the points (exactly round(fraction * n))."""
    n = len(series)
    drop = np.zeros(n, dtype=bool)
    drop[rng.choice(n, size=int(round(fraction * n)), replace=False)] = True
    out = series.select(~drop)
    out.meta = _log(series.meta, {"channel": "drop_points", "fraction": fraction, "dropped": np.flatnonzero(drop)})
    return out


def _events_out(events: EventList, time: np.ndarray, keep: np.ndarray, time_error: np.ndarray | None,
                record: dict[str, Any]) -> EventList:
    amp = None if events.amplitude is None else events.amplitude[keep].copy()
    exposure = None if events.exposure is None else list(events.exposure)
    return EventList(time, amp, time_error, exposure, events.source_id, events.time_unit, _log(events.meta, record),
                     events.provenance)


def jitter_times(events: EventList, sigma: float, rng: np.random.Generator) -> EventList:
    """Add Gaussian timing noise of standard deviation sigma to every event (events re-sorted afterwards)."""
    keep = np.ones(len(events), dtype=bool)
    time = events.time + sigma * rng.standard_normal(len(events))
    te = None if events.time_error is None else np.sqrt(events.time_error ** 2 + sigma ** 2)
    return _events_out(events, time, keep, te, {"channel": "jitter_times", "sigma": sigma})


def drop_events(events: EventList, fraction: float, rng: np.random.Generator) -> EventList:
    """Remove a random `fraction` of the events (exactly round(fraction * n))."""
    n = len(events)
    drop = np.zeros(n, dtype=bool)
    drop[rng.choice(n, size=int(round(fraction * n)), replace=False)] = True
    te = None if events.time_error is None else events.time_error[~drop].copy()
    return _events_out(events, events.time[~drop].copy(), ~drop, te,
                       {"channel": "drop_events", "fraction": fraction, "dropped": np.flatnonzero(drop)})


_TAGS = {"flip_bits": "flip", "erase_bits": "erase", "prepend_random": "prepend", "truncate": "trunc",
         "drop_bits": "dropbits", "add_noise": "noise", "drop_points": "droppts", "jitter_times": "jitter",
         "drop_events": "dropev"}


def impaired(item: ZooItem, channel: Callable[..., Any], **kwargs: Any) -> ZooItem:
    """A copy of item with `channel(item.payload, **kwargs)` applied.

    The name gains a suffix such as "+flip0.03"; truth gains a `channels` list holding each channel's record;
    `prepend_random` adds to `truth["offset"]` (where the original content now starts) and any change of a bit
    stream's length updates `truth["length"]`. The original item is not modified.
    """
    payload = channel(item.payload, **kwargs)
    label = _TAGS.get(channel.__name__, channel.__name__)
    params = [v for k, v in kwargs.items() if k != "rng"]
    suffix = label + ("" if not params else "_".join(f"{p:g}" if isinstance(p, float) else str(p) for p in params))
    truth = copy.deepcopy(item.truth)
    record = copy.deepcopy(payload.meta.get("channels", [{}])[-1]) if isinstance(payload.meta, dict) else {}
    truth["channels"] = list(truth.get("channels", [])) + [record]
    if channel.__name__ == "prepend_random":
        truth["offset"] = int(truth.get("offset", 0)) + int(record["offset"])
    if isinstance(payload, BitStream) and len(payload) != len(item.payload):
        truth["length"] = len(payload)
    elif isinstance(payload, (Series, EventList)) and len(payload) != len(item.payload):
        truth["n_points"] = len(payload)
    return ZooItem(f"{item.name}+{suffix}", item.collection, payload, truth, copy.deepcopy(item.expected),
                   f"{item.description.rstrip('.')}, passed through {channel.__name__}.")
