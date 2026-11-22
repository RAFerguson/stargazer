"""The control collection: signals that no structure test, judged against the right null, may call structured.

They set the false-alarm rate of everything else. Two of them (compressed and encrypted text) are real messages
whose bits are, by design, statistically indistinguishable from noise; an honest tool must report them as
noise-like and stop at content tier 0.
"""
from __future__ import annotations

import hashlib
import math
import zlib

import numpy as np

from stargazer.core import BitStream, ContentTier, EventList, Series, SymbolSequence

from .zoo import ZooItem, bytes_to_bits, english_text, expect, register

_N_BITS = {"small": 20_000, "full": 200_000}
_N_SYMBOLS = {"small": 10_000, "full": 100_000}
# Deflate output is not perfectly random: its bits carry an excess of ones of about 0.2 to 0.5% (Huffman codes
# have whole-bit lengths), which a monobit test starts to detect beyond roughly 15,000 bits, and an excess of runs
# (mean z of about +0.7 at 10,000 bits), so about 3% of texts fail a runs test at the 1% level. The compressed-text
# control is therefore kept at 3,000 characters (about 10,000 bits) at both sizes, and its text is redrawn
# (deterministically, from the item's own generator, at most _MAX_DRAWS times) until the stream passes the monobit
# and runs tests at the 1% level; the number of draws is recorded in truth["draws"]. A scorer that measures the
# false-alarm rate of a monobit or runs test should not use this item, which is conditioned on passing them.
_ZLIB_CHARS = 3_000
_MAX_DRAWS = 50
_CIPHER_CHARS = {"small": 3_000, "full": 15_000}


@register("control/fair_bits")
def fair_bits(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """Independent fair bits. truth: p_one, length, fair_null "iid"."""
    n = _N_BITS[size]
    bits = (rng.random(n) < 0.5).astype(np.uint8)
    return ZooItem("control/fair_bits", "control", BitStream(bits), {"p_one": 0.5, "length": n, "fair_null": "iid"},
                   expect(False, ContentTier.NONE, "noise"), "Independent fair coin flips.")


@register("control/biased_bits_p30")
def biased_bits(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """Independent bits with 30% ones. truth: p_one, length, fair_null "iid" (Bernoulli with fitted p)."""
    n = _N_BITS[size]
    p = 0.3
    bits = (rng.random(n) < p).astype(np.uint8)
    return ZooItem("control/biased_bits_p30", "control", BitStream(bits), {"p_one": p, "length": n, "fair_null": "iid"},
                   expect(False, ContentTier.NONE, "noise"), "Independent biased coin flips with 30% ones.")


def markov_bits(rng: np.random.Generator, n: int, transition: np.ndarray) -> np.ndarray:
    """A two-state Markov chain of length n started from its stationary distribution.

    transition[i, j] is the probability of moving from state i to state j.
    """
    p01, p10 = transition[0, 1], transition[1, 0]
    pi1 = p01 / (p01 + p10)
    u = rng.random(n)
    out = np.empty(n, dtype=np.uint8)
    state = int(u[0] < pi1)
    out[0] = state
    p_one_from = (transition[0, 1], transition[1, 1])
    for i in range(1, n):
        state = int(u[i] < p_one_from[state])
        out[i] = state
    return out


@register("control/markov1_bits")
def markov1_bits(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A first-order Markov bit stream. Structured only against an independent-bits null, so it is judged
    against a Markov-1 null and must not be called structured. truth: transition matrix, stationary p_one,
    fair_null "markov1"."""
    n = _N_BITS[size]
    transition = np.array([[0.8, 0.2], [0.3, 0.7]])
    bits = markov_bits(rng, n, transition)
    truth = {"transition_matrix": transition, "state_order": (0, 1), "p_one_stationary": 0.2 / 0.5, "length": n,
             "fair_null": "markov1"}
    return ZooItem("control/markov1_bits", "control", BitStream(bits), truth, expect(False, ContentTier.NONE, "noise"),
                   "A first-order Markov bit stream, noise relative to a Markov-1 null.")


def _iid_symbols(name: str, rng: np.random.Generator, size: str, alphabet: int) -> ZooItem:
    n = _N_SYMBOLS[size]
    symbols = rng.integers(0, alphabet, size=n)
    return ZooItem(name, "control", SymbolSequence(symbols, alphabet),
                   {"alphabet_size": alphabet, "length": n, "probabilities": np.full(alphabet, 1.0 / alphabet),
                    "fair_null": "iid"},
                   expect(False, ContentTier.NONE, "noise"),
                   f"Independent uniform symbols over an alphabet of {alphabet}.")


@register("control/iid_symbols_4")
def iid_symbols_4(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """Independent uniform symbols over 4 letters. truth: alphabet_size, probabilities, fair_null "iid"."""
    return _iid_symbols("control/iid_symbols_4", rng, size, 4)


@register("control/iid_symbols_16")
def iid_symbols_16(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """Independent uniform symbols over 16 letters. truth: alphabet_size, probabilities, fair_null "iid"."""
    return _iid_symbols("control/iid_symbols_16", rng, size, 16)


def monobit_p(bits: np.ndarray) -> float:
    """Two-sided p-value of the frequency (monobit) test: is the count of ones consistent with p = 1/2?"""
    b = np.asarray(bits, dtype=np.int64)
    return math.erfc(abs(2 * int(b.sum()) - b.size) / math.sqrt(2 * b.size))


def runs_p(bits: np.ndarray) -> float:
    """Two-sided p-value of the Wald-Wolfowitz runs test in the NIST SP 800-22 form (given the fraction of ones)."""
    b = np.asarray(bits, dtype=np.int64)
    n, pi = b.size, float(b.mean())
    if pi in (0.0, 1.0):
        return 0.0
    runs = 1 + int(np.count_nonzero(np.diff(b)))
    return math.erfc(abs(runs - 2 * n * pi * (1 - pi)) / (2 * math.sqrt(2 * n) * pi * (1 - pi)))


@register("control/zlib_text")
def zlib_text(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """English-like text compressed by zlib at level 9, as bits (most significant first), with the 2-byte zlib
    header and the 4-byte Adler-32 trailer stripped: a raw deflate stream. A real message that an honest tool must
    report as noise-like. truth: the text, level, the stripped header and trailer, format "raw_deflate". The text
    length is the same at both sizes, and the text is redrawn until the stream passes the monobit and runs
    tests at the 1% level (see the note above _ZLIB_CHARS); truth["draws"] counts the draws."""
    for draw in range(1, _MAX_DRAWS + 1):
        text = english_text(rng, _ZLIB_CHARS)
        full = zlib.compress(text.encode("ascii"), 9)
        header, body, trailer = full[:2], full[2:-4], full[-4:]
        bits = bytes_to_bits(body)
        if monobit_p(bits) > 0.01 and runs_p(bits) > 0.01:
            break
    else:  # pragma: no cover - fifty failures in a row would mean deflate changed
        raise RuntimeError("no compressed text passed the monobit and runs tests")
    truth = {"text": text, "level": 9, "format": "raw_deflate", "wbits": -15, "stripped_header": header.hex(),
             "stripped_trailer": trailer.hex(), "known_bias": "excess of ones of about 0.2 to 0.5%",
             "bit_order": "msb_first", "length": 8 * len(body), "draws": draw,
             "conditioned_on": "monobit and runs tests passing at the 1% level", "fair_null": "iid"}
    return ZooItem("control/zlib_text", "control", BitStream(bits), truth,
                   expect(False, ContentTier.NONE, "designed"),
                   "English-like text compressed with zlib level 9, headers stripped; must read as noise.")


def sha256_ctr_keystream(key: bytes, n_bytes: int) -> bytes:
    """Keystream: SHA-256(key || counter), counter an 8-byte big-endian integer from 0, blocks concatenated."""
    blocks = [hashlib.sha256(key + i.to_bytes(8, "big")).digest() for i in range((n_bytes + 31) // 32)]
    return b"".join(blocks)[:n_bytes]


@register("control/encrypted_text")
def encrypted_text(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """English-like text XORed with a SHA-256 counter-mode keystream, as bits (most significant first). truth:
    the text, the key (hex), the keystream construction."""
    text = english_text(rng, _CIPHER_CHARS[size])
    key = rng.bytes(32)
    plain = np.frombuffer(text.encode("ascii"), dtype=np.uint8)
    stream = np.frombuffer(sha256_ctr_keystream(key, plain.size), dtype=np.uint8)
    cipher = (plain ^ stream).tobytes()
    truth = {"text": text, "key": key.hex(), "keystream": "sha256(key || uint64_be(counter)), counter from 0",
             "bit_order": "msb_first", "length": 8 * len(cipher), "fair_null": "iid"}
    return ZooItem("control/encrypted_text", "control", BitStream(bytes_to_bits(cipher)), truth,
                   expect(False, ContentTier.NONE, "designed"),
                   "English-like text XORed with a SHA-256 counter-mode keystream; must read as noise.")


@register("control/white_noise_series")
def white_noise_series(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """White Gaussian noise about a constant magnitude at irregular survey sampling (the natural collection's
    seasonal cadence). truth: mean_mag, the per-point errors are the true noise sigmas, fair_null "white"."""
    from .natural import survey_times

    time = survey_times(rng, n_seasons=7 if size == "small" else 10)
    error = rng.uniform(0.01, 0.03, time.size)
    value = 17.0 + error * rng.standard_normal(time.size)
    truth = {"mean_mag": 17.0, "error_range": (0.01, 0.03), "fair_null": "white"}
    return ZooItem("control/white_noise_series", "control",
                   Series(time, value, error, band="g", source_id="control/white_noise_series", kind="mag",
                          time_standard="MJD"),
                   truth, expect(False, ContentTier.NONE, "noise"),
                   "White noise about a constant magnitude at irregular seasonal sampling.")


@register("control/poisson_events")
def poisson_events(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A homogeneous Poisson process of events. truth: rate (per second), duration, fair_null "poisson"."""
    rate, duration = 1.0, 2_000.0 if size == "small" else 20_000.0
    n = int(rng.poisson(rate * duration))
    time = np.sort(rng.uniform(0.0, duration, n))
    events = EventList(time, exposure=[(0.0, duration)], source_id="control/poisson_events", time_unit="s")
    return ZooItem("control/poisson_events", "control", events,
                   {"rate": rate, "duration": duration, "fair_null": "poisson"},
                   expect(False, ContentTier.NONE, "noise"),
                   "A homogeneous Poisson event list at one event per second.")
