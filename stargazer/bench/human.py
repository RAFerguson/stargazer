"""The human-made collection: engineered signals whose framing and coding a tool should find.

Content tiers. Framed and coded items stop at tier 3 (LAYOUT): their frame length, sync word, field boundaries,
line code or error-correcting code are recoverable with compression gain and held-out prediction. The plain 8N1
text is given tier 4 (ANCHORED). Judgement: the ASCII table is shared knowledge between the human sender and any
Stargazer decoder, the candidate set of standard character codes is small (ASCII, EBCDIC, Baudot, a few others),
so a match of decoded bytes to English words survives a trials correction by an enormous margin. It is anchoring
to a human convention, not to universal knowledge; that is exactly what separates a human-made signal from a
designed one. (The coded items carry the same text, and a decoder that undoes the code could anchor them too;
they are kept at 3 because the specification scores them on layout and coding recovery.)
"""
from __future__ import annotations

import numpy as np

from stargazer.core import BitStream, ContentTier, EventList

from .zoo import ZooItem, english_text, expect, int_to_bits, register, text_bits

_TEXT_CHARS = {"small": 600, "full": 3_000}


def frame_8n1(text: str) -> np.ndarray:
    """Asynchronous serial frames: per character one start bit 0, eight data bits least significant first, one
    stop bit 1."""
    codes = np.frombuffer(text.encode("ascii"), dtype=np.uint8)
    data = np.unpackbits(codes[:, None], axis=1, bitorder="little")
    frames = np.concatenate([np.zeros((codes.size, 1), np.uint8), data, np.ones((codes.size, 1), np.uint8)], axis=1)
    return frames.ravel()


_8N1_TRUTH = {"frame_len": 10, "start_bit": 0, "stop_bit": 1, "data_bits": 8, "parity": "none",
              "bit_order": "lsb_first", "fields": {"start": (0, 1), "data": (1, 9), "stop": (9, 10)}}


@register("human/ascii_8n1")
def ascii_8n1(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """ASCII text as asynchronous serial 8N1 frames, back to back. truth: frame_len 10, start/stop bits, bit
    order, field boundaries, the text."""
    text = english_text(rng, _TEXT_CHARS[size])
    truth = dict(_8N1_TRUTH, text=text, length=10 * len(text))
    return ZooItem("human/ascii_8n1", "human", BitStream(frame_8n1(text)), truth,
                   expect(True, ContentTier.ANCHORED, "human"), "ASCII text sent as 8N1 asynchronous serial frames.")


def manchester(bits: np.ndarray) -> np.ndarray:
    """Manchester code, G. E. Thomas convention: 1 -> 10, 0 -> 01."""
    bits = np.asarray(bits, dtype=np.uint8)
    return np.stack([bits, 1 - bits], axis=1).ravel()


@register("human/ascii_manchester")
def ascii_manchester(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """ASCII text (8 bits per character, most significant first) Manchester-coded. truth: the convention
    (1 -> 10, 0 -> 01), chip pairs start at bit 0, the data bits and the text."""
    text = english_text(rng, _TEXT_CHARS[size])
    data = text_bits(text)
    truth = {"line_code": "manchester", "mapping": {1: (1, 0), 0: (0, 1)}, "chips_per_bit": 2, "phase": 0,
             "char_bits": 8, "bit_order": "msb_first", "data": data, "text": text, "length": 2 * data.size}
    return ZooItem("human/ascii_manchester", "human", BitStream(manchester(data)), truth,
                   expect(True, ContentTier.LAYOUT, "human"), "ASCII text, Manchester-coded (1 -> 10, 0 -> 01).")


def crc16_ccitt(bits: np.ndarray, init: int = 0xFFFF, poly: int = 0x1021) -> int:
    """CRC-16-CCITT (poly 0x1021, init 0xFFFF, no reflection, no final XOR) over a bit sequence, MSB first."""
    crc = init
    for b in np.asarray(bits, dtype=int):
        top = ((crc >> 15) & 1) ^ int(b)
        crc = ((crc << 1) & 0xFFFF) ^ (poly if top else 0)
    return crc


SYNC_WORD = 0x1ACFFC1D
TELEMETRY_FIELDS = {"sync": (0, 32), "counter": (32, 48), "payload": (48, 480), "crc": (480, 496)}


def telemetry_frames(payloads: list[np.ndarray]) -> np.ndarray:
    """Frames of 496 bits: sync word 0x1ACFFC1D, 16-bit counter from 0, 432 payload bits, CRC-16-CCITT over the
    counter and payload (bits 32 to 479); every field most significant bit first."""
    sync = int_to_bits(SYNC_WORD, 32)
    frames = []
    for i, payload in enumerate(payloads):
        body = np.concatenate([int_to_bits(i % 65536, 16), payload])
        frames.append(np.concatenate([sync, body, int_to_bits(crc16_ccitt(body), 16)]))
    return np.concatenate(frames)


def _telemetry_truth(n_frames: int) -> dict:
    return {"frame_len": 496, "sync_word": SYNC_WORD, "sync_len": 32, "fields": dict(TELEMETRY_FIELDS),
            "counter_start": 0, "crc": {"name": "CRC-16-CCITT", "poly": 0x1021, "init": 0xFFFF, "reflect": False,
                                        "xorout": 0, "covers": (32, 480)},
            "bit_order": "msb_first", "n_frames": n_frames, "length": 496 * n_frames}


def _n_frames(size: str) -> int:
    return 40 if size == "small" else 400


@register("human/telemetry_text")
def telemetry_text(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """Repeated telemetry frames whose 432-bit payloads carry 54 ASCII characters of text each. truth: frame
    length, sync word, field boundaries, CRC definition, the text."""
    n = _n_frames(size)
    text = english_text(rng, 54 * n)
    bits = text_bits(text).reshape(n, 432)
    truth = dict(_telemetry_truth(n), text=text, payload_kind="ascii_text")
    return ZooItem("human/telemetry_text", "human", BitStream(telemetry_frames(list(bits))), truth,
                   expect(True, ContentTier.ANCHORED, "human"),
                   "Telemetry frames (sync 0x1ACFFC1D, counter, 432-bit text payload, CRC-16) back to back.")


@register("human/telemetry_random")
def telemetry_random(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """The telemetry envelope with random payload bits: framing without readable content. truth: as for
    human/telemetry_text, with the payload bits instead of a text."""
    n = _n_frames(size)
    payload = (rng.random((n, 432)) < 0.5).astype(np.uint8)
    truth = dict(_telemetry_truth(n), payload=payload, payload_kind="random")
    return ZooItem("human/telemetry_random", "human", BitStream(telemetry_frames(list(payload))), truth,
                   expect(True, ContentTier.LAYOUT, "human"),
                   "Telemetry frames with random payloads: the envelope without readable content.")


HAMMING_P = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=np.uint8)
HAMMING_G = np.concatenate([np.eye(4, dtype=np.uint8), HAMMING_P], axis=1)        # 4 x 7, systematic
HAMMING_H = np.concatenate([HAMMING_P.T, np.eye(3, dtype=np.uint8)], axis=1)      # 3 x 7


@register("human/hamming74_text")
def hamming74_text(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """ASCII text (most significant bit first) cut into 4-bit blocks, each sent as the 7-bit codeword d G mod 2
    of a systematic Hamming(7,4) code. truth: n, k, G, H, the data bits and the text."""
    text = english_text(rng, _TEXT_CHARS[size])
    data = text_bits(text)
    code = (data.reshape(-1, 4) @ HAMMING_G % 2).astype(np.uint8).ravel()
    truth = {"code": "hamming(7,4)", "n": 7, "k": 4, "generator": HAMMING_G.copy(), "check": HAMMING_H.copy(),
             "systematic": True, "data": data, "text": text, "bit_order": "msb_first", "length": code.size}
    return ZooItem("human/hamming74_text", "human", BitStream(code), truth, expect(True, ContentTier.LAYOUT, "human"),
                   "ASCII text protected by a systematic Hamming(7,4) code.")


PARITY_G = np.concatenate([np.eye(7, dtype=np.uint8), np.ones((7, 1), np.uint8)], axis=1)  # 7 x 8
PARITY_H = np.ones((1, 8), dtype=np.uint8)


@register("human/parity87_text")
def parity87_text(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """ASCII text as 7-bit characters (most significant first) each followed by an even-parity bit: an (8,7)
    code. truth: n, k, G, H, the data bits and the text."""
    text = english_text(rng, _TEXT_CHARS[size])
    data = text_bits(text).reshape(-1, 8)[:, 1:].ravel()      # drop the always-zero top bit of ASCII
    code = (data.reshape(-1, 7) @ PARITY_G % 2).astype(np.uint8).ravel()
    truth = {"code": "even_parity(8,7)", "n": 8, "k": 7, "generator": PARITY_G.copy(), "check": PARITY_H.copy(),
             "systematic": True, "data": data, "char_bits": 7, "text": text, "bit_order": "msb_first",
             "length": code.size}
    return ZooItem("human/parity87_text", "human", BitStream(code), truth, expect(True, ContentTier.LAYOUT, "human"),
                   "Seven-bit ASCII text with an even-parity bit per character.")


LFSR_POLY = (15, 14, 0)  # x^15 + x^14 + 1


def lfsr_sequence(seed: np.ndarray, n: int) -> np.ndarray:
    """The PRBS for x^15 + x^14 + 1: s[0..14] = seed (not all zero), s[i] = s[i-14] XOR s[i-15]. Period 32767."""
    seed = np.asarray(seed, dtype=np.uint8)
    if seed.size != 15 or not seed.any():
        raise ValueError("seed must be 15 bits, not all zero")
    s = np.zeros(max(n, 15), dtype=np.uint8)
    s[:15] = seed
    # The recurrence reaches back 14 steps, so it can be filled in vectorised chunks of 14.
    i = 15
    while i < n:
        j = min(i + 14, n)
        s[i:j] = s[i - 14:j - 14] ^ s[i - 15:j - 15]
        i = j
    return s[:n]


def _lfsr_truth(seed: np.ndarray) -> dict:
    return {"polynomial": LFSR_POLY, "polynomial_text": "x^15 + x^14 + 1", "degree": 15,
            "recurrence": "s[i] = s[i-14] xor s[i-15]", "seed": seed.copy(), "period": 2 ** 15 - 1}


def _lfsr_seed(rng: np.random.Generator) -> np.ndarray:
    while True:
        seed = (rng.random(15) < 0.5).astype(np.uint8)
        if seed.any():
            return seed


@register("human/lfsr_scrambled_text")
def lfsr_scrambled_text(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """ASCII text (most significant bit first) XORed with the x^15 + x^14 + 1 PRBS starting at stream bit 0.
    truth: polynomial, recurrence, seed (the first 15 PRBS bits), the data and the text."""
    text = english_text(rng, _TEXT_CHARS[size])
    data = text_bits(text)
    seed = _lfsr_seed(rng)
    bits = data ^ lfsr_sequence(seed, data.size)
    truth = dict(_lfsr_truth(seed), data=data, text=text, bit_order="msb_first", length=bits.size)
    return ZooItem("human/lfsr_scrambled_text", "human", BitStream(bits), truth,
                   expect(True, ContentTier.LAYOUT, "human"), "ASCII text scrambled by the x^15 + x^14 + 1 PRBS.")


@register("human/prbs15")
def prbs15(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """The bare x^15 + x^14 + 1 PRBS: random-looking to every test except linear complexity. truth: polynomial,
    recurrence, seed."""
    seed = _lfsr_seed(rng)
    n = 20_000 if size == "small" else 2 * (2 ** 15 - 1)
    truth = dict(_lfsr_truth(seed), length=n)
    return ZooItem("human/prbs15", "human", BitStream(lfsr_sequence(seed, n)), truth,
                   expect(True, ContentTier.LAYOUT, "human"), "The bare x^15 + x^14 + 1 pseudo-random bit sequence.")


@register("human/ook_pulses")
def ook_pulses(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """On-off keying as an EventList: slot k spans [k, k+1) s and holds a pulse at its centre (5 ms Gaussian
    jitter) when bit k of the 8N1-framed text is 1. truth: slot length, slot centre offset, jitter, the bits, the
    text, the framing."""
    text = english_text(rng, 60 if size == "small" else 600)
    bits = frame_8n1(text)
    slot, offset, jitter = 1.0, 0.5, 0.005
    index = np.flatnonzero(bits)
    time = index * slot + offset + jitter * rng.standard_normal(index.size)
    events = EventList(time, amplitude=np.ones(index.size), exposure=[(0.0, bits.size * slot)],
                       source_id="human/ook_pulses", time_unit="s")
    truth = {"slot_length": slot, "slot_offset": offset, "jitter": jitter, "bits": bits, "n_slots": int(bits.size),
             "framing": dict(_8N1_TRUTH), "text": text}
    return ZooItem("human/ook_pulses", "human", events, truth, expect(True, ContentTier.ANCHORED, "human"),
                   "On-off-keyed pulses in 1 s slots carrying 8N1-framed ASCII text.")

