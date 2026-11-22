"""Tests for the geometric layout finder (stargazer.decode.layout).

Calibration streams are built here from first principles; recovery and robustness use the benchmark zoo (only the
tests read an item's truth). Every random draw comes from a named, fixed-seed generator.
"""
from __future__ import annotations

import math
import time

import numpy as np
import pytest

from stargazer.bench import build_zoo, erase_bits, flip_bits, impaired, prepend_random, truncate
from stargazer.core.rng import rng_for
from stargazer.core.types import BitStream
from stargazer.decode.layout import (CONTEXTS, baseline_bits, factorisations, fold_scan, folded_bits, kt_bits,
                                     render, to_text, volume_scan)

SEED = 20261122


@pytest.fixture(scope="module")
def zoo() -> dict:
    return {item.name: item for item in build_zoo(seed=0, collections=("message", "control", "human"))}


def _sequential_kt(bits: np.ndarray, contexts: np.ndarray) -> float:
    """Code length of the sequential adaptive code: each bit predicted with (count + 1/2) / (total + 1) in its
    context, counts updated after coding."""
    counts: dict[int, list[int]] = {}
    total = 0.0
    for b, c in zip(bits.tolist(), contexts.tolist()):
        n = counts.setdefault(c, [0, 0])
        total -= math.log2((n[b] + 0.5) / (n[0] + n[1] + 1.0))
        n[b] += 1
    return total


def _markov(rng: np.random.Generator, n: int, stay0: float, stay1: float) -> np.ndarray:
    u = rng.random(n)
    out = np.empty(n, dtype=np.uint8)
    state = int(u[0] < (1 - stay0) / (2 - stay0 - stay1))
    for i in range(n):
        out[i] = state
        if u[i] >= (stay0 if state == 0 else stay1):
            state = 1 - state
    return out


# ------------------------------------------------------------------------------------------------ code lengths

def test_kt_bits_matches_sequential_code() -> None:
    rng = rng_for(SEED, "kt")
    for n, p in [(1, 0.5), (17, 0.5), (300, 0.1), (1000, 0.7)]:
        bits = (rng.random(n) < p).astype(np.uint8)
        direct = _sequential_kt(bits, np.zeros(n, dtype=np.int64))
        assert kt_bits(int((bits == 0).sum()), int(bits.sum())) == pytest.approx(direct, abs=1e-8)
    assert kt_bits(0, 0) == pytest.approx(0.0, abs=1e-12)
    assert kt_bits(1, 0) == pytest.approx(1.0)
    vec = kt_bits(np.array([0, 3, 10]), np.array([0, 5, 2]))
    assert isinstance(vec, np.ndarray) and vec.shape == (3,)


def test_context_models_are_real_sequential_codes() -> None:
    """The bincount totals equal the bit-by-bit adaptive code, for the baseline and for a 2-D context."""
    rng = rng_for(SEED, "seqcode")
    bits = (rng.random(600) < 0.35).astype(np.uint8)
    pad = np.concatenate([np.zeros(16, np.uint8), bits])
    best = min(_sequential_kt(bits, sum(pad[16 - j:16 - j + bits.size].astype(np.int64) << (j - 1)
                                         for j in range(1, m + 1)) if m else np.zeros(bits.size, np.int64))
               for m in range(9))
    assert baseline_bits(bits)[0] == pytest.approx(best + math.log2(9), abs=1e-6)
    w = 7
    sh = {o: pad[16 - o:16 - o + bits.size].astype(np.int64) for o in (1, 2, w - 1, w, w + 1, 2 * w)}
    full = sh[1] | sh[w + 1] << 1 | sh[w] << 2 | sh[w - 1] << 3 | sh[2] << 4 | sh[2 * w] << 5
    costs = []
    for keep in CONTEXTS.values():
        ctx = sum(((full >> k) & 1) << j for j, k in enumerate(keep))
        costs.append(_sequential_kt(bits, ctx))
    assert folded_bits(bits, w) == pytest.approx(min(costs) + 2.0, abs=1e-6)


def test_erased_bits_are_not_coded() -> None:
    rng = rng_for(SEED, "erase")
    bits = (rng.random(500) < 0.5).astype(np.uint8)
    known = rng.random(500) < 0.7
    garbage = np.where(known, bits, 1 - bits)          # what erased positions hold must not matter
    assert baseline_bits(bits, known)[0] == pytest.approx(baseline_bits(garbage, known)[0])
    assert folded_bits(bits, 9, known) == pytest.approx(folded_bits(garbage, 9, known))


# ------------------------------------------------------------------------------------------------ calibration

def test_calibration_fair_and_markov_streams() -> None:
    fair = np.array([fold_scan(BitStream((rng_for(SEED, "fair", i).random(2000) < 0.5).astype(np.uint8))).gain_bits
                     for i in range(200)])
    markov = np.array([fold_scan(BitStream(_markov(rng_for(SEED, "markov", i), 2000, 0.9, 0.7))).gain_bits
                       for i in range(100)])
    for gains in (fair, markov):
        assert gains.max() < 20.0
        assert np.mean(gains >= 5.0) <= 0.03
    # As measured: fair maximum about -8 bits, Markov maximum about -5 bits; none reaches 5 bits.
    assert fair.max() < 5.0 and markov.max() < 5.0


# ------------------------------------------------------------------------------------------------ recovery

def test_raster_23x73(zoo: dict) -> None:
    item = zoo["message/raster_23x73"]
    r = fold_scan(item.payload)
    assert r.width == item.truth["width"] == 23
    assert r.detected and r.gain_bits > 60 and r.p_bound < 1e-18
    assert r.complete and r.n_rows == 73 and r.header_bits == 0 and r.seam == 0
    assert r.factor_note.startswith("1679 = 23 x 73: both prime")
    assert r.gain_at(73) < 0                         # the transposed reading does not fold
    assert np.array_equal(render(item.payload, r.width), item.truth["image"])


def test_raster_127x127(zoo: dict) -> None:
    item = zoo["message/raster_127x127"]
    r = fold_scan(item.payload)
    assert r.width == 127 and r.detected and r.gain_bits > 1000 and r.header_bits == 0


def test_framed_raster_header(zoo: dict) -> None:
    item = zoo["message/raster_23x73_framed"]
    r = fold_scan(item.payload)
    assert r.width == 23 and r.detected
    assert abs(r.header_bits - item.truth["header_len"]) <= 4          # measured: exactly 40
    assert r.seam == item.truth["offset"] % 23                         # 40 mod 23 = 17
    picture = render(item.payload, 23, offset=r.header_bits, height=73)
    assert np.mean(picture == item.truth["image"]) > 0.95


def test_raster_3d(zoo: dict) -> None:
    item = zoo["message/raster3d_11x13x17"]
    r = fold_scan(item.payload)
    assert r.width == 17 and r.detected
    v = volume_scan(item.payload, r.width)
    assert v.height == 13 and v.detected and v.depth_layers == 11 and v.plane_bits == 221


# ------------------------------------------------------------------------------------------------ robustness

@pytest.mark.parametrize("label, channel, kwargs, detectable", [
    ("flip3", flip_bits, {"rate": 0.03}, True),
    ("flip10", flip_bits, {"rate": 0.10}, False),
    ("prepend17", prepend_random, {"n": 17}, True),
    ("trunc1200", truncate, {"n": 1200}, True),
    ("erase20", erase_bits, {"fraction": 0.20}, True),
])
def test_raster_robustness(zoo: dict, label: str, channel, kwargs: dict, detectable: bool) -> None:
    if channel is not truncate:
        kwargs = dict(kwargs, rng=rng_for(SEED, "robust", label))
    item = impaired(zoo["message/raster_23x73"], channel, **kwargs)
    r = fold_scan(item.payload)
    if detectable:
        assert r.detected and r.width == 23, (label, r.width, r.gain_bits)
    else:
        # 10% flips on a 1,679-bit picture: the best gain is a few bits at most; honestly not detectable.
        assert not r.detected and r.gain_bits < 20, (label, r.gain_bits)
    if label == "prepend17":
        assert abs(r.header_bits - 17) <= 3 and abs(r.seam - 17) <= 3      # measured: 16 and 16
    if label == "erase20":
        img = render(item.payload, 23)
        assert np.ma.isMaskedArray(img) and "?" in to_text(img)


# ------------------------------------------------------------------------------------------------ controls

@pytest.mark.parametrize("name", ["control/fair_bits", "control/markov1_bits", "control/zlib_text",
                                  "control/encrypted_text"])
def test_controls_not_detected(zoo: dict, name: str) -> None:
    r = fold_scan(zoo[name].payload)
    assert not r.detected and r.gain_bits < 0 and r.p_bound == 1.0


# ------------------------------------------------------------------------------------------------ human items

def test_human_items(zoo: dict) -> None:
    ascii_r = fold_scan(zoo["human/ascii_8n1"].payload)
    telem_r = fold_scan(zoo["human/telemetry_text"].payload)
    rand_r = fold_scan(zoo["human/telemetry_random"].payload)
    # The best width is the frame length, but against an order-8 baseline that already models the text these
    # do not fold profitably (measured gains about -530 and -3360 bits), so they are not detected.
    assert ascii_r.width % 10 == 0 and not ascii_r.detected
    assert telem_r.width % 496 == 0 and not telem_r.detected
    # Random payloads leave only the framing: the sync word and counter fold at the frame length.
    assert rand_r.width == 496 and rand_r.detected


# ------------------------------------------------------------------------------------------------ utilities

def test_factorisations_and_render() -> None:
    f = factorisations(1679)
    assert f["prime_factors"] == [23, 73] and f["pairs"] == [(23, 73), (73, 23)] and f["semiprime"]
    assert factorisations(97)["prime"] and factorisations(97)["pairs"] == []
    s = BitStream(np.array([1, 0, 0, 1, 1, 1, 0], dtype=np.uint8))
    img = render(s, 3)
    assert img.shape == (2, 3) and to_text(img) == "#..\n###"
    assert render(s, 3, offset=1).tolist() == [[0, 0, 1], [1, 1, 0]]


def test_surrogate_check(zoo: dict) -> None:
    r = fold_scan(zoo["message/raster_23x73"].payload, n_surrogates=9, rng=rng_for(SEED, "surrogates"))
    assert r.surrogate_p == pytest.approx(0.1)
    assert math.isnan(fold_scan(zoo["message/raster_23x73"].payload).surrogate_p)


def test_speed_20k_bits() -> None:
    bits = (rng_for(SEED, "speed").random(20_000) < 0.5).astype(np.uint8)
    t0 = time.perf_counter()
    r = fold_scan(BitStream(bits))
    assert time.perf_counter() - t0 < 5.0          # measured about 0.4 s (4,095 widths)
    assert r.widths.size == 4095
