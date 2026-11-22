"""Tests for the benchmark zoo: determinism, real ground truth, natural statistics, controls and channels."""
from __future__ import annotations

import copy
import math
import time
import zlib

import numpy as np
import pytest

from stargazer.bench import (GENERATORS, ZooItem, add_noise, build_zoo, drop_bits, drop_events, drop_points,
                             erase_bits, flip_bits, impaired, jitter_times, prepend_random, truncate)
from stargazer.bench.controls import sha256_ctr_keystream
from stargazer.bench.natural import simulate_drw
from stargazer.core import BitStream, ContentTier, EventList, Series, SymbolSequence, rng_for


@pytest.fixture(scope="module")
def zoo() -> dict[str, ZooItem]:
    return {item.name: item for item in build_zoo(0)}


# --------------------------------------------------------------------------------------------- helpers

def _equal(a, b) -> bool:
    """Deep equality over dicts, sequences, numpy arrays and scalars."""
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.asarray(a).shape == np.asarray(b).shape and np.array_equal(np.asarray(a), np.asarray(b))
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(_equal(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return isinstance(b, (list, tuple)) and len(a) == len(b) and all(_equal(x, y) for x, y in zip(a, b))
    return a == b


def _payload_arrays(p) -> list[np.ndarray]:
    if isinstance(p, BitStream):
        return [p.bits, np.array([]) if p.known is None else p.known]
    if isinstance(p, Series):
        return [p.time, p.value, p.error]
    if isinstance(p, EventList):
        return [p.time] + [x for x in (p.amplitude, p.time_error) if x is not None]
    if isinstance(p, SymbolSequence):
        return [p.symbols]
    raise TypeError(type(p))


def _monobit_p(bits: np.ndarray) -> float:
    b = bits.astype(np.int64)
    s = 2 * int(b.sum()) - b.size
    return math.erfc(abs(s) / math.sqrt(2 * b.size))


def _runs_p(bits: np.ndarray) -> float:
    b = bits.astype(np.int64)
    n, pi = b.size, float(b.mean())
    runs = 1 + int(np.count_nonzero(np.diff(b)))
    return math.erfc(abs(runs - 2 * n * pi * (1 - pi)) / (2 * math.sqrt(2 * n) * pi * (1 - pi)))


def _bits_to_int(bits) -> int:
    out = 0
    for b in np.asarray(bits, dtype=int):
        out = 2 * out + int(b)
    return out


def _crc16(bits) -> int:
    crc = 0xFFFF
    for b in np.asarray(bits, dtype=int):
        top = ((crc >> 15) & 1) ^ b
        crc = ((crc << 1) & 0xFFFF) ^ (0x1021 if top else 0)
    return crc


def _ascii_msb(bits: np.ndarray) -> str:
    return np.packbits(np.asarray(bits, dtype=np.uint8)).tobytes().decode("ascii")


def _decode_8n1(bits: np.ndarray) -> str:
    frames = np.asarray(bits).reshape(-1, 10)
    assert np.all(frames[:, 0] == 0) and np.all(frames[:, 9] == 1)
    return np.packbits(frames[:, 1:9], axis=1, bitorder="little").ravel().tobytes().decode("ascii")


def _lfsr(seed: np.ndarray, n: int) -> np.ndarray:
    s = list(int(x) for x in seed)
    while len(s) < n:
        s.append(s[-14] ^ s[-15])
    return np.array(s[:n], dtype=np.uint8)


def _unary_runs(bits: np.ndarray) -> list[int]:
    return [len(run) for run in "".join(map(str, bits.tolist())).split("0")]


# --------------------------------------------------------------------------------------------- the zoo itself

def test_build_is_deterministic_unique_and_fast():
    start = time.perf_counter()
    first = build_zoo(0)
    elapsed = time.perf_counter() - start
    second = build_zoo(0)
    assert elapsed < 10.0
    names = [item.name for item in first]
    assert len(names) == len(set(names))
    assert set(names) == set(GENERATORS)
    for a, b in zip(first, second):
        assert a.name == b.name and type(a.payload) is type(b.payload)
        for x, y in zip(_payload_arrays(a.payload), _payload_arrays(b.payload)):
            assert np.array_equal(x, y), a.name
        assert _equal(a.truth, b.truth), a.name
        assert a.expected == b.expected
    other = {item.name: item for item in build_zoo(1)}
    seed0_bits = first[names.index("control/fair_bits")].payload.bits
    assert not np.array_equal(other["control/fair_bits"].payload.bits, seed0_bits)


def test_expected_keys_and_collection_rules(zoo):
    for item in zoo.values():
        assert {"structured", "max_content_tier", "origin"} <= set(item.expected), item.name
        assert item.name.startswith(item.collection + "/")
        assert isinstance(item.expected["max_content_tier"], int)
        assert item.description and "\n" not in item.description
        if item.collection == "control":
            assert item.expected["structured"] is False, item.name
            assert item.expected["max_content_tier"] == ContentTier.NONE
        if item.collection == "natural":
            assert item.expected["structured"] is True and item.expected["origin"] == "natural"
            assert item.expected["max_content_tier"] <= 2 and item.truth["fair_null"], item.name
        if item.collection in ("human", "message"):
            assert item.expected["structured"] is True
            assert item.expected["origin"] == ("human" if item.collection == "human" else "designed")
    assert zoo["control/markov1_bits"].truth["fair_null"] == "markov1"
    assert zoo["natural/pulsar_nulling_bits"].truth["fair_null"] == "markov1"


def test_collection_filter_and_sizes():
    msgs = build_zoo(0, collections=("message",))
    assert msgs and all(item.collection == "message" for item in msgs)
    for item in build_zoo(0):
        if isinstance(item.payload, BitStream):
            assert len(item.payload) <= 36_000, item.name
    with pytest.raises(ValueError):
        build_zoo(0, size="huge")


# --------------------------------------------------------------------------------------------- human ground truth

def test_8n1_decodes(zoo):
    item = zoo["human/ascii_8n1"]
    assert item.truth["frame_len"] == 10
    assert _decode_8n1(item.payload.bits) == item.truth["text"]


def test_manchester_decodes(zoo):
    item = zoo["human/ascii_manchester"]
    pairs = item.payload.bits.reshape(-1, 2)
    assert np.all(pairs[:, 0] != pairs[:, 1])
    assert _ascii_msb(pairs[:, 0]) == item.truth["text"]


@pytest.mark.parametrize("name", ["human/telemetry_text", "human/telemetry_random"])
def test_telemetry_frames_and_crc(zoo, name):
    item = zoo[name]
    t = item.truth
    frames = item.payload.bits.reshape(-1, t["frame_len"])
    f = t["fields"]
    payloads = []
    for i, frame in enumerate(frames):
        assert _bits_to_int(frame[slice(*f["sync"])]) == t["sync_word"] == 0x1ACFFC1D
        assert _bits_to_int(frame[slice(*f["counter"])]) == i
        assert _crc16(frame[slice(*t["crc"]["covers"])]) == _bits_to_int(frame[slice(*f["crc"])])
        payloads.append(frame[slice(*f["payload"])])
    if t["payload_kind"] == "ascii_text":
        assert _ascii_msb(np.concatenate(payloads)) == t["text"]
    else:
        assert np.array_equal(np.array(payloads), t["payload"])


def test_hamming_corrects_single_errors(zoo):
    item = zoo["human/hamming74_text"]
    g, h = item.truth["generator"], item.truth["check"]
    assert (g @ h.T % 2 == 0).all()
    blocks = item.payload.bits.reshape(-1, item.truth["n"]).copy()
    rng = np.random.default_rng(11)
    blocks[np.arange(len(blocks)), rng.integers(0, 7, len(blocks))] ^= 1
    syndrome = blocks @ h.T % 2
    for row, s in zip(blocks, syndrome):
        pos = int(np.flatnonzero((h.T == s).all(axis=1))[0])
        row[pos] ^= 1
    data = blocks[:, :item.truth["k"]].ravel()
    assert np.array_equal(data, item.truth["data"])
    assert _ascii_msb(data) == item.truth["text"]


def test_parity_code(zoo):
    item = zoo["human/parity87_text"]
    blocks = item.payload.bits.reshape(-1, 8)
    assert np.all(blocks.sum(axis=1) % 2 == 0)
    chars = np.concatenate([np.zeros((len(blocks), 1), np.uint8), blocks[:, :7]], axis=1)
    assert _ascii_msb(chars.ravel()) == item.truth["text"]


def test_lfsr_descrambles(zoo):
    item = zoo["human/lfsr_scrambled_text"]
    t = item.truth
    assert t["polynomial"] == (15, 14, 0)
    clear = item.payload.bits ^ _lfsr(t["seed"], len(item.payload))
    assert _ascii_msb(clear) == t["text"]
    bare = zoo["human/prbs15"]
    assert np.array_equal(bare.payload.bits, _lfsr(bare.truth["seed"], len(bare.payload)))


def test_ook_slots_decode(zoo):
    item = zoo["human/ook_pulses"]
    t = item.truth
    slots = np.floor(item.payload.time / t["slot_length"]).astype(int)
    bits = np.zeros(t["n_slots"], dtype=np.uint8)
    bits[slots] = 1
    assert np.array_equal(bits, t["bits"])
    assert _decode_8n1(bits) == t["text"]


# --------------------------------------------------------------------------------------------- message ground truth

def _check_raster_elements(img: np.ndarray, elements: dict) -> None:
    for key in ("numbers_1_to_10", "atomic_numbers"):
        e = elements[key]
        r0, r1 = e["rows"]
        assert [_bits_to_int(img[r0:r1, c]) for c in e["columns"]] == e["values"]
        assert all(img[e["marker_row"], c] == 1 for c in e["columns"])


@pytest.mark.parametrize("name", ["message/raster_23x73", "message/raster_127x127"])
def test_rasters_reshape(zoo, name):
    item = zoo[name]
    t = item.truth
    img = item.payload.bits.reshape(t["height"], t["width"])
    assert np.array_equal(img, t["image"])
    if name.endswith("23x73"):
        assert (t["width"], t["height"]) == (23, 73)
        _check_raster_elements(img, t["elements"])
        for e in t["elements"].values():
            (r0, r1), (c0, c1) = e["rows"], e["cols"]
            assert img[r0:r1, c0:c1].any()
    else:
        font = t["elements"]["text"]["font"]
        for line in t["elements"]["text"]["lines"]:
            for i, ch in enumerate(line["text"]):
                c = line["col"] + i * 6
                assert np.array_equal(img[line["row"]:line["row"] + 7, c:c + 5], font[ch])


def test_framed_raster(zoo):
    item = zoo["message/raster_23x73_framed"]
    t = item.truth
    bits = item.payload.bits
    assert t["header_len"] == 40 and t["footer_len"] == 24
    assert np.array_equal(bits[-24:], np.tile([1, 0], 12))
    img = bits[t["header_len"]:len(bits) - t["footer_len"]].reshape(t["height"], t["width"])
    assert np.array_equal(img, zoo["message/raster_23x73"].truth["image"])


def test_raster3d(zoo):
    item = zoo["message/raster3d_11x13x17"]
    vol = item.payload.bits.reshape(item.truth["dims"])
    assert np.array_equal(vol, item.truth["volume"])
    assert vol[0, 0, :].all() and vol[-1, -1, :].all() and vol[5, 6, 8] == 1 and vol[0, 6, 8] == 0


def test_unary_sequences(zoo):
    count = zoo["message/unary_count_1_40"]
    assert _unary_runs(count.payload.bits) == count.truth["sequence"] == list(range(1, 41))
    primes = zoo["message/unary_primes_300"]
    runs = _unary_runs(primes.payload.bits)
    assert runs == primes.truth["sequence"]
    assert runs == [n for n in range(2, 300) if all(n % d for d in range(2, int(n ** 0.5) + 1))]


def test_binary_count(zoo):
    item = zoo["message/binary_count_255"]
    t = item.truth
    w = t["field_width"]
    bits = np.concatenate([item.payload.bits, np.zeros(w, np.uint8)]).reshape(-1, 2 * w)
    assert [_bits_to_int(r[:w]) for r in bits] == list(range(1, 256)) == t["sequence"]
    assert all(_bits_to_int(r[w:]) == t["delimiter"] for r in bits)
    assert t["delimiter"] not in t["sequence"]


def test_arithmetic_primer_all_true(zoo):
    item = zoo["message/arithmetic_primer"]
    t = item.truth
    tokens = [_bits_to_int(r) for r in item.payload.bits.reshape(-1, t["token_width"])]
    decode = {v: k for k, v in t["code_book"].items()}
    statements = []
    for i in range(0, len(tokens), 6):
        a, op, b, eq, c, end = tokens[i:i + 6]
        assert decode[eq] == "=" and decode[end] == ";"
        sym = decode[op]
        assert c == (a + b if sym == "+" else a * b)
        statements.append((a, sym, b, c))
    assert statements == t["statements"] and len(statements) == 400


# --------------------------------------------------------------------------------------------- natural statistics

def test_drw_lag_one_correlation_even_sampling():
    tau, dt, n = 300.0, 10.0, 5_000
    x = simulate_drw(np.arange(n) * dt, 0.15, tau, rng_for(3, "drw-test"))
    a = math.exp(-dt / tau)
    r = float(np.corrcoef(x[:-1], x[1:])[0, 1])
    # Sampling error of the lag-one estimate of an AR(1) process, with a small allowance for its O(1/n) bias.
    assert abs(r - a) < 4 * math.sqrt((1 - a ** 2) / n) + 5.0 / n


def test_drw_item_innovations_are_white(zoo):
    item = zoo["natural/drw_quasar"]
    t = item.truth
    x = t["model"] - t["mean_mag"]
    a = np.exp(-np.diff(item.payload.time) / t["tau"])
    z = (x[1:] - a * x[:-1]) / (t["sigma"] * np.sqrt(1 - a ** 2))
    n = z.size
    assert abs(z.mean()) < 4 / math.sqrt(n)
    assert abs(z.var() - 1) < 4 * math.sqrt(2 / n)
    seasons = item.payload.seasons()
    assert seasons.max() + 1 == 7
    assert all(30 <= np.sum(seasons == s) <= 45 for s in range(7))
    assert np.all((item.payload.error >= 0.01) & (item.payload.error <= 0.03))


def test_pulsar_transitions_match_matrix(zoo):
    item = zoo["natural/pulsar_nulling"]
    t = item.truth
    s = t["on_off"].astype(int)
    for state, p_stay in ((1, 0.9), (0, 0.7)):
        frm = s[:-1] == state
        n = int(frm.sum())
        p_hat = float(np.mean(s[1:][frm] == state))
        assert abs(p_hat - p_stay) < 3 * math.sqrt(p_stay * (1 - p_stay) / n)
    assert np.allclose(t["transition_matrix"], [[0.7, 0.3], [0.1, 0.9]])
    assert len(item.payload) == s.sum()
    residual = item.payload.time - t["t0"] - t["pulse_index"] * t["period"]
    assert np.abs(residual).max() < 6 * t["jitter"] and abs(residual.std() - t["jitter"]) < 0.1 * t["jitter"]
    bits = zoo["natural/pulsar_nulling_bits"]
    assert np.array_equal(bits.payload.bits, bits.truth["on_off"])


def test_burster_inside_exposure_and_windows(zoo):
    item = zoo["natural/repeating_burster"]
    t = item.truth
    times = item.payload.time
    assert len(times) == t["n_bursts"] == 400
    exposure = np.array(item.payload.exposure)
    inside = ((times[:, None] >= exposure[:, 0]) & (times[:, None] <= exposure[:, 1])).any(axis=1)
    assert inside.all()
    phase = (times - t["window_start"]) % t["period"]
    assert np.all((phase <= t["window"] + 1e-9) | (phase >= t["period"] - 1e-9))
    covered = float(np.sum(exposure[:, 1] - exposure[:, 0])) / t["duration"]
    assert abs(covered - 0.6) < 1e-9
    # Clustered: within one active stretch, waiting times have a coefficient of variation above the Poisson value 1.
    gaps = np.diff(times)
    within = gaps[gaps < 1.0]
    assert t["shape"] == 0.6 and within.std() / within.mean() > 1.0


def test_eclipsing_binary_truth(zoo):
    item = zoo["natural/eclipsing_binary"]
    t = item.truth
    assert len(item.payload) == 90 * 48
    frac = float(np.mean(t["in_eclipse"]))
    assert abs(frac - t["duration"] / t["period"]) < 0.005
    dip = item.payload.value[~t["in_eclipse"]].mean() - item.payload.value[t["in_eclipse"]].mean()
    assert abs(dip - t["depth"]) < 0.002


# --------------------------------------------------------------------------------------------- controls

@pytest.mark.parametrize("name", ["control/zlib_text", "control/encrypted_text"])
def test_compressed_and_encrypted_look_random(zoo, name):
    bits = zoo[name].payload.bits
    assert _monobit_p(bits) > 0.01
    assert _runs_p(bits) > 0.01
    raw = np.packbits(bits).tobytes()
    assert len(zlib.compress(raw, 9)) >= 0.99 * len(raw)


def test_controls_carry_their_content(zoo):
    z = zoo["control/zlib_text"]
    body = np.packbits(z.payload.bits).tobytes()
    assert zlib.decompress(body, -15).decode("ascii") == z.truth["text"]
    e = zoo["control/encrypted_text"]
    cipher = np.packbits(e.payload.bits).tobytes()
    stream = sha256_ctr_keystream(bytes.fromhex(e.truth["key"]), len(cipher))
    assert bytes(c ^ k for c, k in zip(cipher, stream)).decode("ascii") == e.truth["text"]
    biased = zoo["control/biased_bits_p30"].payload.bits
    assert abs(biased.mean() - 0.3) < 4 * math.sqrt(0.21 / biased.size)


# --------------------------------------------------------------------------------------------- channels

def test_flip_and_erase_rates(zoo):
    stream = zoo["control/fair_bits"].payload
    n = len(stream)
    flipped = flip_bits(stream, 0.03, np.random.default_rng(1))
    rate = float(np.mean(flipped.bits != stream.bits))
    assert abs(rate - 0.03) < 4 * math.sqrt(0.03 * 0.97 / n)
    erased = erase_bits(stream, 0.1, np.random.default_rng(2))
    assert abs(float(np.mean(~erased.known)) - 0.1) < 4 * math.sqrt(0.09 / n)
    assert stream.known is None


def test_impaired_copies_and_updates_truth(zoo):
    item = zoo["message/raster_23x73"]
    before_bits = item.payload.bits.copy()
    before_truth = copy.deepcopy(item.truth)
    noisy = impaired(item, flip_bits, rate=0.03, rng=np.random.default_rng(5))
    assert noisy.name == "message/raster_23x73+flip0.03"
    assert noisy.truth["channels"][-1]["channel"] == "flip_bits"
    shifted = impaired(noisy, prepend_random, n=17, rng=np.random.default_rng(6))
    assert shifted.name == "message/raster_23x73+flip0.03+prepend17"
    assert shifted.truth["offset"] == 17 and shifted.truth["length"] == 1679 + 17
    assert np.array_equal(shifted.payload.bits[17:], noisy.payload.bits)
    cut = impaired(item, truncate, n=1000)
    assert cut.truth["length"] == 1000 and len(cut.payload) == 1000
    dropped = impaired(item, drop_bits, fraction=0.01, rng=np.random.default_rng(7))
    assert len(dropped.payload) == 1679 - 17 and len(dropped.truth["channels"][0]["deleted"]) == 17
    assert np.array_equal(item.payload.bits, before_bits)
    assert _equal(item.truth, before_truth) and "channels" not in item.truth
    assert item.name == "message/raster_23x73" and item.payload.meta == {}


def test_series_and_event_channels(zoo):
    series = zoo["natural/drw_quasar"].payload
    noisy = add_noise(series, 0.05, np.random.default_rng(1))
    assert len(noisy) == len(series) and np.all(noisy.error > series.error)
    fewer = drop_points(series, 0.2, np.random.default_rng(2))
    assert len(fewer) == len(series) - round(0.2 * len(series))
    events = zoo["natural/pulsar_nulling"].payload
    before = events.time.copy()
    jittered = jitter_times(events, 0.01, np.random.default_rng(3))
    assert len(jittered) == len(events) and not np.array_equal(jittered.time, events.time)
    thinned = drop_events(events, 0.5, np.random.default_rng(4))
    assert len(thinned) == len(events) - round(0.5 * len(events))
    assert np.array_equal(events.time, before)
    item = impaired(zoo["natural/pulsar_nulling"], drop_events, fraction=0.5, rng=np.random.default_rng(4))
    assert item.name.endswith("+dropev0.5") and item.truth["n_points"] == len(thinned)
