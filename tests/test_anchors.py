"""Tests for the anchor search (stargazer.decode.anchors).

Each anchor must find its benchmark target with a corrected p-value of at most 1e-6 and the right parameters
(checked against the item's truth here, never inside the module). The false-alarm calibration runs the whole
library on 300 fair and 100 Markov streams and on the zoo controls; the cross-talk test prints which anchors fire
on which human and message items. The heavy searches run once, in a small process pool, and are shared by the
tests. Run with -s to see the printed calibration numbers, cross-talk table and timings.
"""
from __future__ import annotations

import multiprocessing as mp
import time
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache

import numpy as np
import pytest

from stargazer.bench import build_zoo, flip_bits, impaired, prepend_random
from stargazer.bench.controls import markov_bits
from stargazer.bench.zoo import english_text, text_bits
from stargazer.core import BitStream, rng_for
from stargazer.decode import anchors as A
from stargazer.decode.anchors import LIBRARY, AnchorHit, AnchorReport, search

SEED = 20261122
ALPHA = 1e-3
TARGET_P = 1e-6
N_FAIR, N_MARKOV, CAL_BITS = 300, 100, 4_000
MARKOV_T = np.array([[0.9, 0.1], [0.3, 0.7]])          # stay probabilities 0.9 (zeros) and 0.7 (ones)


# ------------------------------------------------------------------------------------------------ shared data

@lru_cache(maxsize=1)
def _zoo() -> dict:
    return {it.name: it for it in build_zoo(SEED) if isinstance(it.payload, BitStream)}


def _calibration_stream(kind: str, i: int) -> BitStream:
    rng = rng_for(SEED, "anchors-calibration", kind, i)
    if kind == "fair":
        return BitStream((rng.random(CAL_BITS) < 0.5).astype(np.uint8))
    return BitStream(markov_bits(rng, CAL_BITS, MARKOV_T))


def _job(job: tuple[str, str, int]) -> tuple[str, list[tuple[str, float, float]]]:
    """Run search on one stream; return (label, [(anchor, corrected p, log10 p), ...])."""
    kind, name, i = job
    stream = _zoo()[name].payload if kind == "zoo" else _calibration_stream(kind, i)
    rep = search(stream, rng_for(SEED, "anchors-search", kind, name, i))
    label = name if kind == "zoo" else f"{kind}/{i}"
    return label, [(h.anchor, pc, h.evidence["log10_p"]) for h, pc in zip(rep.hits, rep.p_corrected)]


@lru_cache(maxsize=1)
def _batch() -> dict[str, list[tuple[str, float, float]]]:
    """Every heavy search, once: calibration streams and all bit-stream zoo items."""
    jobs = [("fair", "", i) for i in range(N_FAIR)] + [("markov", "", i) for i in range(N_MARKOV)]
    jobs += [("zoo", name, 0) for name in sorted(_zoo())]
    t0 = time.perf_counter()
    try:
        with ProcessPoolExecutor(max_workers=4, mp_context=mp.get_context("forkserver")) as pool:
            out = dict(pool.map(_job, jobs, chunksize=8))
    except (OSError, RuntimeError):                     # pragma: no cover - no process pool available
        out = dict(map(_job, jobs))
    print(f"\n[anchors] {len(jobs)} searches in {time.perf_counter() - t0:.1f} s (4 processes)")
    return out


def _hit(report: AnchorReport, anchor: str) -> tuple[AnchorHit, float]:
    for h, pc in zip(report.hits, report.p_corrected):
        if h.anchor == anchor:
            return h, pc
    raise AssertionError(f"{anchor} missing from report")


# ------------------------------------------------------------------------------------------------ the library

def test_library_is_declared_and_always_reports():
    assert list(LIBRARY) == ["run_length_sequence", "field_counting", "constants", "text_code",
                             "functional_relations", "repeated_block"]
    rng = rng_for(SEED, "tiny")
    for n in (0, 5, 40, 300):
        rep = search(BitStream((rng.random(n) < 0.5).astype(np.uint8)), rng_for(SEED, "tiny-search", n))
        assert rep.n_anchors == 6 and len(rep.hits) == 6 and rep.best is None
        assert sorted(h.anchor for h in rep.hits) == sorted(LIBRARY)
        for h, pc in zip(rep.hits, rep.p_corrected):
            assert isinstance(h, AnchorHit) and h.sentence and 0 < h.p_value <= 1
            assert pc == pytest.approx(min(1.0, 6 * h.p_value))
            assert h.family in ("mathematics", "human encoding")
    sub = search(_zoo()["message/unary_count_1_40"].payload, anchors=["run_length_sequence", "constants"])
    assert sub.n_anchors == 2 and sub.best is not None and sub.best.anchor == "run_length_sequence"
    with pytest.raises(ValueError):
        search(BitStream(np.zeros(10, np.uint8)), anchors=["astrology"])


def test_erased_bits_are_read_as_zero():
    item = _zoo()["message/unary_count_1_40"]
    bits = item.payload.bits.copy()
    known = np.ones(bits.size, dtype=bool)
    known[-50:] = False                                  # erase the tail; the anchor must still run
    rep = search(BitStream(bits, known), anchors=["run_length_sequence"])
    assert rep.best is not None and rep.best.parameters["length"] >= 35


# ------------------------------------------------------------------------------------------------ targets

def test_run_length_sequence_finds_unary_counting_and_primes():
    for name, sequence in (("message/unary_count_1_40", "naturals"), ("message/unary_primes_300", "primes")):
        item = _zoo()[name]
        rep = search(item.payload, rng_for(SEED, name))
        h, pc = _hit(rep, "run_length_sequence")
        assert pc <= TARGET_P and rep.best is h
        truth = item.truth["sequence"]
        assert h.parameters["sequence"] == sequence and h.parameters["polarity"] == "ones"
        assert h.parameters["length"] == len(truth) and h.parameters["first_term"] == truth[0]
        assert list(h.evidence["run_lengths"]) == truth
        assert h.evidence["permutation_fraction"] == 0.0
        print(f"\n[run_length_sequence] {name}: corrected p {pc:.2e}; {h.sentence}")


def test_field_counting_finds_binary_count():
    item = _zoo()["message/binary_count_255"]
    rep = search(item.payload, rng_for(SEED, "count"))
    h, pc = _hit(rep, "field_counting")
    assert pc <= TARGET_P
    p = h.parameters
    assert p["bit_order"] == item.truth["bit_order"]
    # Either the truth's reading (8-bit numbers with a delimiter field between: stride 2) or the equivalent 16-bit
    # reading whose fields are "delimiter, number" (stride 1, starting half-way into a number pair).
    assert (p["field_bits"], p["stride"], p["offset"] % p["field_bits"]) in {(8, 2, 0), (16, 1, 8)}
    stretch = list(h.evidence["stretch_values"])
    seq = item.truth["sequence"]
    assert h.evidence["longest_stretch_steps"] >= 250
    assert stretch[0] in seq and stretch[-1] == seq[-1]
    print(f"\n[field_counting] corrected p {pc:.2e}; {h.sentence}")


def test_constants_finds_embedded_pi():
    pi = A._constant_expansions()["pi"][1]
    assert "".join(map(str, pi[:12])) == "110010010000"                  # 11.0010010000... = pi
    assert "".join(map(str, A._constant_expansions()["e"][1][:12])) == "101011011111"     # 10.1011011111...
    assert "".join(map(str, A._constant_expansions()["sqrt2"][1][:12])) == "101101010000"  # 1.01101010000...
    assert "".join(map(str, A._constant_expansions()["golden_ratio"][1][:12])) == "110011110001"  # 1.10011110001
    rng = rng_for(SEED, "pi-item")
    for start_in_pi in (0, int(rng.integers(100, 15_000))):
        bits = (rng.random(20_000) < 0.5).astype(np.uint8)
        pos = int(rng.integers(0, 19_000))
        bits[pos:pos + 300] = pi[start_in_pi:start_in_pi + 300]
        rep = search(BitStream(bits), rng_for(SEED, "pi-search", start_in_pi))
        h, pc = _hit(rep, "constants")
        assert pc <= TARGET_P and rep.best is h
        p = h.parameters
        assert p["constant"] == "pi" and p["length"] >= 300
        # The match may extend by chance by a few bits on either side; its diagonal must be the embedding's.
        assert p["stream_position"] - p["constant_position"] == pos - start_in_pi
        assert p["stream_position"] <= pos and p["stream_position"] + p["length"] >= pos + 300
        assert p["natural_reading_start"] == (start_in_pi == 0 and p["stream_position"] == pos)
        print(f"\n[constants] corrected p {pc:.2e}; {h.sentence}")


def test_text_code_finds_serial_and_plain_text():
    item = _zoo()["human/ascii_8n1"]
    rep = search(item.payload, rng_for(SEED, "8n1"))
    h, pc = _hit(rep, "text_code")
    assert pc <= TARGET_P and rep.best is h
    assert h.parameters == {"reading": "8n1", "offset": 0, "bit_order": "lsb_first"}
    assert h.evidence["preview"][:60] == item.truth["text"][:60]
    assert h.evidence["space_most_common"]
    rng = rng_for(SEED, "plain-text")
    text = english_text(rng, 400)
    for shift in (0, 3):
        bits = np.concatenate([(rng.random(shift) < 0.5).astype(np.uint8), text_bits(text)])
        rep = search(BitStream(bits), rng_for(SEED, "plain", shift))
        h, pc = _hit(rep, "text_code")
        assert pc <= TARGET_P and rep.best is h
        assert h.parameters == {"reading": "8bit", "offset": shift, "bit_order": "msb_first"}
        assert h.evidence["preview"][:40] == text[:40]
        print(f"\n[text_code] plain text shifted {shift}: corrected p {pc:.2e}; {h.sentence}")


def test_text_code_on_hidden_text_is_reported():
    """Text under a line code, an error-correcting code or a scrambler. The scrambled text must not read as text;
    Manchester and Hamming(7,4) leave enough of each character in place that a reading passes (see the report)."""
    for name in ("human/ascii_manchester", "human/hamming74_text", "human/lfsr_scrambled_text"):
        h, pc = _hit(search(_zoo()[name].payload, rng_for(SEED, name), anchors=["text_code"]), "text_code")
        print(f"\n[text_code] {name}: corrected p {pc:.2e}; {h.sentence} preview {h.evidence['preview'][:40]!r}")
        if name == "human/lfsr_scrambled_text":
            assert pc > ALPHA


def _primer_reading_ok(reading: dict, truth: dict, offset: int) -> None:
    codes = truth["code_book"]
    assert reading["field_bits"] == truth["token_width"] == 8
    assert reading["bit_order"] == truth["bit_order"]
    assert reading["record_fields"] == len(truth["statement_layout"]) == 6
    assert (reading["bit_offset"] - offset) % 48 == 0
    assert reading["operator_column"] == 1
    assert reading["operators"] == {codes["+"]: "+", codes["x"]: "*"}
    assert reading["delimiters"] == {3: codes["="], 5: codes[";"]}
    assert reading["operands"] == (0, 2) and reading["result"] == 4


def test_functional_relations_reads_the_primer():
    base = _zoo()["message/arithmetic_primer"]
    cases = [(base, 1.0), (impaired(base, prepend_random, n=13, rng=rng_for(SEED, "prepend")), 1.0),
             (impaired(base, flip_bits, rate=0.002, rng=rng_for(SEED, "flip")), 0.85),
             (impaired(impaired(base, prepend_random, n=13, rng=rng_for(SEED, "prepend2")), flip_bits, rate=0.002,
                       rng=rng_for(SEED, "flip2")), 0.85)]
    for item, min_held in cases:
        t0 = time.perf_counter()
        rep = search(item.payload, rng_for(SEED, item.name), anchors=["functional_relations"])
        h, pc = _hit(rep, "functional_relations")
        dt = time.perf_counter() - t0
        assert pc <= TARGET_P
        assert h.reading is not None
        _primer_reading_ok(h.reading, item.truth, int(item.truth.get("offset", 0)))
        held = h.reading["held_out"]
        assert held["reading_found_in_first_half"] and held["n"] >= 150
        assert held["correct"] >= min_held * held["n"]
        print(f"\n[functional_relations] {item.name} ({dt:.1f} s): corrected p {pc:.2e}; held out "
              f"{held['correct']} of {held['n']}; {h.sentence}")
    # The example sentence of the brief, with this item's own counts.
    h, _ = _hit(search(base.payload, anchors=["functional_relations"]), "functional_relations")
    assert h.sentence.startswith("In records of six 8-bit fields, field 5 equals field 1 plus field 3 whenever "
                                 "field 2 is 252 (")
    assert "and field 1 times field 3 whenever field 2 is 253 (" in h.sentence


def test_repeated_block_finds_the_beacon():
    rng = rng_for(SEED, "beacon")
    block = (rng.random(300) < 0.5).astype(np.uint8)
    bits = np.tile(block, 12) ^ (rng.random(3_600) < 0.02).astype(np.uint8)
    rep = search(BitStream(bits), rng_for(SEED, "beacon-search"))
    h, pc = _hit(rep, "repeated_block")
    assert pc <= TARGET_P and rep.best is h
    assert h.parameters["block_length"] == 300
    assert h.evidence["full_repeats"] == 12
    assert np.array_equal(h.evidence["consensus_block"], block)
    assert 0.94 <= h.evidence["agreement"] <= 0.98
    print(f"\n[repeated_block] corrected p {pc:.2e}; {h.sentence}")


# ------------------------------------------------------------------------------------------------ calibration

def test_false_alarm_calibration():
    """No fair, Markov or control stream may reach `best`, and few may even come near 0.05."""
    batch = _batch()
    controls = [n for n in _zoo() if n.startswith("control/")]
    groups = {"fair": [f"fair/{i}" for i in range(N_FAIR)], "markov": [f"markov/{i}" for i in range(N_MARKOV)],
              "controls": controls}
    all_min = []
    for group, labels in groups.items():
        mins = np.array([batch[lab][0][1] for lab in labels])
        all_min.extend(mins.tolist())
        q = np.quantile(mins, [0.0, 0.01, 0.05, 0.25, 0.5])
        print(f"\n[calibration] {group:8s} n={mins.size:3d}  smallest corrected p: min {q[0]:.3g}, 1% {q[1]:.3g}, "
              f"5% {q[2]:.3g}, 25% {q[3]:.3g}, median {q[4]:.3g};  fraction < 0.05: {np.mean(mins < 0.05):.3f};  "
              f"fraction <= alpha: {np.mean(mins <= ALPHA):.3f}")
        culprits = {}
        for lab in labels:
            anchor, pc, _ = batch[lab][0]
            if pc < 0.2:
                culprits[anchor] = culprits.get(anchor, 0) + 1
        print(f"[calibration] {group:8s} anchors giving the smallest p when below 0.2: {culprits}")
        assert np.all(mins > ALPHA), f"{group}: a stream reached best"
    all_min = np.array(all_min)
    assert np.mean(all_min < 0.05) <= 0.05


def test_cross_talk_table():
    """Which anchors fire (corrected p <= alpha) on which human and message items. Rasters must not fire the
    arithmetic or counting anchors; every target fires its own anchor."""
    batch = _batch()
    names = [n for n in sorted(_zoo()) if n.startswith(("human/", "message/", "natural/"))]
    short = {"run_length_sequence": "runs", "field_counting": "count", "constants": "const", "text_code": "text",
             "functional_relations": "arith", "repeated_block": "repeat"}
    lines = ["\n[cross-talk] log10 corrected p (blank: > 1e-3; '*' marks anchors that fire)",
             f"{'item':32s}" + "".join(f"{s:>9s}" for s in short.values())]
    for name in names:
        res = {a: pc for a, pc, _ in batch[name]}
        cells = []
        for a in LIBRARY:
            pc = res[a]
            cells.append(f"{max(np.log10(pc), -300):8.0f}*" if pc <= ALPHA else " " * 9)
        lines.append(f"{name:32s}" + "".join(cells))
    print("\n".join(lines))
    for name in names:
        res = {a: pc for a, pc, _ in batch[name]}
        if "raster" in name:
            assert res["field_counting"] > ALPHA and res["functional_relations"] > ALPHA, name
    targets = {"message/unary_count_1_40": "run_length_sequence",
               "message/unary_primes_300": "run_length_sequence", "message/binary_count_255": "field_counting",
               "message/arithmetic_primer": "functional_relations", "human/ascii_8n1": "text_code"}
    for name, anchor in targets.items():
        assert {a: pc for a, pc, _ in batch[name]}[anchor] <= TARGET_P, name
    # Natural bits and text without arithmetic must not claim arithmetic or counting.
    for name in ("natural/pulsar_nulling_bits", "human/telemetry_text", "human/telemetry_random",
                 "control/zlib_text"):
        res = {a: pc for a, pc, _ in batch[name]}
        assert res["functional_relations"] > ALPHA and res["field_counting"] > ALPHA, name


def test_timing():
    rng = rng_for(SEED, "timing")
    stream = BitStream((rng.random(20_000) < 0.5).astype(np.uint8))
    search(stream)                                       # warm the caches (constants, tables)
    t0 = time.perf_counter()
    rep = search(stream)
    dt = time.perf_counter() - t0
    per = {}
    for name, fn in LIBRARY.items():
        t1 = time.perf_counter()
        fn(stream, np.random.default_rng(0))
        per[name] = round(time.perf_counter() - t1, 3)
    print(f"\n[timing] search on 20,000 random bits: {dt:.2f} s; per anchor {per}")
    assert rep.best is None
    assert dt < 8.0                                      # a few seconds, with headroom for a loaded machine
