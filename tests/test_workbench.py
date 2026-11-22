"""The decoding workbench on signals with known answers."""
import numpy as np
import pytest

from stargazer.bench import build_zoo, impaired, prepend_random
from stargazer.bench import score
from stargazer.core import BitStream, rng_for
from stargazer.decode import workbench
from stargazer.structure.battery import battery, markov_order_bic

FAST = workbench.Settings(n_surrogates=49)


@pytest.fixture(scope="module")
def zoo():
    return {item.name: item for item in build_zoo(seed=0)}


def test_markov_order_selection():
    rng = np.random.default_rng(0)
    fair = (rng.random(5000) < 0.5).astype(int)
    assert markov_order_bic(fair, 2) == 0
    sticky = np.zeros(5000, dtype=int)
    for i in range(1, 5000):
        sticky[i] = sticky[i - 1] if rng.random() < 0.85 else 1 - sticky[i - 1]
    assert markov_order_bic(sticky, 2) == 1


def test_battery_sees_framing_but_not_memory(zoo):
    rng = np.random.default_rng(1)
    sticky = np.zeros(4000, dtype=int)
    for i in range(1, 4000):
        sticky[i] = sticky[i - 1] if rng.random() < 0.85 else 1 - sticky[i - 1]
    assert battery(sticky, 2, rng_for(1, "a"), n_surrogates=49).p_value > 0.02
    framed_text = np.asarray(zoo["human/ascii_8n1"].payload.bits, dtype=int)
    result = battery(framed_text, 2, rng_for(1, "b"), n_surrogates=49)
    assert result.p_value <= 0.02 and result.null.startswith("markov order")


@pytest.mark.parametrize("name", ["control/fair_bits", "control/markov1_bits", "control/zlib_text"])
def test_controls_stay_at_tier_zero(zoo, name):
    reading = workbench.analyse(zoo[name].payload, name=name, settings=FAST)
    assert reading.tier == 0 and "not proof" in reading.statement
    assert len(reading.tried) >= 9


def test_raster_is_found_and_drawn(zoo):
    item = zoo["message/raster_23x73"]
    reading = workbench.analyse(item.payload, name=item.name, settings=FAST)
    assert reading.tier == 3 and reading.layout["width"] == 23 and reading.layout["rows"] == 73
    assert reading.picture.count("\n") == 72
    shifted = impaired(item, prepend_random, n=17, rng=rng_for(2, "prepend"))
    assert workbench.analyse(shifted.payload, name=shifted.name, settings=FAST).layout["width"] == 23


def test_counting_and_primes_are_anchored(zoo):
    for name, word in (("message/unary_count_1_40", "1, 2, 3"), ("message/unary_primes_300", "prime")):
        reading = workbench.analyse(zoo[name].payload, name=name, settings=FAST)
        assert reading.tier == 4
        assert any(word in e.sentence for e in reading.evidence)


def test_arithmetic_primer_gets_a_tested_reading(zoo):
    reading = workbench.analyse(zoo["message/arithmetic_primer"].payload, name="primer", settings=FAST)
    assert reading.tier == 5
    top = max(reading.evidence, key=lambda e: e.rung)
    assert top.detail["reading"]["operators"] == {252: "+", 253: "*"}
    assert "200 of 200" in top.sentence


def test_shift_register_explains_everything(zoo):
    reading = workbench.analyse(zoo["human/prbs15"].payload, name="prbs", settings=FAST)
    assert reading.tier == 3 and reading.layout == {"kind": "shift register", "length": 15,
                                                    "polynomial": "x^15 + x^14 + 1"}
    assert not any(t["test"].startswith("anchor") for t in reading.tried)


def test_text_is_read_only_when_it_is_text(zoo):
    plain = workbench.analyse(zoo["human/ascii_8n1"].payload, name="8n1", settings=FAST)
    assert plain.tier == 4 and any("It begins" in e.sentence for e in plain.evidence)
    coded = workbench.analyse(zoo["human/hamming74_text"].payload, name="hamming", settings=FAST)
    assert coded.tier == 3


def test_natural_signals_are_not_overclaimed(zoo):
    eclipse = workbench.analyse(zoo["natural/eclipsing_binary"].payload, name="eb", settings=FAST)
    assert eclipse.tier == 1 and eclipse.layout["kind"] == "repeating cycle"
    quasar = workbench.analyse(zoo["natural/drw_quasar"].payload, name="q", settings=FAST)
    assert quasar.tier == 0
    pulsar = workbench.analyse(zoo["natural/pulsar_nulling"].payload, name="psr", settings=FAST)
    assert pulsar.tier == 1 and pulsar.layout["period"] == pytest.approx(0.714, rel=1e-3)


def test_short_streams_are_refused():
    reading = workbench.analyse(BitStream(np.ones(40, dtype=np.uint8)), name="tiny", settings=FAST)
    assert reading.tier == 0 and "too few" in reading.notes[0]


def test_markdown_and_scorecard(zoo):
    reading = workbench.analyse(zoo["message/raster_23x73"].payload, name="raster", settings=FAST)
    text = workbench.to_markdown(reading)
    assert text.startswith("# Reading of raster") and "Everything that was tried" in text and "```" in text
    card = score.run(items=[zoo["control/fair_bits"], zoo["message/unary_count_1_40"]], workers=1, settings=FAST)
    assert card["summary"]["overclaims"] == [] and card["summary"]["message"]["reached"] == 1
    assert "| control/fair_bits |" in score.to_markdown(card)
