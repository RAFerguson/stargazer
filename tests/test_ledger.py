"""The ledgers: hash chain, predictions, trials."""
import json

import pytest

from stargazer.ledger.predictions import HashChain, PredictionLedger
from stargazer.ledger.trials import TrialsLedger


def test_hash_chain_detects_tampering(tmp_path):
    chain = HashChain(tmp_path / "chain.jsonl")
    for i in range(5):
        chain.append({"type": "note", "i": i})
    assert chain.verify() and len(chain.entries()) == 5
    lines = (tmp_path / "chain.jsonl").read_text().splitlines()
    entry = json.loads(lines[2])
    entry["i"] = 99
    lines[2] = json.dumps(entry, sort_keys=True)
    (tmp_path / "chain.jsonl").write_text("\n".join(lines) + "\n")
    assert not chain.verify()
    (tmp_path / "chain.jsonl").write_text("\n".join(lines[:2] + lines[3:]) + "\n")
    assert not chain.verify()


def test_predictions_commit_reveal_score(tmp_path):
    ledger = PredictionLedger(tmp_path / "predictions.jsonl")
    sealed = ledger.commit("object A will brighten", due="next release", metadata={"run": "x"})
    opened = ledger.commit_open("object B keeps its period", metadata={"run": "x"})
    assert not ledger.reveal(sealed["entry_id"], "something else", sealed["salt"])["valid"]
    assert ledger.reveal(sealed["entry_id"], "object A will brighten", sealed["salt"])["valid"]
    ledger.score(sealed["entry_id"], "it brightened", True, null_correct=False)
    ledger.score(opened["entry_id"], "it did not", False)
    with pytest.raises(ValueError):
        ledger.score(opened["entry_id"], "again", True)
    with pytest.raises(KeyError):
        ledger.score("0" * 64, "nothing", True)
    rate = ledger.hit_rate({"run": "x"})
    assert rate["n"] == 2 and rate["hits"] == 1 and rate["ci95"][0] < 0.5 < rate["ci95"][1]
    assert ledger.verify()


def test_trials_ledger_counts_everything(tmp_path):
    ledger = TrialsLedger(tmp_path / "trials.jsonl", alpha=0.05)
    with pytest.raises(KeyError):
        ledger.record("0" * 64, n_objects=1, n_hypotheses=1, summary={})
    registered = ledger.preregister("test-search", {"threshold": 1e-3}, planned_hypotheses=100)
    first_level = ledger.next_level()
    miss = ledger.test("object 1 is periodic", 0.2, search_id=registered["hash"])
    hit = ledger.test("object 2 is periodic", 1e-9, search_id=registered["hash"])
    assert miss["level"] == pytest.approx(first_level) and not miss["discovery"] and hit["discovery"]
    ledger.record(registered["hash"], n_objects=100, n_hypotheses=58_000, summary={"promoted": 1})
    summary = ledger.summary()
    assert summary["hypotheses_examined"] == 58_000 and summary["discoveries"] == 1 and summary["formal_tests"] == 2
    assert summary["chain_intact"] and summary["searches_registered"] == 1
    assert TrialsLedger(tmp_path / "trials.jsonl").next_level() == pytest.approx(ledger.next_level())
