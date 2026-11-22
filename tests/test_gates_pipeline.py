"""Vetting gates, the dossier, and the blind pipeline end to end on simulated objects."""
import json

import numpy as np
import pytest

from stargazer.bench.natural import simulate_drw, survey_times
from stargazer.core import Series, rng_for
from stargazer.detect import periodicity as per
from stargazer.pipelines import quasar
from stargazer.report import dossier
from stargazer.trust import inject
from stargazer.vet import gates, quality

FAST = per.SearchConfig(max_decoys=150, coherence_templates=60)


def _bands(seed, amplitude=0.0, period=230.0):
    rng = rng_for(seed, "pipeline-object")
    tg, tr = survey_times(rng, nights=(25, 40)), survey_times(rng, nights=(30, 45))
    both = np.sort(np.concatenate([tg, tr]))
    signal = simulate_drw(both, 0.12, 200.0, rng) + amplitude * np.sin(2 * np.pi * both / period)
    bands = []
    for name, time, error, scale in (("g", tg, 0.04, 1.0), ("r", tr, 0.03, 0.9)):
        value = 18.5 + scale * signal[np.searchsorted(both, time)] + error * rng.standard_normal(time.size)
        bands.append(Series(time, value, np.full(time.size, error), name, f"sim:{seed}:{name}"))
    return bands


def test_calendar_gate():
    assert not gates.calendar_gate(365.0, 2400.0).passed
    assert not gates.calendar_gate(182.0, 2400.0).passed
    assert gates.calendar_gate(230.0, 2400.0).passed
    assert not gates.calendar_gate(29.5, 2400.0).passed
    from stargazer.detect.calendar import calendar_periods
    assert abs(calendar_periods()["synodic month"] - 29.53) < 0.01


def test_gates_pass_a_strong_real_signal():
    bands = _bands(3, amplitude=0.2)
    merged = per.merge_bands(bands)
    result = per.search(merged, rng_for(1, "gates"), FAST)
    assert result.period_days == pytest.approx(230.0, rel=0.05)
    results = {g.gate: g for g in gates.run_gates(merged, bands, result, FAST)}
    assert set(results) == {"calendar", "window", "seasons", "bands", "holdout", "coherence"}
    for name in ("calendar", "window", "seasons", "bands", "holdout"):
        assert results[name].passed, results[name].sentence
    assert gates.verdict(list(results.values()))[0] in (True, False)
    assert all(isinstance(g.sentence, str) and g.sentence for g in results.values())


def test_quality_rules():
    rng = np.random.default_rng(0)
    time = np.sort(np.concatenate([np.arange(100.0), np.arange(100.0) + 0.1]))
    series = Series(time, 0.08 * rng.standard_normal(200), np.full(200, 0.04))
    scale = quality.error_scale(quality.pair_residuals(series))
    assert scale == pytest.approx(2.0, rel=0.25)
    assert quality.rescale_errors(series, 0.5).error[0] == 0.04          # errors are never shrunk
    spiked = series.with_values(np.where(np.arange(200) == 7, 5.0, series.value))
    assert quality.clip_outliers(spiked).meta["n_clipped"] == 1


def test_dossier_wording_rules():
    dossier.check_wording("A periodic component with a period of 200 days.", 0)
    with pytest.raises(dossier.WordingError):
        dossier.check_wording("This message has been decoded.", 0)
    with pytest.raises(dossier.WordingError):
        dossier.check_wording("An extraterrestrial intelligence did this.", 5)
    dossier.check_wording("The message layout is a raster.", 3)


def test_blind_pipeline_end_to_end(tmp_path):
    objects = {f"sim:{i:02d}": _bands(40 + i, amplitude=0.25 if i < 2 else 0.0) for i in range(8)}
    config = quasar.RunConfig(name="test-run", seed=4, workers=1, plant_fraction=0.5, deep_decoys=400,
                              search=per.SearchConfig(max_decoys=200, coherence_templates=60))
    summary = quasar.run(objects, config, tmp_path / "run", ledger_dir=tmp_path / "ledger")
    assert summary["searched"] == {"real": 8, "decoy": 8, "plant": 4}
    assert summary["promoted"]["decoy"] == 0
    names = {row["name"] for row in summary["top_real"][:2]}
    assert names == {"sim:00", "sim:01"}
    import gzip
    rows = json.loads(gzip.open(tmp_path / "run" / "results.json.gz", "rt").read())
    assert len(rows) == 20 and {r["kind"] for r in rows} == {"real", "decoy", "plant"}
    table = (tmp_path / "run" / "results.csv").read_text().splitlines()
    assert len(table) == 21 and table[0].startswith("label,kind,name")
    proof = json.loads((tmp_path / "run" / "sealed_proof.json").read_text())
    from stargazer.trust.blind import SealedSet
    assert SealedSet.verify(proof)
    assert proof["commitment"] == (tmp_path / "run" / "sealed_commitment.txt").read_text().strip()
    from stargazer.ledger.trials import TrialsLedger
    ledger = TrialsLedger(tmp_path / "ledger" / "trials.jsonl").summary()
    assert ledger["chain_intact"] and ledger["searches_recorded"] == 1 and ledger["objects_examined"] == 8
    assert summary["dossiers"] and all(open(path).read().startswith("# sim:") for path in summary["dossiers"])


def test_a_registered_run_must_match_its_registration(tmp_path):
    config = quasar.RunConfig(name="registered", seed=1, workers=1)
    entry = quasar.register(config, tmp_path / "ledger", sample=quasar.SampleConfig(max_objects=5),
                            description="test", code_version="abc")
    assert entry["config"]["sample"]["max_objects"] == 5 and entry["config"]["code_version"] == "abc"
    changed = quasar.RunConfig(name="registered", seed=1, workers=1, promote_p=0.5)
    with pytest.raises(ValueError):
        quasar.run({}, changed, tmp_path / "run", ledger_dir=tmp_path / "ledger", search_id=entry["hash"])
    with pytest.raises(KeyError):
        quasar.run({}, config, tmp_path / "run", ledger_dir=tmp_path / "ledger", search_id="0" * 64)
