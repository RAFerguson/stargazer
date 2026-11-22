"""Dossiers and figures."""
import json

import pytest

from stargazer.report import dossier

SEARCH = {"period_days": 230.0, "cycles_covered": 10.4, "n_points": 400, "span_days": 2400.0, "amplitude": 0.8,
          "harmonic_fraction": 0.05, "n_frequencies": 280, "flags": [],
          "decoy": {"p": 1 / 1001, "draws": 1000, "exceedances": 0}, "classical": {"white_noise_fap": 1e-40}}
GATES = [{"gate": "calendar", "passed": True, "sentence": "Clear of the calendar.", "detail": {}},
         {"gate": "holdout", "passed": False, "sentence": "The early data prefer another period.", "detail": {}},
         {"gate": "bands", "passed": None, "sentence": "No comparison possible.", "detail": {}}]


def test_dossier_states_status_evidence_and_limits(tmp_path):
    record = dossier.periodicity_dossier("gaia:1", SEARCH, GATES, promoted=False, p_final=2e-3,
                                         info={"RA (deg)": 1.0}, context={"run": "test"})
    assert record["reality_tier"] == 0 and record["failed_gates"] == ["holdout"]
    assert "gate not evaluable: bands" in record["not_checked"]
    assert "No prediction was recorded" in record["what_would_change_this"]
    path = dossier.write(record, tmp_path)
    text = path.read_text()
    assert text.startswith("# gaia:1") and "1 in 1,001" in text and "## Not checked" in text
    assert json.loads((tmp_path / "gaia_1.json").read_text())["name"] == "gaia:1"
    promoted = dossier.periodicity_dossier("gaia:2", SEARCH, GATES[:1], promoted=True, p_final=1e-4)
    assert promoted["reality_tier"] == 1 and "recorded in the prediction ledger" in promoted["what_would_change_this"]


def test_figures_render(tmp_path):
    pytest.importorskip("matplotlib")
    import numpy as np

    from stargazer.report import figures

    rng = np.random.default_rng(0)
    rows = [{"kind": kind, "searched": True, "p_raw": float(p)} for kind in ("real", "decoy")
            for p in rng.random(300)]
    cells = [{"label": shape, "amplitude_mag": a, "period_days": p, "trials": 10, "recovered": int(10 * a * 4)}
             for shape in ("sinusoid", "sawtooth", "flare") for a in (0.05, 0.1, 0.2) for p in (150.0, 300.0, 600.0)]
    for theme in ("light", "dark"):
        figures.calibration_figure(rows, tmp_path / f"c-{theme}.png", theme)
        figures.completeness_figure(cells, tmp_path / f"k-{theme}.png", theme)
        assert (tmp_path / f"c-{theme}.png").stat().st_size > 20_000
        assert (tmp_path / f"k-{theme}.png").stat().st_size > 20_000
