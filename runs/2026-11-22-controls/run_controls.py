"""Control objects, searched one at a time with the registered periodicity configuration and the same gates.

Chosen before looking at their results with this code:

* PG 1553+113 (4FGL J1555.7+1111): a blazar with a published claim of a roughly 2.2-year cycle in Fermi-LAT
  gamma rays. A positive control of the "famous candidate" kind: the question is not whether the claim is
  right but what this pipeline says when given the same public data.
* KIC 8462852 (Boyajian's star): deep, irregular dips in Kepler photometry. A strange object that is not known
  to be periodic: the pipeline should not turn strangeness into a period.

Usage: python run_controls.py register | run SEARCH_ID
"""
import json
import pathlib
import sys

import numpy as np

from stargazer.core import Series, rng_for
from stargazer.detect import periodicity
from stargazer.ledger.trials import TrialsLedger
from stargazer.report import dossier
from stargazer.sources import fermi, kepler
from stargazer.vet import gates

HERE = pathlib.Path(__file__).resolve().parent
LEDGER = HERE.parent.parent / "ledger"
CONFIG = periodicity.SearchConfig(max_decoys=5000)
DEFINITION = {"search": CONFIG.__dict__, "seed": 20261122, "hit_p": 0.01, "promote_p": 1e-3,
              "objects": ["4FGL J1555.7+1111 (Fermi LCR weekly, as magnitudes)",
                          "KIC 8462852 (Kepler long cadence, PDCSAP, relative flux, nightly means)"]}


def as_magnitudes(series: Series) -> Series:
    """Flux to magnitudes relative to the median, so that multiplicative variability becomes additive."""
    scale = np.median(series.value)
    return Series(series.time, -2.5 * np.log10(series.value / scale), 1.0857 * series.error / series.value,
                  series.band, series.source_id, "mag", series.time_standard, dict(series.meta), series.provenance)


def load() -> dict[str, Series]:
    blazar = fermi.fetch_lightcurve("4FGL J1555.7+1111", cadence="weekly")
    blazar = blazar.select(blazar.value > 0)
    star = kepler.fetch_lightcurve(8462852)
    return {"PG 1553+113 (Fermi-LAT weekly)": as_magnitudes(blazar), "KIC 8462852 (Kepler)": star}


def main() -> None:
    ledger = TrialsLedger(LEDGER / "trials.jsonl")
    if sys.argv[1] == "register":
        entry = ledger.preregister("controls-2026-11-22", DEFINITION, description=__doc__)
        print("search id", entry["hash"])
        return
    search_id = sys.argv[2]
    registered = next(e for e in ledger.entries() if e["hash"] == search_id)
    assert registered["config"] == dossier._jsonable(DEFINITION), "configuration differs from the registration"
    summary, n_frequencies = {}, 0
    for name, series in load().items():
        result = periodicity.search(series, rng_for(DEFINITION["seed"], "controls", name), CONFIG)
        checks = [g.to_dict() for g in gates.run_gates(series, [series], result, CONFIG)]
        passed, failed = gates.verdict([gates.GateResult(**g) for g in checks])
        promoted = bool(result.p <= DEFINITION["promote_p"] and passed)
        record = dossier.periodicity_dossier(name, result.to_dict(), checks, promoted=promoted, p_final=result.p,
                                             info={"points": len(series), "span (days)": round(series.span)},
                                             context={"run": "controls-2026-11-22", "search id": search_id,
                                                      "population recalibration": "none (single objects)"})
        if promoted:
            formal = ledger.test(f"controls: {name} is periodic at {result.period_days:.1f} d", result.p,
                                 search_id=search_id)
            record["project_wide_test"] = {k: formal[k] for k in ("p", "level", "discovery")}
        dossier.write(record, HERE / "dossiers")
        n_frequencies += result.n_frequencies
        summary[name] = {"period_days": result.period_days, "p": result.p, "statistic": result.statistic,
                         "decoys": result.decoy.draws, "exceedances": result.decoy.exceedances,
                         "coherence_p": result.coherence_p, "white_noise_fap": result.classical["white_noise_fap"],
                         "failed_gates": failed, "promoted": promoted, "flags": result.flags,
                         "noise": result.noise, "nights": result.n_points}
        print(name, json.dumps(summary[name], default=float))
    dossier.dump_json(summary, HERE / "summary.json")
    ledger.record(search_id, n_objects=len(summary), n_hypotheses=n_frequencies, summary=summary)


if __name__ == "__main__":
    main()
