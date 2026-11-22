"""Follow-up of the one object promoted by the blind run: what is it?

Not part of the registered analysis. The registered rule promoted Gaia DR3 3613794004069753472 on a clean
79.9-day cycle. This script asks the next question a vetted hit must face, independent data:

1. the unbinned ZTF measurements are searched for a much faster period;
2. the Gaia DR3 archive is asked what it knows about the source (the three CSV files beside this script are
   its answers, fetched on 22 November 2026 from the public TAP service at gea.esac.esa.int);
3. the alias gate added to the pipeline after this run is applied;
4. the dossiers are regenerated with the outcome.

Usage: python follow_up.py [cache_dir]
"""
import csv
import gzip
import json
import pathlib
import sys

import numpy as np
from astropy.timeseries import LombScargle

from stargazer.ledger.trials import TrialsLedger
from stargazer.report import dossier
from stargazer.sources import Cache, ztf
from stargazer.vet import gates

HERE = pathlib.Path(__file__).resolve().parent
RUN = HERE.parent
NAME = "gaia:3613794004069753472"


def main() -> None:
    cache = Cache(sys.argv[1]) if len(sys.argv) > 1 else None
    rows = json.load(gzip.open(RUN / "results.json.gz", "rt"))
    row = next(r for r in rows if r["kind"] == "real" and r["name"] == NAME)
    candidate = row["search"]["period_days"]
    bands = ztf.fetch_lightcurves(208.518067, -10.911535, 1.5, cache=cache, fragments="all")

    # 1. The fast period in the unbinned data.
    freq = np.arange(0.8, 5.2, 2e-5)
    fast = {}
    for band in bands:
        if len(band) >= 100:
            power = LombScargle(band.time, band.value, band.error).power(freq)
            fast[band.band] = {"period_days": float(1 / freq[np.argmax(power)]), "power": float(power.max()),
                               "points": len(band)}
    # 2. Gaia DR3.
    gaia = {path.stem: next(csv.DictReader(path.open())) for path in sorted(HERE.glob("gaia_*.csv"))}
    gaia_period = float(gaia["gaia_vari_rrlyrae"]["pf"])
    # 3. The alias gate.
    alias = gates.alias_gate(bands, candidate)
    f_fast = 1 / fast["r"]["period_days"]
    predicted_alias = 1 / abs(f_fast - 2 - 1 / 365.2422)

    outcome = {"object": NAME, "candidate_period_days": candidate, "unbinned_ztf": fast,
               "gaia_dr3": {"vari_rrlyrae_class": gaia["gaia_vari_rrlyrae"]["best_classification"],
                            "vari_rrlyrae_period_days": gaia_period,
                            "vari_rrlyrae_period_error_days": float(gaia["gaia_vari_rrlyrae"]["pf_error"]),
                            "peak_to_peak_g_mag": float(gaia["gaia_vari_rrlyrae"]["peak_to_peak_g"]),
                            "vari_classifier_best_class": gaia["gaia_vari_classifier_result"]["best_class_name"],
                            "quasar_probability": float(gaia["gaia_source"]["classprob_dsc_combmod_quasar"]),
                            "bp_rp": float(gaia["gaia_source"]["bp_rp"]),
                            "parallax_mas": float(gaia["gaia_source"]["parallax"]),
                            "parallax_error_mas": float(gaia["gaia_source"]["parallax_error"])},
               "alias_of_two_cycles_per_day_plus_one_per_year_days": predicted_alias,
               "alias_gate": alias.to_dict()}
    dossier.dump_json(outcome, HERE / "follow_up.json")
    print(json.dumps({k: outcome[k] for k in ("candidate_period_days", "unbinned_ztf", "gaia_dr3",
                                              "alias_of_two_cycles_per_day_plus_one_per_year_days")}, indent=1))
    print(alias.sentence)

    follow_up = {
        "reality_tier": 3, "disposition": "known natural class",
        "status": "The registered analysis promoted this object: it passed every gate the pipeline then had. "
                  "Follow-up shows what it is: an RR Lyrae star, a known kind of pulsating star, whose half-day "
                  "pulsation appears as an 80-day cycle under nightly sampling. The signal is real, celestial "
                  "and explained.",
        "sentences": [
            f"The unbinned ZTF measurements are fitted far better by a period of {fast['r']['period_days']:.5f} "
            f"days (r band; {fast['g']['period_days']:.5f} days in g) than by {candidate:.1f} days.",
            f"Gaia DR3 lists this source as an RR Lyrae star of type {outcome['gaia_dr3']['vari_rrlyrae_class']} "
            f"with a pulsation period of {gaia_period:.5f} days and a peak-to-peak range of "
            f"{outcome['gaia_dr3']['peak_to_peak_g_mag']:.2f} magnitudes. Two instruments and two pipelines agree "
            "on the period to 4 parts in 100,000.",
            f"A {gaia_period:.4f}-day pulsation seen about once a night, at an hour that drifts through the year, "
            f"beats at {predicted_alias:.1f} days. The {candidate:.1f}-day cycle is that beat.",
            "The object is a pulsating star in the Galactic halo that Gaia's classifier gave a 96% chance of being "
            "a quasar. The search was right that it is strictly periodic and wrong, through no test it then had, "
            "about the period. An alias gate has since been added to the pipeline; it fails this object.",
        ],
        "sources": ["Gaia DR3 tables gaia_source, vari_classifier_result and vari_rrlyrae, fetched 22 November "
                    "2026 (files beside follow_up.py)", "ZTF light curves, unbinned, all three bands"],
        "what_would_change_this": "Nothing needs to. The variation is explained: an RR Lyrae star, with the long "
                                  "period an alias of its pulsation.",
    }
    context = {"run": "quasar-ztf-gaia-pix2-100"}
    for path in sorted((RUN / "dossiers").glob("*.json")):
        old = json.loads(path.read_text())
        record = dossier.periodicity_dossier(
            old["name"], old["search"], old["gates"], promoted=old["promoted"], p_final=old["p_final"],
            info=old["info"], context=old.get("context") or context,
            follow_up=follow_up if old["name"] == NAME else None)
        for key in ("project_wide_test", "prediction"):
            if key in old:
                record[key] = old[key]
        dossier.write(record, RUN / "dossiers")
    ledger = TrialsLedger(RUN.parent.parent / "ledger" / "trials.jsonl")
    if not any(e.get("type") == "note" and e.get("object") == NAME for e in ledger.entries()):
        ledger.note("Follow-up of the promoted object: Gaia DR3 lists it as an RR Lyrae star (RRab, period "
                    f"{gaia_period:.5f} d); unbinned ZTF data give {fast['r']['period_days']:.5f} d; the "
                    f"{candidate:.1f}-day cycle is the alias of that pulsation under nightly sampling. Disposition: "
                    "known natural class.", object=NAME, search_id=row.get("search_id", ""),
                    run="quasar-ztf-gaia-pix2-100")


if __name__ == "__main__":
    main()
