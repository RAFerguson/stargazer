"""A second blind search on the same quasars: does any light curve sit on discrete levels?

A source whose brightness is keyed between a few fixed levels is the simplest coded signal a slowly sampled
light curve could carry. The test is strict (see `stargazer.decode.symbolise`): nearly every nightly point
must lie within its error bar of one of a few levels that are at least eight error bars apart, and red noise
simulated from the object's own fit must almost never look as discrete. Every object's decoy twin (a
simulation of its fitted red noise at its own sampling) is put through the same test under hidden labels, so
the false-alarm count is measured.

Usage: python run_levels.py register | run SEARCH_ID [engineering]
"""
import json
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from stargazer.core import rng_for
from stargazer.decode import symbolise
from stargazer.detect.periodicity import merge_bands
from stargazer.ledger.trials import TrialsLedger
from stargazer.pipelines import quasar
from stargazer.report import dossier
from stargazer.sources import Cache
from stargazer.trust import blind

HERE = pathlib.Path(__file__).resolve().parent
LEDGER = HERE.parent.parent / "ledger"
SEED = 20261122
SAMPLE = quasar.SampleConfig(healpix_order=2, healpix_pixel=100)
PREPARE = quasar.RunConfig(name="levels-ztf-gaia-pix2-100", seed=SEED)
DEFINITION = {"sample": SAMPLE.to_dict(), "quality_rules": PREPARE.to_dict(), "seed": SEED,
              "test": "decode.symbolise.discreteness_test on the merged nightly light curve, default settings "
                      "(up to 6 levels, at least 8 error bars apart, at most 10% of points off-level, 200 "
                      "red-noise decoys)", "alpha": 0.005,
              "rule": "an object is a hit if the test reports discrete levels at probability <= alpha; the same "
                      "test runs on each object's decoy twin under hidden labels"}


def analyse(subject: blind.Subject) -> dict:
    nightly = merge_bands(subject.bands).binned(1.0)
    result = symbolise.discreteness_test(nightly, rng_for(SEED, "levels", subject.name), alpha=DEFINITION["alpha"])
    return {"label": subject.name, "nights": len(nightly), "k": int(result.k), "detected": bool(result.detected),
            "p_value": float(result.p_value) if np.isfinite(result.p_value) else None,
            "off_level_fraction": float(result.off_level_fraction), "separation": float(result.separation),
            "levels": [float(v) for v in result.levels], "note": result.note}


def main() -> None:
    ledger = TrialsLedger(LEDGER / "trials.jsonl")
    if sys.argv[1] == "register":
        entry = ledger.preregister(PREPARE.name, DEFINITION, description=__doc__)
        print("search id", entry["hash"])
        return
    search_id = sys.argv[2]
    engineering = len(sys.argv) > 3 and sys.argv[3] == "engineering"
    sample = quasar.SampleConfig(healpix_order=4, healpix_pixel=600, max_objects=48) if engineering else SAMPLE
    if not engineering:
        registered = next(e for e in ledger.entries() if e["hash"] == search_id)
        assert registered["config"] == dossier._jsonable(DEFINITION), "configuration differs from the registration"
    cache = Cache(sys.argv[4]) if len(sys.argv) > 4 else None
    started = time.time()
    objects, info, report = quasar.fetch_sample(sample, cache=cache)
    subjects, preparation = quasar.prepare(objects, PREPARE)
    twins = [quasar.decoy_twin(s, rng_for(SEED, "levels-decoy", s.name), PREPARE) for s in subjects]
    sealed = blind.seal(subjects + twins, rng_for(SEED, "levels-seal"))
    out = HERE if not engineering else pathlib.Path(sys.argv[5])
    out.mkdir(parents=True, exist_ok=True)
    (out / "sealed_commitment.txt").write_text(sealed.commitment + "\n")
    with ProcessPoolExecutor(max_workers=4) as pool:
        analysed = list(pool.map(analyse, sealed.subjects, chunksize=4))
    key = sealed.unseal()
    rows = [{**row, **key[row["label"]]} for row in analysed]
    counts = {kind: {"searched": sum(r["kind"] == kind for r in rows),
                     "hits": sum(r["kind"] == kind and r["detected"] for r in rows),
                     "qualified_before_decoys": sum(r["kind"] == kind and r["k"] > 1 for r in rows)}
              for kind in ("real", "decoy")}
    summary = {"sample": report, "preparation": preparation, "counts": counts,
               "hits": [r for r in rows if r["kind"] == "real" and r["detected"]],
               "sealed_commitment": sealed.commitment, "seconds": time.time() - started}
    dossier.dump_json(summary, out / "summary.json")
    dossier.dump_json(sealed.proof(), out / "sealed_proof.json")
    dossier.dump_json(sorted(rows, key=lambda r: (r["kind"], r["name"])), out / "results.json")
    print(json.dumps({"counts": counts, "seconds": round(summary["seconds"])}, indent=1))
    for r in summary["hits"]:
        print("HIT", r["name"], r["k"], r["levels"], r["p_value"], info.get(r["name"]))
    if not engineering:
        ledger.record(search_id, n_objects=counts["real"]["searched"], n_hypotheses=counts["real"]["searched"],
                      summary={"counts": counts, "sealed_commitment": sealed.commitment})


if __name__ == "__main__":
    main()
