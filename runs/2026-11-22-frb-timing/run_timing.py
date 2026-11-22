"""Do the bursts of two repeating fast radio bursts keep time?

Two questions per source, each answered against surrogates that keep the burst rate and its slow changes and
destroy only fine timing (see `stargazer.detect.intervals`):

* clock: do the arrival times fall on the ticks of a regular clock, with some ticks empty?
* spacing: are the gaps between bursts whole multiples of one unit?

Sources, chosen before looking at their timing with this code:

* FRB 20220912A, 696 bursts in 68 Nançay scans, moved to the Solar System barycentre, with the scan log as
  exposure;
* FRB 20220529, 643 FAST bursts, already barycentric, without exposure information.

Four tests in all. A test is a hit if its probability is at or below 0.01 / 4 (Bonferroni over the four),
with 999 surrogates each. Period range: the default of the tests, a tenth of the median gap to twice it.
FRB 20180916B is left out: its compiled observing windows are incomplete, and a timing test without honest
exposure would not mean much.

Usage: python run_timing.py register | run SEARCH_ID | engineering
"""
import json
import pathlib
import sys

import numpy as np

from stargazer.core import EventList, rng_for
from stargazer.detect import intervals
from stargazer.ledger.trials import TrialsLedger
from stargazer.report import dossier
from stargazer.sources import frb

HERE = pathlib.Path(__file__).resolve().parent
LEDGER = HERE.parent.parent / "ledger"
SEED = 20261122
DEFINITION = {"sources": [{"name": "FRB 20220912A", "barycentre": True}, {"name": "FRB 20220529", "barycentre": False}],
              "tests": ["intervals.clock_test", "intervals.spacing_test"], "n_surrogates": 999,
              "alpha": 0.01, "hit_rule": "probability <= alpha / 4", "seed": SEED,
              "period_range": "default_period_range: 0.1 to 2 median gaps"}


def probability(value) -> float:
    return float(getattr(value, "p", value))


def analyse(name: str, events: EventList) -> dict:
    clock = intervals.clock_test(events, rng_for(SEED, "frb-clock", name), n_surrogates=DEFINITION["n_surrogates"],
                                 alpha=DEFINITION["alpha"] / 4)
    spacing = intervals.spacing_test(events, rng_for(SEED, "frb-spacing", name),
                                     n_surrogates=DEFINITION["n_surrogates"], alpha=DEFINITION["alpha"] / 4)
    gaps = np.diff(events.time)
    return {"bursts": len(events), "span_days": float(events.span / 86400.0) if hasattr(events, "span")
            else float((events.time[-1] - events.time[0]) / 86400.0), "median_gap_s": float(np.median(gaps)),
            "clock": {"p": probability(clock.p_value), "detected": bool(clock.detected), "period_s": float(clock.period),
                      "searched_s": [float(clock.min_period), float(clock.max_period)],
                      "occupancy": float(clock.occupancy)},
            "spacing": {"p": probability(spacing.p_value), "detected": bool(spacing.detected),
                        "quantum_s": float(spacing.quantum),
                        "searched_s": [float(spacing.min_quantum), float(spacing.max_quantum)]}}


def synthetic_like(events: EventList, rng: np.random.Generator) -> EventList:
    """Poisson bursts with the same count, spread over the same exposure (or the same span): used only to check
    that the code runs, so that the real answer is not seen before registration."""
    if events.exposure:
        windows = np.asarray(events.exposure, dtype=float)
        lengths = windows[:, 1] - windows[:, 0]
        which = rng.choice(len(windows), size=len(events), p=lengths / lengths.sum())
        time = np.sort(windows[which, 0] + rng.random(len(events)) * lengths[which])
    else:
        time = np.sort(rng.uniform(events.time[0], events.time[-1], len(events)))
    return EventList(time, exposure=events.exposure, source_id="synthetic", time_unit="s")


def main() -> None:
    ledger = TrialsLedger(LEDGER / "trials.jsonl")
    mode = sys.argv[1]
    if mode == "register":
        entry = ledger.preregister("frb-timing-2026-11-22", DEFINITION, description=__doc__)
        print("search id", entry["hash"])
        return
    loaded = {s["name"]: frb.fetch_bursts(s["name"], barycentre=s["barycentre"]) for s in DEFINITION["sources"]}
    if mode == "engineering":
        for name, events in loaded.items():
            fake = synthetic_like(events, rng_for(1, "engineering", name))
            print(name, "(synthetic stand-in)", json.dumps(analyse(name + " synthetic", fake), default=float))
        return
    search_id = sys.argv[2]
    registered = next(e for e in ledger.entries() if e["hash"] == search_id)
    assert registered["config"] == dossier._jsonable(DEFINITION), "configuration differs from the registration"
    summary = {name: analyse(name, events) for name, events in loaded.items()}
    for name, events in loaded.items():
        summary[name]["time_standard"] = events.meta.get("time_standard", "")
        summary[name]["citation"] = events.provenance.citation if events.provenance else ""
    dossier.dump_json(summary, HERE / "summary.json")
    print(json.dumps(summary, indent=1, default=float))
    hits = sum(summary[n][t]["detected"] for n in summary for t in ("clock", "spacing"))
    ledger.record(search_id, n_objects=len(summary), n_hypotheses=4,
                  summary={"hits": int(hits), "probabilities": {n: {t: summary[n][t]["p"] for t in ("clock", "spacing")}
                                                                for n in summary}})


if __name__ == "__main__":
    main()
