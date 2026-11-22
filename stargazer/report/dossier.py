"""Dossiers: one self-contained record per candidate, written from evidence and nothing else.

Every sentence in a dossier is produced from a measured number by a fixed template. The tier a candidate has
reached decides which words may appear: `check_wording` refuses text that uses a word above its tier (for
example "message" before any layout has been recovered) or any phrase on the never-list. A dossier also says
what was not checked and what observation would change its conclusion, because a candidate that cannot be
wrong is not a result.
"""
from __future__ import annotations

import json
import pathlib
import re
from typing import Any

import numpy as np

from ..core import ContentTier, Disposition, RealityTier, content_phrase, reality_phrase
from ..core.tiers import FORBIDDEN_BELOW, NEVER


class WordingError(ValueError):
    """Text claims more than the evidence tier allows."""


def check_wording(text: str, content_tier: ContentTier | int) -> None:
    """Raise if the text uses a word reserved for a higher content tier, or a phrase that is never allowed."""
    lowered = text.lower()
    for phrase in NEVER:
        if phrase in lowered:
            raise WordingError(f'the phrase "{phrase}" is never used in an output')
    for word, tier in FORBIDDEN_BELOW.items():
        if int(content_tier) < int(tier) and re.search(rf"\b{word}\b", lowered):
            raise WordingError(f'the word "{word}" needs content tier {int(tier)}; this output is at '
                               f"{int(content_tier)}")


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def dump_json(obj: Any, path: str | pathlib.Path) -> None:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(obj), indent=1, sort_keys=True) + "\n")


def _probability(p: float) -> str:
    if not np.isfinite(p):
        return "not computed"
    if p >= 0.01:
        return f"{p:.2f}"
    return f"1 in {1 / p:,.0f}" if p > 1e-9 else f"below 1 in a billion ({p:.1e})"


def periodicity_dossier(name: str, search: dict, gates: list[dict], *, promoted: bool, p_final: float,
                        info: dict | None = None, context: dict | None = None) -> dict:
    """The record for one object of a periodicity search.

    `search` is `PeriodicityResult.to_dict()`, `gates` the list of `GateResult.to_dict()`, `p_final` the
    probability after population recalibration, `info` anything known about the object (coordinates,
    magnitudes, catalogue identifiers) and `context` the run it came from (name, configuration hash, numbers
    searched). Reality tier is 1 when promoted, 0 for an unpromoted hit; content tier is 0: a period is not
    content.
    """
    failed = [g["gate"] for g in gates if g["passed"] is False]
    untested = [g["gate"] for g in gates if g["passed"] is None]
    tier = RealityTier.VETTED if promoted else RealityTier.HIT
    decoy = search.get("decoy") or {}
    classical = search.get("classical") or {}
    shape = search.get("harmonic_fraction", float("nan"))
    sentences = [
        f"The best periodic model has a period of {search['period_days']:.1f} days and covers "
        f"{search['cycles_covered']:.1f} cycles in {search['n_points']} nightly points over "
        f"{search['span_days']:.0f} days; its fundamental has an amplitude of {search['amplitude']:.2f} in the "
        "standardised units of the merged light curve (1 is the object's own scatter)"
        + (f", and {100 * shape:.0f}% of its power is in higher harmonics." if np.isfinite(shape) else "."),
        f"Against red noise fitted to this object, with the search over {search['n_frequencies']} trial periods "
        f"included, the probability of a periodic model this good is "
        f"{_probability(decoy.get('p', float('nan')))} ({decoy.get('exceedances', 0)} of "
        f"{decoy.get('draws', 0)} simulated nulls matched it).",
        f"After recalibration against the whole population searched, the probability is {_probability(p_final)}.",
    ]
    for flag in search.get("flags", []):
        sentences.append(f"Note: {flag}.")
    if classical:
        sentences.append("For contrast, a test that assumes white noise would have given "
                         f"{_probability(classical.get('white_noise_fap', float('nan')))}; that number is not "
                         "evidence, because quasar noise is not white.")
    sentences += [g["sentence"] for g in gates]
    if promoted:
        status = ("This object passed every vetting gate that could be evaluated. It is a vetted hit: a periodic "
                  "component that is not explained by red noise as modelled, by the calendar or by the sampling. "
                  "That is not yet a confirmed signal.")
    elif failed:
        status = "This object is not promoted: it failed " + ", ".join(failed) + "."
    else:
        status = "This object is not promoted: its probability does not reach the pre-registered threshold."
    not_checked = ["other surveys of the same object (independent confirmation)",
                   "neighbouring stars on the same detector (local systematics)",
                   "noise models beyond a single damped random walk"] + [f"gate not evaluable: {g}" for g in untested]
    would_change = ("New data are the test. If the period is real, a sinusoid of this period and phase, fitted "
                    "now, will fit data taken later better than red noise alone; if it is a red-noise bump, it "
                    "will not. The prediction is recorded in the prediction ledger before those data exist.")
    record = {"name": name, "kind": "periodicity", "reality_tier": int(tier), "reality": reality_phrase(tier),
              "content_tier": int(ContentTier.NONE), "content": content_phrase(ContentTier.NONE),
              "disposition": Disposition.OPEN.value, "promoted": bool(promoted), "p_final": p_final,
              "sentences": sentences, "status": status, "failed_gates": failed, "not_checked": not_checked,
              "what_would_change_this": would_change, "search": search, "gates": gates, "info": info or {},
              "context": context or {}}
    check_wording(" ".join(sentences + [status, would_change]), ContentTier.NONE)
    return record


def to_markdown(record: dict) -> str:
    """A dossier as Markdown: conclusion first, then evidence, then limits."""
    info = record.get("info", {})
    lines = [f"# {record['name']}", "",
             f"**Status:** {record['reality']} (reality tier {record['reality_tier']} of 4); "
             f"{record['content']} (content tier {record['content_tier']} of 5); disposition: "
             f"{record['disposition']}.", "", record["status"], ""]
    if info:
        lines += ["## Object", ""] + [f"- {key}: {value}" for key, value in info.items()] + [""]
    lines += ["## Evidence", ""] + [f"- {s}" for s in record["sentences"]] + [""]
    lines += ["## Not checked", ""] + [f"- {s}" for s in record["not_checked"]] + [""]
    lines += ["## What would change this", "", record["what_would_change_this"], ""]
    context = record.get("context", {})
    if context:
        lines += ["## Run", ""] + [f"- {key}: {value}" for key, value in context.items()] + [""]
    return "\n".join(lines)


def write(record: dict, directory: str | pathlib.Path) -> pathlib.Path:
    """Write `<name>.json` and `<name>.md` into `directory`; returns the Markdown path."""
    directory = pathlib.Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", record["name"])
    dump_json(record, directory / f"{stem}.json")
    path = directory / f"{stem}.md"
    path.write_text(to_markdown(record))
    return path
