"""The scorecard: run the decoding workbench on every benchmark signal and compare with the known answer.

Each benchmark item states the highest content tier a correct and honest tool could justify for it. The
workbench's tier is then one of three things:

* **reached**: exactly the justified tier;
* **short**: below it. Something recoverable was missed. This costs sensitivity, not honesty;
* **overclaim**: above it. The tool said more than the signal supports. This is the error that matters.

The controls (random, compressed and encrypted streams, white noise, Poisson events) and the natural items
(red noise, variable stars, an eclipsing binary, a nulling pulsar, a bursting source) are where overclaims
would show. The human and message items are where shortfalls would.
"""
from __future__ import annotations

import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, field

from ..decode import workbench
from .zoo import ZooItem, build_zoo


@dataclass
class ItemScore:
    name: str
    collection: str
    kind: str
    expected_tier: int
    tier: int
    verdict: str                 # "reached", "short" or "overclaim"
    seconds: float
    found: list[str] = field(default_factory=list)       # tests that produced evidence
    layout: dict = field(default_factory=dict)
    statement: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def score_item(item: ZooItem, seed: int = 0, settings: workbench.Settings = workbench.Settings()) -> ItemScore:
    started = time.time()
    reading = workbench.analyse(item.payload, name=item.name, seed=seed, settings=settings)
    expected = int(item.expected["max_content_tier"])
    verdict = "reached" if reading.tier == expected else ("short" if reading.tier < expected else "overclaim")
    return ItemScore(item.name, item.collection, reading.kind, expected, reading.tier, verdict,
                     time.time() - started, [f"{e.test} (tier {e.rung})" for e in reading.evidence],
                     reading.layout, reading.statement)


def _score_packed(args: tuple) -> ItemScore:
    return score_item(*args)


def run(seed: int = 0, size: str = "small", collections: tuple[str, ...] | None = None, workers: int = 4,
        settings: workbench.Settings = workbench.Settings(), items: list[ZooItem] | None = None) -> dict:
    """Score the whole zoo (or the given items). Returns {"items": [...], "summary": {...}}."""
    items = items if items is not None else build_zoo(seed=seed, collections=collections, size=size)
    jobs = [(item, seed, settings) for item in items]
    if workers > 1 and len(jobs) > 2:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            scores = list(pool.map(_score_packed, jobs))
    else:
        scores = [_score_packed(job) for job in jobs]
    summary: dict = {}
    for collection in sorted({s.collection for s in scores}):
        mine = [s for s in scores if s.collection == collection]
        summary[collection] = {verdict: sum(s.verdict == verdict for s in mine)
                               for verdict in ("reached", "short", "overclaim")}
        summary[collection]["items"] = len(mine)
    summary["overclaims"] = [s.name for s in scores if s.verdict == "overclaim"]
    summary["short"] = [f"{s.name} (tier {s.tier} of {s.expected_tier})" for s in scores if s.verdict == "short"]
    summary["seed"] = seed
    summary["size"] = size
    return {"items": [s.to_dict() for s in scores], "summary": summary}


def to_markdown(card: dict) -> str:
    lines = ["| Signal | Kind | Justified tier | Tier given | Verdict | Evidence found |", "|---|---|---|---|---|---|"]
    for s in card["items"]:
        lines.append(f"| {s['name']} | {s['kind']} | {s['expected_tier']} | {s['tier']} | {s['verdict']} | "
                     f"{'; '.join(s['found']) or 'none'} |")
    lines += ["", "| Collection | Items | Reached | Short | Overclaim |", "|---|---|---|---|---|"]
    for collection, counts in card["summary"].items():
        if isinstance(counts, dict):
            lines.append(f"| {collection} | {counts['items']} | {counts['reached']} | {counts['short']} | "
                         f"{counts['overclaim']} |")
    return "\n".join(lines) + "\n"
