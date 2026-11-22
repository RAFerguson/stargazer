"""The prediction ledger: forecasts committed before the outcome is known.

Adapted from the SeaQuest project (same author), module seaquest/tools/ledger.py at commit 064d667. Changes:
the last hash is found by reading only the end of the file, so appending no longer rescans the whole ledger;
a commitment can be scored only once; and entries carry a `due` field naming the data release or date at
which the forecast can be checked.

A forecast is committed as a salted hash in a hash-chained file. Later it is revealed and scored. The chain
makes any after-the-fact edit detectable, and committing the file to a public Git history gives every entry
an independent timestamp. A pattern is believed when it predicts data it has not seen; this file is the
record of whether it did.
"""
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import secrets
from datetime import datetime, timezone
from typing import Any

GENESIS = "0" * 64


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


class HashChain:
    """An append-only JSON-lines file in which every entry names the hash of the one before it."""

    def __init__(self, path: str | os.PathLike):
        self.path = pathlib.Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.touch(exist_ok=True)

    def entries(self) -> list[dict]:
        with self.path.open() as handle:
            return [json.loads(line) for line in handle if line.strip()]

    def _last_hash(self) -> str:
        with self.path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            if size == 0:
                return GENESIS
            handle.seek(max(0, size - 65536))
            lines = [line for line in handle.read().splitlines() if line.strip()]
        return json.loads(lines[-1])["hash"]

    def append(self, entry: dict[str, Any]) -> dict:
        entry = dict(entry)
        entry["prev"] = self._last_hash()
        entry["utc"] = _now()
        entry["hash"] = _digest(json.dumps(entry, sort_keys=True))
        with self.path.open("a") as handle:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")
        return entry

    def verify(self) -> bool:
        """True if no entry has been altered, removed or reordered."""
        previous = GENESIS
        for entry in self.entries():
            body = {k: v for k, v in entry.items() if k != "hash"}
            if entry.get("prev") != previous or _digest(json.dumps(body, sort_keys=True)) != entry.get("hash"):
                return False
            previous = entry["hash"]
        return True


class PredictionLedger(HashChain):
    """Commit, reveal and score forecasts."""

    def commit(self, prediction: str, *, due: str = "", metadata: dict[str, Any] | None = None) -> dict:
        """Commit a forecast without disclosing it. Keep the returned salt: it is needed to reveal."""
        salt = secrets.token_hex(16)
        entry = self.append({"type": "commit", "commitment": _digest(salt + prediction), "due": due,
                             "metadata": metadata or {}})
        return {"entry_id": entry["hash"], "salt": salt, "utc": entry["utc"]}

    def commit_open(self, prediction: str, *, due: str = "", metadata: dict[str, Any] | None = None) -> dict:
        """Commit a forecast in the clear. Use when there is no reason to hide it, only to date it."""
        entry = self.append({"type": "commit", "commitment": _digest(prediction), "prediction": prediction,
                             "due": due, "metadata": metadata or {}})
        return {"entry_id": entry["hash"], "salt": "", "utc": entry["utc"]}

    def reveal(self, entry_id: str, prediction: str, salt: str) -> dict:
        target = next((e for e in self.entries() if e["hash"] == entry_id and e["type"] == "commit"), None)
        if target is None:
            raise KeyError("unknown commitment")
        valid = _digest(salt + prediction) == target["commitment"]
        return self.append({"type": "reveal", "commit_id": entry_id, "prediction": prediction, "valid": valid})

    def score(self, entry_id: str, outcome: str, correct: bool, *, null_correct: bool | None = None) -> dict:
        """Record the outcome. `null_correct` records whether the null model's forecast would also have passed,
        so that the ledger can show how often a claim beat the null and not merely how often it was right."""
        entries = self.entries()
        if not any(e["hash"] == entry_id and e["type"] == "commit" for e in entries):
            raise KeyError("unknown commitment")
        if any(e["type"] == "score" and e["commit_id"] == entry_id for e in entries):
            raise ValueError("this commitment has already been scored")
        return self.append({"type": "score", "commit_id": entry_id, "outcome": outcome, "correct": bool(correct),
                            "null_correct": null_correct})

    def hit_rate(self, metadata_filter: dict[str, Any] | None = None) -> dict[str, Any]:
        """Scored forecasts, hits, and a Wilson 95% interval on the hit rate."""
        entries = self.entries()
        commits = {e["hash"]: e for e in entries if e["type"] == "commit"}
        disclosed = {h for h, e in commits.items() if "prediction" in e}
        disclosed |= {e["commit_id"] for e in entries if e["type"] == "reveal" and e["valid"]}
        scores = [e for e in entries if e["type"] == "score" and e["commit_id"] in disclosed]
        if metadata_filter:
            scores = [s for s in scores
                      if all(commits[s["commit_id"]]["metadata"].get(k) == v for k, v in metadata_filter.items())]
        n = len(scores)
        hits = sum(1 for s in scores if s["correct"])
        if n == 0:
            return {"n": 0, "hits": 0, "hit_rate": float("nan"), "ci95": [float("nan"), float("nan")]}
        p, z = hits / n, 1.96
        denominator = 1 + z * z / n
        centre = (p + z * z / (2 * n)) / denominator
        half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / denominator
        return {"n": n, "hits": hits, "hit_rate": p, "ci95": [centre - half, centre + half]}
