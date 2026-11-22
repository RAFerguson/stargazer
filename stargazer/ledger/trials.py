"""The trials ledger: a permanent count of everything that has been tested.

The commonest way for a search like this to fool itself is to forget how many places it has looked. This
ledger is the cure. Before a search touches data it registers what it is about to do and the hash of its full
configuration. Afterwards it records how many hypotheses it tested and the smallest probabilities it met.
Promotions to "discovery" go through `test`, which applies online false-discovery-rate control across the
whole history of the project, not just the current run.

The file is hash-chained (see `ledger.predictions.HashChain`): an entry cannot be altered or dropped without
breaking the chain, and committing the file to a public repository dates every entry independently.
"""
from __future__ import annotations

import os
from typing import Any

from ..core import content_hash
from ..trust.onlinefdr import OnlineFDR
from .predictions import HashChain


class TrialsLedger(HashChain):
    def __init__(self, path: str | os.PathLike, alpha: float = 0.05):
        super().__init__(path)
        self.alpha = alpha

    # -- before the data -----------------------------------------------------------------------------------
    def preregister(self, search: str, config: dict[str, Any], *, description: str = "",
                    planned_hypotheses: int | None = None) -> dict:
        """Declare a search before running it. Returns the entry; keep its `hash` as the search id."""
        return self.append({"type": "preregister", "search": search, "config": config,
                            "config_hash": content_hash(config), "description": description,
                            "planned_hypotheses": planned_hypotheses})

    # -- after the data ------------------------------------------------------------------------------------
    def record(self, search_id: str, *, n_objects: int, n_hypotheses: int, summary: dict[str, Any]) -> dict:
        """Record what a registered search did: how many objects, how many hypotheses in all (objects times
        trial parameters), and a summary of the outcome, including when the outcome is nothing."""
        if not any(e["hash"] == search_id and e["type"] == "preregister" for e in self.entries()):
            raise KeyError("no such registered search; register before running")
        return self.append({"type": "record", "search_id": search_id, "n_objects": int(n_objects),
                            "n_hypotheses": int(n_hypotheses), "summary": summary})

    def test(self, hypothesis: str, p: float, *, search_id: str = "", detail: dict[str, Any] | None = None) -> dict:
        """Put one hypothesis to the project-wide test. The level is fixed by the past, not by `p`."""
        state = self._state()
        reject, level = state.test(float(p))
        return self.append({"type": "test", "hypothesis": hypothesis, "p": float(p), "level": level,
                            "discovery": bool(reject), "search_id": search_id, "detail": detail or {}})

    def note(self, text: str, **fields: Any) -> dict:
        return self.append({"type": "note", "text": text, **fields})

    # -- reading -------------------------------------------------------------------------------------------
    def _state(self) -> OnlineFDR:
        decisions = [bool(e["discovery"]) for e in self.entries() if e["type"] == "test"]
        return OnlineFDR.from_history(decisions, self.alpha)

    def next_level(self) -> float:
        return self._state().next_level()

    def summary(self) -> dict[str, Any]:
        entries = self.entries()
        records = [e for e in entries if e["type"] == "record"]
        tests = [e for e in entries if e["type"] == "test"]
        return {"searches_registered": sum(e["type"] == "preregister" for e in entries),
                "searches_recorded": len(records),
                "objects_examined": sum(e["n_objects"] for e in records),
                "hypotheses_examined": sum(e["n_hypotheses"] for e in records),
                "formal_tests": len(tests),
                "discoveries": sum(bool(e["discovery"]) for e in tests),
                "alpha": self.alpha,
                "next_test_level": self._state().next_level(),
                "chain_intact": self.verify()}
