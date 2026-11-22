"""Blind analysis: the search does not know which of its inputs are real.

Three kinds of subject go into a sealed set under meaningless labels, in random order:

* **real** objects, exactly as observed;
* **decoy twins**: for each real object, a simulation of its fitted null at the same times with the same error
  bars. Nothing in a decoy is real, so anything the pipeline promotes from a decoy is a false alarm;
* **plants**: real objects with a known artificial signal added. The fraction of plants the pipeline promotes
  is its sensitivity, measured in the presence of the real noise.

The whole pipeline, vetting gates included, runs on the sealed set. Only when the list of promoted labels is
final is the set unsealed. Then

* promoted decoys estimate the number of false alarms among the promoted real series (the target–decoy
  estimate used in mass spectrometry: about one false real promotion for each promoted decoy),
* promoted plants give the completeness as a function of the planted signal, and hence upper limits,
* and no threshold or gate can have been tuned, knowingly or not, on the real answer.

A salted hash of the label table is recorded before the analysis, so the unsealing can be checked afterwards.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

from ..core import Series


@dataclass
class Subject:
    """One thing to analyse: an object's bands, and (hidden once sealed) what it really is."""

    name: str                                  # the real object's name; an opaque label once sealed
    bands: list[Series]
    kind: str = "real"                         # "real", "decoy" or "plant"; "sealed" once sealed
    parameters: dict = field(default_factory=dict)   # for plants: what was planted


@dataclass
class SealedSet:
    subjects: list[Subject]              # what the pipeline sees: opaque names, kind "sealed", no parameters
    commitment: str                      # sha256 of the salted key; publish before analysing
    _key: dict[str, dict] = field(default_factory=dict, repr=False)
    _salt: str = field(default="", repr=False)

    def unseal(self) -> dict[str, dict]:
        """The label table: label -> {"kind", "name", "parameters"}. Call only when the analysis is final."""
        return self._key

    def proof(self) -> dict:
        """What a third party needs to check that the table was fixed before the analysis."""
        return {"salt": self._salt, "key": self._key, "commitment": self.commitment}

    @staticmethod
    def verify(proof: dict) -> bool:
        text = proof["salt"] + json.dumps(proof["key"], sort_keys=True)
        return hashlib.sha256(text.encode()).hexdigest() == proof["commitment"]


def seal(subjects: list[Subject], rng: np.random.Generator) -> SealedSet:
    """Hide what each subject is: random order, random labels, nothing in the data that tells the kinds apart.

    Band names are kept (the pipeline needs them); survey object identifiers, provenance and all other metadata
    are stripped from the sealed copies, because they would give the kind away.
    """
    order = rng.permutation(len(subjects))
    key: dict[str, dict] = {}
    sealed: list[Subject] = []
    for index in order:
        subject = subjects[index]
        label = bytes(rng.integers(0, 256, 8, dtype=np.uint8)).hex()
        key[label] = {"kind": subject.kind, "name": subject.name, "parameters": dict(subject.parameters)}
        bands = [Series(b.time, b.value, b.error, b.band, label, b.kind, b.time_standard) for b in subject.bands]
        sealed.append(Subject(label, bands, "sealed"))
    salt = bytes(rng.integers(0, 256, 16, dtype=np.uint8)).hex()
    commitment = hashlib.sha256((salt + json.dumps(key, sort_keys=True)).encode()).hexdigest()
    return SealedSet(sealed, commitment, key, salt)


def decoy_estimate(n_real_promoted: int, n_decoy_promoted: int, n_real: int, n_decoy: int) -> dict:
    """What promoted decoys say about promoted real series.

    Returns the estimated number of false alarms among the real promotions, the estimated false-discovery
    fraction (with the customary +1 so that it is never optimistically zero), and a 95% upper limit on the
    pipeline's per-object false-alarm probability from the decoys alone.
    """
    scale = n_real / n_decoy if n_decoy else float("nan")
    expected_false = n_decoy_promoted * scale
    fdr = min(1.0, (n_decoy_promoted + 1) * scale / max(n_real_promoted, 1)) if n_decoy else float("nan")
    from scipy.stats import beta

    upper = float(beta.ppf(0.95, n_decoy_promoted + 1, n_decoy - n_decoy_promoted)) if n_decoy else float("nan")
    return {"real_promoted": n_real_promoted, "decoy_promoted": n_decoy_promoted, "n_real": n_real,
            "n_decoy": n_decoy, "expected_false_among_real": expected_false, "estimated_fdr": fdr,
            "false_alarm_probability_upper95": upper}
