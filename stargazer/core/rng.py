"""Seeded randomness. Every stochastic step draws from a generator named after what it is for, so results do
not change when an unrelated step is added or reordered."""
from __future__ import annotations

import hashlib

import numpy as np


def rng_for(seed: int, *names: object) -> np.random.Generator:
    """A generator determined by the run seed and a path of names, for example rng_for(7, "decoys", source_id)."""
    digest = hashlib.sha256("/".join(str(n) for n in names).encode()).digest()
    words = [int.from_bytes(digest[i:i + 4], "little") for i in range(0, 16, 4)]
    return np.random.default_rng(np.random.SeedSequence([int(seed) & 0xFFFFFFFF, *words]))
