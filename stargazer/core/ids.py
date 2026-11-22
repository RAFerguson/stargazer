"""Content hashing, so that every result can name exactly the data and settings it came from."""
from __future__ import annotations

import hashlib
import json
import pathlib
from typing import Any

import numpy as np


def _plain(obj: Any) -> Any:
    """Reduce numpy scalars and arrays to plain Python so JSON output is canonical."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return {"__ndarray__": array_hash(obj), "shape": list(obj.shape), "dtype": str(obj.dtype)}
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float):
        return float(repr(obj)) if np.isfinite(obj) else repr(obj)
    if isinstance(obj, pathlib.Path):
        return str(obj)
    return obj


def content_hash(obj: Any) -> str:
    """SHA-256 of the canonical JSON form of a plain object (dicts sorted, arrays replaced by their own hash)."""
    return hashlib.sha256(json.dumps(_plain(obj), sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def array_hash(array: np.ndarray) -> str:
    """SHA-256 of an array's bytes in a fixed layout and byte order."""
    a = np.ascontiguousarray(array)
    if a.dtype.byteorder == ">":
        a = a.byteswap().view(a.dtype.newbyteorder("<"))
    return hashlib.sha256(a.tobytes()).hexdigest()


def file_hash(path: str | pathlib.Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()
