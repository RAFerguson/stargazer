"""The benchmark zoo: items with known answers, the registry of their generators, and `build_zoo`.

Every item carries its payload (what a tool sees), its `truth` (everything a scorer needs, including every number
that defines the item) and its `expected` verdict (what a correct, honest tool should conclude). Generators are
registered by item name; each takes a numpy Generator and an optional `size` ("small" or "full") and returns one
`ZooItem`. `build_zoo` gives every generator its own named generator, `rng_for(seed, name)`, so adding or removing an
item never changes any other item.
"""
from __future__ import annotations

import importlib
from dataclasses import dataclass
from typing import Any, Callable, Union

import numpy as np

from stargazer.core import BitStream, ContentTier, EventList, Series, SymbolSequence, rng_for

Payload = Union[Series, EventList, BitStream, SymbolSequence]
COLLECTIONS: tuple[str, ...] = ("natural", "human", "message", "control")
ORIGINS: tuple[str, ...] = ("natural", "human", "designed", "noise")
SIZES: tuple[str, ...] = ("small", "full")
_MODULES = ("natural", "human", "messages", "controls")


@dataclass
class ZooItem:
    """One benchmark signal with its ground truth.

    name: unique and stable, "<collection>/<label>", for example "message/raster_23x73".
    collection: "natural" | "human" | "message" | "control".
    payload: what a tool under test receives.
    truth: everything a scorer needs (parameters, noise-free curves, decoded content, layouts).
    expected: always holds `structured` (bool), `max_content_tier` (int, a ContentTier value) and `origin`
        ("natural" | "human" | "designed" | "noise").
    description: one sentence.
    """

    name: str
    collection: str
    payload: Payload
    truth: dict[str, Any]
    expected: dict[str, Any]
    description: str

    def __post_init__(self) -> None:
        if self.collection not in COLLECTIONS:
            raise ValueError(f"unknown collection {self.collection!r}")
        missing = {"structured", "max_content_tier", "origin"} - set(self.expected)
        if missing:
            raise ValueError(f"expected is missing {sorted(missing)}")
        if self.expected["origin"] not in ORIGINS:
            raise ValueError(f"unknown origin {self.expected['origin']!r}")
        if not 0 <= int(self.expected["max_content_tier"]) <= 5:
            raise ValueError("max_content_tier must lie in 0..5")


Generator = Callable[..., ZooItem]
GENERATORS: dict[str, Generator] = {}


def register(name: str) -> Callable[[Generator], Generator]:
    """Decorator: register a generator under its item name (which must start with its collection)."""
    collection = name.split("/", 1)[0]
    if collection not in COLLECTIONS:
        raise ValueError(f"item name {name!r} does not start with a known collection")

    def wrap(fn: Generator) -> Generator:
        if name in GENERATORS:
            raise ValueError(f"duplicate generator name {name!r}")
        GENERATORS[name] = fn
        fn.item_name = name  # type: ignore[attr-defined]
        return fn

    return wrap


def expect(structured: bool, max_content_tier: ContentTier | int, origin: str) -> dict[str, Any]:
    """The `expected` dict in its canonical form."""
    return {"structured": bool(structured), "max_content_tier": int(max_content_tier), "origin": origin}


def _load_generators() -> None:
    """Import the collection modules so that their generators register themselves."""
    for module in _MODULES:
        importlib.import_module(f"stargazer.bench.{module}")


def build_zoo(seed: int = 0, collections: tuple[str, ...] | None = None, size: str = "small") -> list[ZooItem]:
    """Build every registered item (or those of the given collections), deterministically, in name order.

    "small" keeps payloads modest (bit streams up to about 30,000 bits, series of a few hundred points) so the
    whole zoo builds in a few seconds; "full" uses larger payloads where the item allows it.
    """
    if size not in SIZES:
        raise ValueError(f"size must be one of {SIZES}")
    if collections is not None:
        unknown = set(collections) - set(COLLECTIONS)
        if unknown:
            raise ValueError(f"unknown collections {sorted(unknown)}")
    _load_generators()
    items = []
    for name in sorted(GENERATORS):
        if collections is not None and name.split("/", 1)[0] not in collections:
            continue
        item = GENERATORS[name](rng_for(seed, name), size=size)
        if item.name != name:
            raise RuntimeError(f"generator registered as {name!r} produced {item.name!r}")
        items.append(item)
    return items


# ---------------------------------------------------------------------------------------------------------------
# Shared helpers: built-in English-like text and bit conversions.
# ---------------------------------------------------------------------------------------------------------------

WORDS: tuple[str, ...] = (
    "the", "of", "and", "to", "a", "in", "is", "it", "that", "was", "for", "on", "are", "with", "as", "at", "be",
    "this", "have", "from", "or", "one", "had", "by", "word", "but", "not", "what", "all", "were", "we", "when",
    "your", "can", "said", "there", "use", "an", "each", "which", "she", "do", "how", "their", "if", "will", "up",
    "other", "about", "out", "many", "then", "them", "these", "so", "some", "her", "would", "make", "like", "him",
    "into", "time", "has", "look", "two", "more", "write", "go", "see", "number", "no", "way", "could", "people",
    "my", "than", "first", "water", "been", "call", "who", "oil", "its", "now", "find", "long", "down", "day",
    "did", "get", "come", "made", "may", "part", "star", "light", "sky", "night", "signal", "radio", "distant",
    "world", "small", "large", "between", "under", "never", "always", "again", "together", "answer", "question",
)


def english_text(rng: np.random.Generator, n_chars: int) -> str:
    """English-like ASCII text of exactly n_chars characters: Zipf-weighted words from WORDS in sentences."""
    weights = 1.0 / np.arange(1, len(WORDS) + 1) ** 0.9
    weights /= weights.sum()
    parts: list[str] = []
    length = 0
    while length < n_chars:
        n_words = int(rng.integers(5, 15))
        words = [WORDS[i] for i in rng.choice(len(WORDS), size=n_words, p=weights)]
        if n_words > 8 and rng.random() < 0.4:
            cut = int(rng.integers(3, n_words - 2))
            words[cut] = words[cut] + ","
        sentence = " ".join(words)
        sentence = sentence[0].upper() + sentence[1:] + ". "
        parts.append(sentence)
        length += len(sentence)
    return "".join(parts)[:n_chars]


def bytes_to_bits(data: bytes, msb_first: bool = True) -> np.ndarray:
    """Bits of a byte string, eight per byte, most significant first unless msb_first is False."""
    arr = np.frombuffer(bytes(data), dtype=np.uint8)
    return np.unpackbits(arr, bitorder="big" if msb_first else "little")


def bits_to_bytes(bits: np.ndarray, msb_first: bool = True) -> bytes:
    """Inverse of bytes_to_bits (the bit count must be a multiple of eight)."""
    bits = np.asarray(bits, dtype=np.uint8)
    if bits.size % 8:
        raise ValueError("bit count must be a multiple of 8")
    return np.packbits(bits, bitorder="big" if msb_first else "little").tobytes()


def int_to_bits(value: int, width: int) -> np.ndarray:
    """The width-bit binary form of a non-negative integer, most significant bit first."""
    if value < 0 or value >= 1 << width:
        raise ValueError(f"{value} does not fit in {width} bits")
    return np.array([(value >> (width - 1 - i)) & 1 for i in range(width)], dtype=np.uint8)


def bits_to_int(bits: np.ndarray) -> int:
    """Inverse of int_to_bits."""
    out = 0
    for b in np.asarray(bits, dtype=int):
        out = (out << 1) | int(b)
    return out


def text_bits(text: str) -> np.ndarray:
    """The 8-bit ASCII codes of text, most significant bit first."""
    return bytes_to_bits(text.encode("ascii"))

