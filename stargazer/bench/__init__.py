"""Benchmark collections with known answers: natural, human-made, message and control signals, their generators,
and channel impairments. Importing this package registers every generator in GENERATORS."""

from . import controls, human, messages, natural  # noqa: F401  (imported to register their generators)
from .channels import (add_noise, drop_bits, drop_events, drop_points, erase_bits, flip_bits, impaired, jitter_times,
                       prepend_random, truncate)
from .zoo import COLLECTIONS, GENERATORS, ORIGINS, ZooItem, build_zoo, english_text, expect, register

__all__ = [
    "COLLECTIONS", "GENERATORS", "ORIGINS", "ZooItem", "add_noise", "build_zoo", "controls", "drop_bits",
    "drop_events", "drop_points", "english_text", "erase_bits", "expect", "flip_bits", "human", "impaired",
    "jitter_times", "messages", "natural", "prepend_random", "register", "truncate",
]
