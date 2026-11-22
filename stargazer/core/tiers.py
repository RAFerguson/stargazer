"""The two ladders every output is placed on, and the wording each rung permits.

The reality ladder says how much scrutiny a signal's existence and origin have survived. The content ladder
says how much of what it carries has been recovered. They are independent: a pulsar is explained (low on the
first) and can still be analysed on the second.
"""
from __future__ import annotations

from enum import Enum, IntEnum


class RealityTier(IntEnum):
    HIT = 0          # a statistic crossed its threshold in one dataset
    VETTED = 1       # not an artefact of the instrument, the reduction or the sampling
    CONFIRMED = 2    # recurs in independent data, or a committed prediction came true
    CELESTIAL = 3    # tied to a position on the sky, not to the Earth, a satellite or the telescope
    UNEXPLAINED = 4  # the natural models tried do not account for it


class ContentTier(IntEnum):
    NONE = 0         # content not analysed, or nothing found
    SYMBOLS = 1      # discrete states recovered, stable across methods
    STRUCTURE = 2    # sequence structure beyond what the fitted natural model produces
    LAYOUT = 3       # framing, dimensions or coding recovered, with compression gain and held-out prediction
    ANCHORED = 4     # part of the content matches shared knowledge, corrected for everything tried
    HYPOTHESIS = 5   # a proposed reading, labelled as a hypothesis


class Disposition(str, Enum):
    """What a candidate turned out to be. Everything that fails a rung still gets one."""

    OPEN = "open"
    ARTEFACT = "artefact"
    HUMAN_MADE = "human-made"
    KNOWN_NATURAL = "known natural class"
    NEW_NATURAL = "candidate new natural phenomenon"
    UNEXPLAINED = "unexplained"


_REALITY_PHRASE = {
    RealityTier.HIT: "a hit",
    RealityTier.VETTED: "a vetted hit",
    RealityTier.CONFIRMED: "a confirmed signal",
    RealityTier.CELESTIAL: "a confirmed celestial signal",
    RealityTier.UNEXPLAINED: "an unexplained signal",
}

_CONTENT_PHRASE = {
    ContentTier.NONE: "no content analysis",
    ContentTier.SYMBOLS: "discrete symbols recovered",
    ContentTier.STRUCTURE: "structure beyond the natural-model null",
    ContentTier.LAYOUT: "layout recovered",
    ContentTier.ANCHORED: "content anchored to shared knowledge",
    ContentTier.HYPOTHESIS: "a hypothesised reading",
}

# Words that no output may use below the tier shown. The renderer enforces this.
FORBIDDEN_BELOW = {
    "message": ContentTier.LAYOUT,
    "decoded": ContentTier.LAYOUT,
    "means": ContentTier.HYPOTHESIS,
    "translation": ContentTier.HYPOTHESIS,
}
NEVER = ("alien", "extraterrestrial intelligence", "artificial origin", "technosignature detected")


def reality_phrase(tier: RealityTier | int) -> str:
    return _REALITY_PHRASE[RealityTier(int(tier))]


def content_phrase(tier: ContentTier | int) -> str:
    return _CONTENT_PHRASE[ContentTier(int(tier))]
