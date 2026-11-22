"""Error control that survives an open-ended search: online false-discovery-rate control.

A project like this never stops testing. Each new source, each new method and each new data release adds
hypotheses, and no one knows in advance how many there will be. Ordinary corrections need the total number of
tests; this one does not. The LORD++ rule (Javanmard & Montanari 2018; Ramdas, Yang, Wainwright & Jordan 2017)
gives every test its own significance level, decided before its p-value is seen, from a fixed budget ("alpha
wealth") that is spent by testing and earned back by discoveries. Whatever the number and order of tests, the
expected fraction of false discoveries among all discoveries stays at or below `alpha`, provided the p-values
are valid and the null p-values are independent of each other. Under dependence it controls a slightly weaker
quantity; the trials ledger records every test so that anyone can recompute with a stricter rule.

The state is fully determined by the sequence of past decisions, so it can be rebuilt from the ledger at any
time and cannot be quietly reset.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# gamma_t proportional to log(max(t, 2)) / (t * exp(sqrt(log t))); the constant makes the infinite sum equal one
# (Javanmard & Montanari 2018, section 3). A sum below one would only make the procedure more conservative.
GAMMA_CONSTANT = 0.07720838


def gamma(t: np.ndarray | int) -> np.ndarray:
    """The spending sequence: the share of wealth allotted to the t-th test after a (re)start, t = 1, 2, ..."""
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    ok = t >= 1
    out[ok] = GAMMA_CONSTANT * np.log(np.maximum(t[ok], 2.0)) / (t[ok] * np.exp(np.sqrt(np.log(t[ok]))))
    return out


@dataclass
class OnlineFDR:
    """LORD++ at level `alpha` with initial wealth `w0` (default alpha / 2)."""

    alpha: float = 0.05
    w0: float | None = None
    n_tests: int = 0
    rejections: list[int] = field(default_factory=list)   # 1-based indices of the tests that were discoveries

    def __post_init__(self) -> None:
        if self.w0 is None:
            self.w0 = self.alpha / 2.0
        if not 0 < self.w0 <= self.alpha:
            raise ValueError("initial wealth must lie in (0, alpha]")

    def next_level(self) -> float:
        """The significance level of the next test. It depends only on past decisions, never on the next p."""
        t = self.n_tests + 1
        level = float(gamma(t)) * self.w0
        if self.rejections:
            first = self.rejections[0]
            level += (self.alpha - self.w0) * float(gamma(t - first))
            if len(self.rejections) > 1:
                level += self.alpha * float(gamma(t - np.array(self.rejections[1:])).sum())
        return level

    def test(self, p: float) -> tuple[bool, float]:
        """Test one hypothesis. Returns (discovery, the level it was tested at)."""
        if not 0.0 <= p <= 1.0:
            raise ValueError("a p-value must lie in [0, 1]")
        level = self.next_level()
        self.n_tests += 1
        reject = p <= level
        if reject:
            self.rejections.append(self.n_tests)
        return reject, level

    @classmethod
    def from_history(cls, decisions: list[bool], alpha: float = 0.05, w0: float | None = None) -> "OnlineFDR":
        """Rebuild the state from the recorded sequence of decisions."""
        state = cls(alpha, w0)
        state.n_tests = len(decisions)
        state.rejections = [i + 1 for i, d in enumerate(decisions) if d]
        return state
