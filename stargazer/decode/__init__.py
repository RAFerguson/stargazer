"""Working out the structure and content of an unknown signal.

`workbench.analyse` is the entry point. The modules it draws on can each be used alone:

* `symbolise`: discrete levels in a light curve, with honest per-point probabilities;
* `layout`: raster and volume geometry from exact code lengths;
* `codes` and `gf2`: frames, sync words, linear block codes and shift registers;
* `anchors`: a fixed library of shared knowledge (counting, primes, arithmetic, constants, character codes).
"""
from . import anchors, codes, gf2, layout, symbolise, workbench
from .workbench import Reading, Settings, analyse, to_markdown

__all__ = ["Reading", "Settings", "analyse", "anchors", "codes", "gf2", "layout", "symbolise", "to_markdown",
           "workbench"]
