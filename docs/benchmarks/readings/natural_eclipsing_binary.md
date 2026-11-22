<!-- What the workbench was given: Box-shaped periodic eclipses (P = 3.7 d, 2% deep, 3 h long) at 30-minute sampling. It was told nothing else. -->

# Reading of natural/eclipsing_binary

**Content tier 1 of 5: discrete symbols recovered.** The signal takes discrete states, so it can be written down as a sequence of symbols. The sequence shows no order beyond what a short-memory or repeating natural process produces.

Input: Series, 4,320 elements.

## Evidence

- (tier 1, discrete levels, p = 0.005) All but 0.3% of the points lie within their error bars of one of 2 levels, which are at least 10 error bars apart; the series changes level 48 times.

## Notes

- The sequence repeats every 355 symbols (99.4% of symbols are explained by one repeating cycle). A repeating cycle is a clock: it carries no information beyond the cycle itself.
- The cycle has 4 changes of state: nothing further to analyse.

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| discrete levels | 0.005 | yes | 2 levels, 10.2 error bars apart at the closest, with 0.3% of points off every level; p = 0.00498 against the drw null (200 decoys) |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
