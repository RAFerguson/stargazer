# Do the bursts of two repeating fast radio bursts keep time?

**Question.** Fast radio bursts are millisecond flashes of radio waves from other galaxies; some sources
repeat. If the bursts of a repeater fell on the ticks of a clock, or if the gaps between them were whole
multiples of one unit, that would be the kind of timing a rotating star or a deliberate transmitter
produces, and it could be turned into a sequence of symbols. Do they?

**Answer.** No. Neither source shows a clock or a unit of spacing in the ranges searched.

| Source | Bursts | Median gap | Clock test (range searched) | Spacing test (range searched) |
|---|---|---|---|---|
| FRB 20220912A, Nançay, 68 scans over 186 days | 696 | 65 s | probability 0.71 (6.6 to 131 s) | probability 0.43 (0.45 to 14.4 s) |
| FRB 20220529, FAST, over 798 days | 643 | 66 s | probability 0.10 (6.7 to 134 s) | probability 0.53 (0.39 to 12.4 s) |

A hit needed a probability at or below 0.0025 (0.01 shared over four tests). The smallest is 0.10.

## How

Both tests come from `stargazer.detect.intervals`. Each compares the data with 999 surrogates that keep the
burst rate and its slow changes, burst storms included, and destroy only fine timing. For the Nançay source
the surrogates respect the telescope's scan log, so no surrogate burst falls where nobody was observing. The
sources, tests, ranges and hit rule were written to the trials ledger and published before either burst
list was tested (search id `e59e9fad…`); the code was checked beforehand on synthetic bursts only.

## What this does and does not mean

- No clock with a period between about 7 and 130 seconds, and no spacing unit between about 0.4 and 14
  seconds, organises these bursts. Faster periods (milliseconds to seconds, where a neutron star's rotation
  would lie) were not searched: the default range of the test is tied to the typical gap between bursts.
- The Nançay times are not corrected for dispersion to infinite frequency, and the FAST list is a selected
  subset without exposure information. Neither matters at the tens-of-seconds scale searched here; both
  would matter for a millisecond search.
- FRB 20180916B, which has a known 16.35-day cycle of activity, was left out on purpose: the public
  compilation's observing windows are incomplete, and a timing test without honest exposure would not mean
  much. A cycle of activity is also a question about rate, not about a clock.

## Data

- FRB 20220912A: Konijn et al. 2024, MNRAS 534, 3331 (doi:10.1093/mnras/stae2296); data deposit
  doi:10.5281/zenodo.13889738 (CC-BY-4.0).
- FRB 20220529: Li et al. 2026, Science 391, 280 (doi:10.1126/science.adq3225); data deposit
  doi:10.5281/zenodo.17446358 (CC-BY-4.0).

## Files

- `run_timing.py`: the registered script.
- `summary.json`: the numbers above, with time standards and full citations.
