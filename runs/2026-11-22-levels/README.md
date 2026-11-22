# Do any of the quasars sit on discrete brightness levels?

**Question.** A source whose brightness is switched between a few fixed levels is the simplest coded signal a
slowly sampled light curve could carry. Does any of the 1,210 quasars in the first sample look like that?

**Answer.** No more than chance. Three of 1,210 real quasars pass the test. So do two of their 1,210 decoy
twins, which are simulations of pure red noise. The real hit rate is the false-alarm rate.

| | Searched | Two or more resolved levels | Passed (probability at or below 0.005) |
|---|---|---|---|
| Real quasars | 1,210 | 3 | 3 |
| Decoy twins | 1,210 | 3 | 2 |

## The test

For each object the g, r and i light curves were merged and averaged by night. The test
(`stargazer.decode.symbolise.discreteness_test`) asks for two to six levels at least eight error bars apart,
with at most 10% of the nightly points further than three error bars from every level, and then asks how
often red noise simulated from the object's own fit looks as discrete (200 simulations). Real objects and
decoy twins went through it under hidden labels. The definition, the settings and the rule were written to
the trials ledger and published before any object was tested (search id `f591a2bb…`).

## The three real objects that passed

| Object | Nights | Levels (standardised units) | Level separation (error bars) | Points off every level |
|---|---|---|---|---|
| gaia:3605204722313354880 | 346 | −0.51, 3.25 | 8.7 | 8.1% |
| gaia:3620681516705461376 | 454 | −0.43, 5.57 | 10.2 | 6.8% |
| gaia:3626277068817735808 | 375 | −0.50, 2.04 | 8.8 | 8.0% |

All three sit at the edge of the definition: the separations are barely above the minimum of eight error
bars and the off-level fractions are close to the limit of 10%. In each, one level holds nearly all the
points and the other a few faint ones. The two decoys that passed look the same (separations 8.2 and 8.5,
off-level fractions 5.4% and 5.9%). A few faint nights against a steady level is what red noise, or a few
poor measurements, produces at this rate.

## What this does and does not mean

- No quasar in this sample shows brightness keyed between resolved levels on timescales of nights to years.
- The test needs levels at least eight error bars apart. Keying shallower than that, or faster than one
  night, would not be seen.
- The smallest probability 200 simulations can give is 1 in 201, so "passed" here means only "rarer than
  1 in 200 for this object's own noise". The decoys show that noise alone produces about two such passes
  in a sample of this size.

## Files

- `run_levels.py`: the registered script.
- `summary.json`, `results.json`: counts and every subject's result.
- `sealed_commitment.txt`, `sealed_proof.json`: the hidden label table's hash, and the table.
