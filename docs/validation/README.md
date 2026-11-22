# Validation runs

Measurements of the search's own behaviour on simulated data. Each can be repeated with the script beside it.

| What | Script | Result (22 November 2026) |
|---|---|---|
| Does the periodicity search call pure red noise significant at the stated rate? 1,000 simulated two-band objects at survey cadence, final configuration | `null_calibration.py 1000 two` | 41 at or below 0.05 (50 expected); 7 at or below 0.01 (10 expected). Before the search was confined to periods of 40 days and longer: 52 and 6 |
| How good is the fitted-tail extrapolation beyond the simulated nulls? 40,000 simulated nulls for one sampling pattern | `tail_extrapolation_check.py 40000` | It overstates the probability by a factor of 4 at 0.01, 10 at 0.001 and 24 at 0.0001 (`tail_extrapolation_check.json`). It errs towards less significance |
| Does the decoding workbench overclaim? 33 benchmark signals, four seeds | `stargazer bench --seed N` | 32 at the justified tier, 1 short, 0 overclaims, for seeds 0 to 3 (`../benchmarks/`) |
| Does the whole blind pipeline promote noise? 270 simulated objects, their decoy twins and 270 planted signals | `stargazer.pipelines.quasar.run` on simulated objects | 0 of 270 pure-noise objects and 0 of 270 decoys promoted |
| Does it behave on real quasars? 47 ZTF quasars in a field outside every registered sample (the engineering check) | `stargazer.pipelines.quasar` on HEALPix order 4, pixel 600 | 2 of 47 at or below 0.05 (2.35 expected); none promoted; 9 of 12 planted signals recovered. This check is what exposed the lunar month |
