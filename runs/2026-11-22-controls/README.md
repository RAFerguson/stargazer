# Control objects, 22 November 2026

Two objects were searched one at a time with the registered periodicity configuration and the same vetting
gates as the blind quasar run. The choice of objects, the configuration and the promotion rule were written
to the trials ledger and published before either object was searched (search id `ef8ac427…`).

| Object | Data | Best period | Probability against red noise | White-noise test, for contrast | Outcome |
|---|---|---|---|---|---|
| PG 1553+113 | Fermi-LAT weekly flux, 936 points, 18.2 years | 773 days | 0.09 (10 of 111 simulated nulls matched it) | 7 × 10⁻³⁶ | not promoted |
| KIC 8462852 | Kepler long cadence, 1,393 nightly means, 4.0 years | 79 days | 1.0 | 0.12 | not promoted |

## PG 1553+113

This blazar has a published claim of a gamma-ray cycle of about 2.2 years (about 800 days). Given the same
public data and no hint, the search returns 773 days as its best period, so it recovers the claimed cycle.
It does not find it significant: red noise fitted to this light curve produces a periodic model at least as
good 9% of the time, once the search over 817 trial periods is counted. The cycle is steady from year to
year (removing any one of 19 years moves its phase by 0.01 cycles at most), but the first 70% of the data,
searched alone, does not find it.

The contrast is the point. The textbook test puts the same peak at 7 × 10⁻³⁶. The difference between that
number and 0.09 is the difference between assuming white noise and fitting the noise the object has.

This is a statement about these data and this noise model, not a verdict on the object. A cycle of this
length will need many more years, or independent data, to settle.

## KIC 8462852 (Boyajian's star)

A strange object, with deep irregular dips, that is not known to be periodic. The search finds nothing: its
best periodic model is matched by every simulated null. The dips are absorbed by a noise model with a
correlation time of 2.7 days. The pipeline does not turn strangeness into a period, which is what this
control was for.

The decoding workbench was also run on both light curves, as an unregistered check. Both stay at content
tier 0: neither sits on discrete levels.

## Files

- `run_controls.py`: the registered script.
- `summary.json`: the numbers above.
- `dossiers/`: one record per object, in JSON and Markdown.
