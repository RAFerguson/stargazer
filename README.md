# Stargazer

Stargazer searches public astronomical data for signals that natural models do not explain, and works out how
much structure an unknown signal carries. It does not assume that anything is there. Its job is to give the
question a fair, independent test and to report the outcome honestly, including when the outcome is nothing.

**Status, 22 November 2026: built and run on real data.** The search engine, the decoding workbench and the
benchmark exist and are tested. Four searches have been registered, run and recorded.

## What happened when it was run

Four searches were registered in the ledger, published, and then run. All four are recorded
([`runs/`](runs/)).

**1. Do any of 1,210 quasars vary with a strict period?** ([report](runs/2026-11-22-quasar-pix2-100/))

| | Searched | Promoted |
|---|---|---|
| Real objects (Gaia quasar candidates with ZTF light curves, one 215-square-degree patch) | 1,210 | 1 |
| Decoy twins (simulations of each object's own noise; anything promoted is a false alarm) | 1,210 | 0 |
| Planted signals (known signals added to real light curves) | 484 | 41 at the right period |

The one object promoted turned out to be a star. It showed a clean 80-day cycle in all three bands and passed
every gate. Follow-up found a period of 0.49622 days in the unbinned data, and Gaia's own catalogue lists
the source as an RR Lyrae star with a period of 0.49624 days: a half-day pulsation seen once a night beats
at 80 days. Gaia's classifier had given it a 96% chance of being a quasar. So the result is: **no quasar in
the sample has a strict period that survives vetting; the search found one genuinely periodic object blind,
and it is a known kind of star; nothing false was promoted.** An alias gate was added afterwards and fails it.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/first-search-calibration-dark.png">
  <img alt="Share of real objects and of decoy twins at or below each probability. Decoy twins follow the line a calibrated null gives. Real objects have an excess at small probabilities." src="docs/figures/first-search-calibration-light.png">
</picture>

The figure shows the other finding. Decoys, which are the noise model, sit on the calibrated line. A few
percent of real quasars do not: slow trends and annual systematics that a single damped random walk lacks.
The population recalibration absorbed them (41 of the 46 strays also failed a vetting gate), at a price in
sensitivity: a 0.2-magnitude sine wave at 150 days is recovered 72% of the time, sawtooth and pulse shapes
rarely. A better noise model is the next thing to build.

**2. What does the search say about famous objects?** ([report](runs/2026-11-22-controls/))
For the blazar PG 1553+113 it returns 773 days as the best period in 18 years of Fermi data, matching the
published claim of about 2.2 years, and puts its probability under red noise at 0.09. A white-noise test
gives 7 × 10⁻³⁶ for the same peak. For Boyajian's star it finds no period. Neither is promoted.

**3. Does any of the same quasars sit on discrete brightness levels?** ([report](runs/2026-11-22-levels/))
Three of 1,210 pass the test, against two of their 1,210 decoy twins: the false-alarm rate.

**4. Do the bursts of two repeating fast radio bursts keep time?** ([report](runs/2026-11-22-frb-timing/))
No clock and no unit of spacing in 696 bursts of FRB 20220912A or 643 of FRB 20220529: the smallest
probability is 0.10 against a threshold of 0.0025.

The ledger now holds 4 searches, 2,424 objects, 383,798 trial hypotheses and one formal discovery, which is
the star's periodicity.

## What it does

**Finds periodic signals without being fooled by red noise.** Quasars flicker in a way that imitates
periodicity, and every optical "periodic quasar" tested on later data has failed. Stargazer fits the red noise
each object has, searches for a periodic signal of any shape (sine, sawtooth, narrow pulse), and measures
significance by running the identical search on simulated nulls. On 1,000 simulated pure-noise objects it
called 4.1% significant at the 5% level and 0.7% at the 1% level ([validation](docs/validation/)). A standard
white-noise test calls almost all of them significant.

**Keeps itself blind.** Every real object goes into the search alongside a decoy twin (a simulation of its own
fitted noise) and, for a fraction, a copy with a known signal planted in it. Labels are hidden until the
analysis is final. Promoted decoys count the false alarms; recovered plants measure the sensitivity. No
threshold can be tuned on the real answer, because the search cannot see it.

**Vets every hit.** Seven gates ask what an impostor would fail: is the period one the Moon or the year would
impose; does the sampling pattern repeat at it; does it survive removing any year of data; do the bands
agree; do the early data alone find it and predict the late data; does the search do worse when the time axis
is warped; is it the beat of a much faster variation against the nightly cycle.

**Says how much of an unknown signal can be read.** The decoding workbench takes a bit stream, a light curve
or a list of event times and climbs a ladder: discrete symbols, structure, layout, anchor, reading. It finds
raster pictures by exact code length, frames and sync words, parity checks, shift registers, counting, primes
and arithmetic, and it reports everything it tried. On 33 benchmark signals with known answers it placed 32 at
exactly the justified tier and overclaimed on none, on the development seed and on three held-out seeds
([scorecard](docs/benchmarks/scorecard.md), [example readings](docs/benchmarks/readings/)).

**Keeps count.** Every search is registered in a hash-chained ledger before it runs and recorded after, so
the number of things ever tested is known, and discoveries are controlled across the whole history of the
project, not one run at a time.

## The two ladders

Nothing is ever reported as "a message". Every result is placed on two independent ladders.

| Reality: has it survived scrutiny? | Content: how much has been recovered? |
|---|---|
| 0 a hit | 0 nothing found |
| 1 a vetted hit | 1 discrete symbols |
| 2 confirmed in independent data, or by a prediction that came true | 2 structure beyond a short-memory process |
| 3 tied to a position on the sky | 3 layout, framing or coding recovered |
| 4 unexplained by the natural models tried | 4 content anchored to shared knowledge |
| | 5 a tested reading, labelled as a hypothesis |

A pulsar is low on the first ladder (explained) and can still be analysed on the second. A raster picture found
in a test file is high on the second and says nothing about the first.

## Using it

```
python -m venv .venv && .venv/bin/pip install -e ".[dev,figures]"

stargazer sources check                        # which archives answer from here
stargazer fetch ztf 196.38754 -10.555389 -o pg1302.json
stargazer search pg1302.json                   # periodicity search and gates on one object
stargazer decode bits.txt                      # content ladder for a file of 0s and 1s
stargazer bench                                # the benchmark scorecard
stargazer quasar register --name my-run        # write a blind run's definition to the ledger ...
stargazer quasar run --name my-run --search-id <id> --out runs/my-run     # ... publish it, then run
stargazer ledger                               # what has been tested so far
pytest                                         # 222 tests; add -m "network or slow" for 15 more
```

## What is here

| Path | What it is |
|---|---|
| `stargazer/nulls`, `detect`, `vet`, `trust` | Red-noise model, periodicity and event-timing searches, vetting gates, blind analysis, error control |
| `stargazer/pipelines/quasar.py` | The blind quasar search, end to end |
| `stargazer/structure`, `decode` | Sequence statistics and the decoding workbench |
| `stargazer/bench` | 33 benchmark signals with known answers, damage channels and the scorecard |
| `stargazer/sources` | ZTF, Gaia, Kepler, Fermi and fast-radio-burst connectors, a HATS catalogue reader, a download cache |
| `stargazer/ledger`, `report` | Trials and prediction ledgers, dossiers, figures |
| [`ledger/`](ledger/) | The project's own record: every search registered and recorded |
| [`runs/`](runs/) | Each run's configuration, sealed-label proof, results table, dossiers and report |
| [`docs/PLAN.md`](docs/PLAN.md), [`docs/INNOVATION_AGENDA.md`](docs/INNOVATION_AGENDA.md) | The plan and the list of tools to invent, written before the code |
| [`docs/research/`](docs/research/) | The research behind the plan: a report, eight topic notes, 516 checked citations |

## What it cannot do

- It cannot tell content from noise when content is compressed or encrypted. Nothing can.
- It cannot read meaning without an anchor. The anchor library is small and fixed: counting, primes,
  arithmetic, four constants, ASCII text, repetition.
- From the ground it cannot find periods near a month, near a year, or shorter than 40 days: the Moon and the
  seasons put their own signals there. An engineering check on real ZTF quasars showed a 29.5-day lunar line
  in 60% of them.
- A period seen for only a few cycles cannot be told from a wandering one, and the search says so by losing
  sensitivity there.
- The noise model is a single damped random walk. Real quasars are more varied. The population recalibration
  and the warped-time test guard against this; they do not remove it.
- Radio data (filterbank spectra, pulsar timing, fast radio bursts) are planned and not yet connected.
- One benchmark item is still missed: text hidden under an additive scrambler.

## Credit

Data: the ZTF light-curve catalogue served by IRSA (doi:[10.26131/IRSA598](https://www.ipac.caltech.edu/doi/10.26131/IRSA598))
and its HATS copy on public cloud storage; Gaia DR3 (ESA/DPAC); the Fermi-LAT Light Curve Repository; Kepler
(NASA/MAST); burst lists from Konijn et al. 2024, Pastor-Marazuela et al. 2021 and Li et al. 2026 (Zenodo). Several sequence-analysis functions and the prediction ledger are adapted from the author's
SeaQuest project; each file says what was changed.
