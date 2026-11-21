# Stargazer innovation agenda: the custom tools to invent

**Status:** draft v0.1, 2026-11-21. Research phase only: nothing here is built. Companion to `docs/PLAN.md`. The plan
describes what can be assembled from published methods and from the owner's private sister project. This document
describes the tools for analysing and translating unknown data that do not exist yet, in software, mathematics and
experimental practice, and says how each would be judged.

---

## 0. How to read this document

Section 1 names the ten walls the research phase found: limits that more engineering on existing methods will not
get through. Sections 2 to 4 describe the inventions, grouped as tools (T), mathematics and statistics (M), and
experimental practice (X). Each item has the same shape: why it is needed, what has to be invented, how to get
there, what success looks like in numbers, and a checkpoint that says when to stop or change course. Where an item
grows out of the sister project, its lineage is named by function. Section 5 bundles the inventions into five grand
challenges. Section 6 is the portfolio view and Section 7 the sequencing.

Necessity tiers:

- **Required**: the project cannot meet its goals without it.
- **Enabling**: the project works without it but is far weaker.
- **Visionary**: a real advance if it works; fund the first step and keep the option open.

Horizons: near (within the first six months of building), mid (six to eighteen months), far (beyond).

No hardware inventions are proposed. Stargazer works on public archives; the nearest thing to hardware is the
follow-up pathway in X6.

---

## 1. The ten walls

1. **The red-noise wall.** Quasars, and many other sources, wander randomly in a way that imitates a few cycles of
   a period. Every optical-survey list of periodic quasars so far has been dominated by false positives, and
   baselines are too short to settle most claims.
2. **The sampling wall.** Survey data are irregular, gappy and noisy, with an error bar on every point. Sequence
   tools (entropy, grammar, decipherment) assume clean, evenly indexed symbols. Applied naively, they measure the
   observing schedule.
3. **The noise-equivalence wall.** A well-compressed or encrypted message is statistically indistinguishable from
   noise. Structure tests can only ever find senders who chose to be found or were careless.
4. **The haystack wall.** Billions of sources, many tests each, and a search that never ends. Without a count of
   every look, no probability means anything; with one, most "discoveries" dissolve.
5. **The contamination wall.** Satellites, radio interference and instrument systematics produce almost all
   anomalies. Human technology is the most common technosignature in the data.
6. **The natural-structure wall.** Natural processes are structured. Pulsar nulling is not random, burst sequences
   have memory, and artificial sign systems with no language behind them pass every published "is it language"
   test. Organised does not mean artificial, and no reference exists that places a signal among known generators.
7. **The format wall.** Going from a physical signal to bits, frames and layout for an unknown format is done by
   hand, differently each time. The communications-engineering tools that exist assume strong human-made signals
   and a person at the controls.
8. **The grounding wall.** Meaning has never been recovered without an outside anchor: a bilingual text, a related
   language, a shared world. A one-way signal offers only mathematics, physics, pictures and its own properties.
9. **The single-shot wall.** The most intriguing events happen once and cannot be re-observed, and a one-way
   sender cannot be questioned. Confirmation needs recurrence or a prediction that comes true.
10. **The access wall.** Public products are reduced and averaged; the time resolution that would carry fast
    modulation is rarely released; archives change and lose funding.

---

## 2. Tools

### T1. Sampling-aware symboliser

**Tier:** Required. **Horizon:** near. **Walls:** 2. **Lineage:** two principles from the sister project:
quantisation chosen by minimum description length under a reversibility constraint, and a test for whether a
signal's states are discrete or graded. The principles carry over; the code does not.

**Why.** Every content test needs a symbol sequence. The sister project's sequences are certain and evenly indexed.
Astronomical samples are uncertain and irregular, and thresholding them makes symbols whose errors and gaps pass
for structure or hide it.

**What to invent.** A symboliser whose output is a lattice, not a string: for each time slot, a probability for
each symbol and for "not observed". Alphabet size, bin edges and slot length are chosen by minimum description
length with the measurement noise inside the likelihood. Every downstream statistic is computed as an expectation
over the lattice, so symbol uncertainty reaches the final number.

**How.** (1) Write the description-length objective with a noise model inside it. (2) Slot irregular samples with a
hidden-Markov forward–backward pass. (3) Propagate by sampling sequences from the lattice. (4) Validate on synthetic
symbol sources observed through real sampling patterns and error bars taken from the archives.

**Success metrics.** On synthetic sources with 2 to 16 symbols seen through ZTF-like and Kepler-like sampling: the
alphabet size is recovered at least 90% of the time at a signal-to-noise ratio of 5 or more; entropy-rate estimates
fall within 0.1 bit per symbol of the truth, with 95% intervals that cover it 90 to 98% of the time; pure red noise
put through the symboliser triggers structure tests at the nominal rate.

**Checkpoint.** If no alphabet larger than two symbols can be recovered at survey cadence and realistic noise,
content analysis is restricted to the fast-timing channels, and the plan says so.

### T2. Natural-model simulator bank

**Tier:** Required. **Horizon:** near. **Walls:** 1, 6.

**Why.** A shuffle is the wrong null for a natural source. It destroys the source's own correlations, so every
quasar and every pulsar looks "structured" against it. The sister project's controls use shuffles and frequency
draws, which is right for text and wrong here.

**What to invent.** One interface, `fit(series) → posterior` and `simulate(posterior, sampling) → series`, behind
which sits a generative model for each natural class: red-noise processes for active galaxies; rotation, pulsation,
flares and eclipses for stars; nulling, mode changing, drifting sub-pulses and giant pulses for pulsars; waiting
times, energies and activity windows for repeating bursts; scintillation. Simulations pass through the same
sampling, noise and symboliser as the data, which makes them the right null for every content test.

**How.** Wrap the maintained Gaussian-process libraries for red noise. Implement the pulse-sequence models from the
literature. Check each class against real catalogued objects with posterior predictive tests.

**Success metrics.** Posterior predictive p-values are uniform on held-out real objects of each class. Content
tests run on real natural objects raise alarms at no more than the nominal rate.

**Checkpoint.** Any class that fails its predictive checks has content tests disabled until it passes.

### T3. Reference-generator bank and placement

**Tier:** Required. **Horizon:** near to mid. **Walls:** 6. **Lineage:** the sister project's
typology bank (nearest neighbours on a handful of statistics).

**Why.** A structure statistic means nothing alone. Its meaning comes from where it falls among known generators,
measured with the same code at the same length and noise.

**What to invent.** A bank of generator classes (natural processes, human protocols and file formats, natural
languages, constructed messages, compressed and encrypted data, noise), each with its structure statistics as a
function of sequence length and noise. A placement routine returns a calibrated distribution over classes with an
explicit "none of these", replacing nearest-neighbour distance.

**Success metrics.** On held-out generators: the right class is ranked first at least 85% of the time at 10³
symbols or more; calibration error is under 5%; an unseen class is flagged as "none of these" with an area under
the curve of 0.9 or better. The bank must also reproduce the known negative result: it cannot separate a real
language from an organised non-language, and reports them as one class.

**Checkpoint.** If placement is no better than inventory size and entropy alone, ship those two numbers and drop
the classifier.

### T4. Universal modulation search

**Tier:** Enabling. **Horizon:** mid. **Walls:** 7, 3.

**Why.** Searches look for one shape at a time: a tone, a pulse, a period. Once a carrier is found, nothing
systematically asks whether anything is impressed on it. Communications engineering has a mature chain for this
(carrier and symbol-rate estimation, timing recovery, blind equalisation, constellation clustering), but its
modulation classifiers only choose among known human modulation types, and its tools are interactive, assume strong
signals and are mostly under the GPL. They can generate hypotheses; they cannot judge them.

**What to invent.** Given a carrier (a tone in a spectrogram, a periodic pulse train, a steady light curve),
demodulate under each family in turn (on–off, amplitude levels, frequency shifts, phase where available, pulse
position and interval, pulse width, polarisation) over a range of symbol clocks. Score each by the description
length saved relative to the fitted natural model, corrected for the number of family-and-clock hypotheses tried.

**How.** Start with pulse trains, where intervals and amplitudes are already measured; the sister project already
analyses interval patterns of this kind. Then narrowband tones in Breakthrough Listen spectrograms, with modulated
tones injected by setigen. Voltage data last.

**Success metrics.** On injected signals the family and symbol clock are recovered at least 90% of the time above a
stated signal-to-noise ratio; the false-alarm rate on the natural and human-made collections is nominal; Voyager
1's telemetry sidebands are reported as modulation on a carrier.

**Checkpoint.** If public spectrogram products are too coarse in time to show any modulation, the tool is limited
to pulse trains and the limit is published.

### T5. Layout finder with calibrated significance

**Tier:** Required for decoding. **Horizon:** near. **Walls:** 7.

**Why.** Finding how a bit stream folds is the classic first step with a designed message; the Arecibo message's
1,679 bits fold into 23 columns. Scanning candidate widths for a drop in complexity was published in 2023 and
ranked the right size first for about 80% of 102 test images. A 2026 follow-up found that the raw measure picks
the wrong width for Arecibo itself (11, not 23) until it is normalised for shape and calibrated against shuffles,
which recovers 23 at a search-corrected probability of 0.005. That work is the starting point. What is still
missing is a reusable tool that also handles unknown offsets, headers and bit order, missing bits, more than two
dimensions, and reports how many layouts fit nearly as well.

**What to invent.** A scan over folding hypotheses (row length, offset, header and footer length, bit order, two
or three dimensions, interleaving), each scored by several compressibility measures and calibrated against
shuffled and Markov-matched surrogates. It tolerates bit errors and missing bits and returns ranked layouts with an
ambiguity count.

**Success metrics.** The Arecibo and Cosmic Call layouts and the 256 × 256 bitmap of "A Sign in Space" are the top
hypothesis from bits alone. With 3% of bits inverted the correct width is in the top three at least 93% of the
time, matching the published figure, and with 10% inverted it beats the published 44%. On random and Markov
surrogates, false layouts are claimed at the nominal rate. Synthetic three-dimensional rasters are recovered.

**Checkpoint.** If error tolerance cannot reach 1%, the tool is documented as needing clean bits and stage 1 of the
workbench must supply them.

### T6. Blind code recogniser

**Tier:** Enabling. **Horizon:** mid. **Walls:** 3, 7.

**Why.** The noise-equivalence wall has a crack. A sender can make the payload look like noise, but a receiver
still has to lock on and correct errors, so framing, synchronisation and error-correcting redundancy must be
there. They are the envelope, and the envelope is detectable even when the letter is not. Methods for recognising
codes blind are published (rank deficiency and low-weight parity checks for block codes, Berlekamp–Massey for shift
registers, parity-check methods for interleavers). The research phase found no maintained open-source package that
implements them.

**What to invent.** Open implementations of blind recognition for linear block and convolutional codes, scramblers
built on shift registers, interleavers and synchronisation words, tolerant of bit errors, each returning a
significance against surrogates.

**Success metrics.** Code length and dimension are identified for Hamming, BCH and Reed–Solomon coded synthetic
streams at bit-error rates up to 10⁻² given a stated number of codewords; shift-register polynomials are recovered;
no code is claimed on the control collection beyond the nominal rate.

**Checkpoint.** If modern near-capacity codes prove undetectable at any realistic length, say so and keep the
classical codes, which a sender wishing to be understood would prefer.

### T7. Anchor search with a declared hypothesis budget

**Tier:** Required for decoding. **Horizon:** near to mid. **Walls:** 8.

**Why.** Meaning can enter only through anchors a sender would share: counting, arithmetic, physical constants,
pictures, the signal's own properties. Hunting for them without counting the hypotheses is numerology.

**What to invent.** A registry of anchor families, each with its explicit number of hypotheses: sequences (counting,
primes, squares, powers of two) crossed with number representations (unary, binary, positional bases, both bit
orders); arithmetic identities that must hold across every instance; ratios tied to physics, such as multiples of
the hydrogen line or the ordering of the lightest elements; shapes in recovered rasters. A search engine returns
matches with significance corrected for the whole registry, and the rate at which the same registry "finds" anchors
in surrogates.

**Success metrics.** Planted anchors in the message collection are found (the numbers one to ten and the atomic
numbers in the Arecibo message, the operators in a Lincos-style message); zero anchors are claimed in 10⁴ surrogate
streams at the declared level.

**Checkpoint.** If the registry cannot be made small enough to leave any power after correction, it is cut back to
counting and arithmetic identities.

### T8. Decoding hypothesis manager

**Tier:** Required. **Horizon:** near to mid. **Walls:** 7, 8. **Lineage:** the sister project's decipherment
pipeline (tiered readings, consistency checks, translatability statement).

**Why.** A decoding is a chain of choices, and analysts keep the chain they like. The sister project's pipeline follows one
chain for text. Nothing keeps rival chains alive, scores them on one scale, and reports how many fit equally well.

**What to invent.** A ledger of decoding hypotheses (modulation, alphabet, layout, code, mapping), each with a
two-part description length and a held-out score, with acceptance only by a declared margin over rivals. It reports
the ambiguity count with every reading and feeds the translatability statement.

**Success metrics.** On the message collection the true chain ranks first in at least 80% of cases; the ambiguity
count is higher where the top hypothesis is wrong; no reading is accepted on the control collection.

**Checkpoint.** If description lengths from different stages cannot be put on one scale, rank within stages and
report the chain as a set of stage-wise winners.

### T9. Benchmark collections and blind-challenge harness

**Tier:** Required. **Horizon:** near. **Walls:** all. **Lineage:** the sister project's two benchmarks: a generated
communication system with ground truth, and known text rendered into an unknown script.

**Why.** Without known answers, no detector or decoder can be trusted, and none can be compared.

**What to invent.** The four collections of plan Section 9.1 (natural, human-made, message, control) with their
generators; channel simulators for dispersion, scattering, scintillation, Doppler drift, radiometer and photon
noise, and sampling masks lifted from real sources; scoring; and a blind mode in which ground truth is sealed under
a hash commitment until scoring.

**Success metrics.** Every collection regenerates from its seeds; every Stargazer detector and decoder has a
scorecard; at least one outside tool is scored on the same harness for comparison.

**Checkpoint.** If the natural collection's simulated members are distinguishable from its real members by a
simple classifier, the generators are not realistic enough to calibrate anything.

### T10. Vetting gates as code, with a contaminant catalogue

**Tier:** Required. **Horizon:** near. **Walls:** 5.

**Why.** Verification frameworks exist as papers and checklists. They are applied by hand, late, and differently
each time.

**What to invent.** The gates of plan Section 7.6 as versioned, executable checks, and a catalogue of contaminant
signatures: calendar and instrument periods for each survey, satellite downlink bands, known interference at each
radio site, spacecraft systematics.

**Success metrics.** On a labelled set of historical false positives, at least 95% are stopped at the correct
gate; no more than 5% of injected true signals are wrongly stopped.

**Checkpoint.** Any gate that stops more than 5% of injected signals is redesigned or demoted to a warning.

### T11. Dossier generator with evidence-constrained text

**Tier:** Required. **Horizon:** near. **Walls:** 4, 5. **Lineage:** the sister project's renderer (sentences tied to
evidence paths, entailment check, ablation test).

**What to invent.** The dossier of plan Section 8.1, generated from a candidate's records alone, with every
sentence bound to the evidence that supports it.

**Success metrics.** Every sentence passes the entailment check; removing an evidence field removes exactly the
sentences that depend on it; a dossier regenerates byte for byte from its manifest.

### T12. Cross-survey stitcher

**Tier:** Required for the quasar module. **Horizon:** near. **Walls:** 1, 10.

**Why.** Only a stitched baseline covers enough cycles of a slow signal. Surveys differ in filters, apertures and
calibration, and a wrong offset between two of them is a step that looks like low-frequency power.

**What to invent.** A joint model in which per-survey offsets, colour terms, host-galaxy dilution and extra error
terms are inferred together with the red-noise process, so their uncertainty flows into the periodicity evidence.

**Success metrics.** On simulated multi-survey light curves the offset posteriors cover the truth at nominal
rates; false-alarm calibration is unchanged from the single-survey case; where real surveys overlap in time, their
corrected measurements agree within errors.

**Checkpoint.** If offsets and long periods cannot be separated for surveys that do not overlap in time, periods
longer than the shortest single-survey baseline are reported as unconstrained.

### T13. Coincidence engine

**Tier:** Enabling. **Horizon:** mid. **Walls:** 5, 9. **Lineage:** the sister project's tests for
transitions between emitters against permutation nulls.

**Why.** A signal seen at once in two unrelated sources, or one event seen by two instruments, is either a
systematic or something important. Surveys observe neighbouring fields on the same nights, so naive coincidence
rates are set by the schedule.

**What to invent.** Coincidence statistics whose nulls are built from the actual observing schedules and shared
systematics of the instruments involved.

**Success metrics.** The false-alarm rate on real unrelated pairs is nominal; planted synchronous events are
recovered at the stated completeness.

### T14. Language-model hypothesis generator, under constraints

**Tier:** Visionary. **Horizon:** mid. **Walls:** 7, 8.

**Why.** Language models are strong pattern proposers and unreliable judges. Research found fluent, unverifiable
"decipherments" produced this way.

**What to invent.** A proposer that suggests decoding chains and anchors to the hypothesis manager and has no other
power. Its proposals are scored like any others.

**Success metrics.** On the message collection, its proposals are accepted more often than enumerated baselines,
with no accepted reading on the control collection.

**Checkpoint.** A single accepted reading on controls disables it.

---

## 3. Mathematics and statistics

### M1. Error-aware, finite-sample estimators of structure

**Tier:** Required. **Horizon:** near to mid. **Walls:** 2, 6. **Lineage:** the sister project's sequence
statistics, which use plug-in estimates.

**Problem.** Entropy and mutual-information estimates are biased upward for short sequences. On random text of 200
symbols over 26 letters, the plug-in "entropy drop" reads 1.8 bits where the truth is zero. Astronomical sequences
are short, and their symbols are uncertain as well.

**Innovation.** Estimators of entropy rate, mutual information by lag and compression complexity that combine the
standard small-sample corrections with symbol uncertainty from the lattice of T1, and report intervals whose
coverage is checked by simulation. Lag is measured in time, with a gap-aware pairing rule, not in index.

**How to achieve it.** Start from the published corrections (Miller–Madow, Grassberger, Nemenman–Shafee–Bialek,
Chao–Shen); add propagation by sampling from the lattice; replace low-order n-gram models with context-tree or
compression-based estimators for binary streams, where structure lives at long range.

**Success metric.** On synthetic Markov sources of order 0 to 3 and lengths 10² to 10⁵, 95% intervals cover the
truth 90 to 98% of the time. The estimators detect unary-coded primes and a repeated 16-bit word with 5% of bits
flipped, both of which a control that looks one symbol back misses.

### M2. Detectability bounds for modulation in red noise

**Tier:** Enabling, and decisive for upper limits. **Horizon:** mid. **Walls:** 1, 2, 3.

**Problem.** Nobody has published how much coded modulation could hide in a quasar light curve, or in any
irregularly sampled red-noise series. Without that, "we found nothing" has no content.

**Innovation.** For a given series (its sampling, its errors, its fitted red noise), the minimum modulation depth
and maximum symbol rate detectable at a stated confidence, for each modulation family of T4: in effect the capacity
of the observation as a channel. Summed over a survey, this is the upper limit the survey supports.

**How to achieve it.** Derive the bound for on–off and amplitude-level keying under a Gaussian-process noise model;
check it against injection and recovery; extend to interval and period modulation.

**Success metric.** Bounds agree with injection–recovery results within a factor of two across families. Stargazer
publishes modulation limits for ZTF quasars and for Kepler stars.

### M3. Error control for a search that never ends

**Tier:** Required. **Horizon:** near. **Walls:** 4.

**Problem.** Fixed thresholds assume a fixed number of tests. Stargazer will keep adding sources, detectors and
data releases for years, and each new look weakens every earlier claim unless the accounting is built for it.

**Innovation.** A trials ledger with error control that stays valid as tests accumulate: false-discovery control by
alpha-investing or e-values, structured over the hierarchy of sources, frequencies, detectors and time, with
allowance for correlated tests.

**Success metric.** In a simulated programme of searches added one after another, the false-discovery rate stays at
or below target. The ledger returns the corrected significance of any past hit immediately.

### M4. A translatability bound for one-way messages

**Tier:** Enabling. **Horizon:** mid. **Walls:** 8. **Lineage:** the sister project's
translatability statement, whose weights are set by hand.

**Problem.** The one rigorous theory of translation without parallel text (Goldwasser, Gruber, Kalai and Paradise
2023) bounds the error by a term that falls with the amount of common ground and the amount of data, and it needs a
prior over plausible meanings. It gives no estimator. A one-way message has no interlocutor, so common ground can
come only from anchors.

**Innovation.** An operational estimator with three measured terms: complexity (description-length gain on
held-out data), common ground (the anchored fraction, from T7), and data volume. Its output is a probability that a
given tier of reading is reachable, with the single most valuable missing requirement, fitted to cases with known
answers in place of hand-set weights.

**Success metric.** The estimator ranks the benchmark messages in their known order of difficulty and predicts
decoding success on held-out benchmark items with a calibrated probability.

### M5. Alignability to shared-knowledge structures

**Tier:** Visionary. **Horizon:** mid. **Walls:** 8. **Lineage:** the sister project's seedless
alignment of two symbol systems (entropic Gromov–Wasserstein).

**Problem.** A message that teaches arithmetic contains symbols whose pattern of co-occurrence mirrors the number
line and the operation tables. No tool tests for that mirror without first guessing a dictionary.

**Innovation.** Seedless alignment between the co-occurrence geometry of a message's symbols and candidate
structures (the number line, arithmetic tables, the ordering of the elements), with adequate permutation nulls,
finite-sample corrections and a breakdown of which symbols align.

**Success metric.** Numerals and operators are identified in Lincos-style synthetic messages above a stated
accuracy. Controls are reported as not alignable.

### M6. Calibrated origin odds

**Tier:** Enabling. **Horizon:** mid. **Walls:** 4, 5.

**Problem.** A dossier holds ten gate outcomes and several model fits. Readers will collapse them into one
impression, usually the wrong one.

**Innovation.** A Bayesian combination of gate outcomes and natural-model evidence into posterior odds for each
disposition, with the priors explicit and their influence shown, mapped onto the probability term of the Rio
Scale.

**Success metric.** On the benchmark collections, whose origins are known, stated odds match observed frequencies
within tolerance.

### M7. Population inference for periodic sources

**Tier:** Enabling. **Horizon:** mid. **Walls:** 1, 4.

**Problem.** Individual periodic-quasar claims keep failing. The better-posed question is what fraction of sources
is periodic at all; one 2025 analysis put it at no more than a few per million.

**Innovation.** A hierarchical model that estimates the periodic fraction, and the distribution of periods and
amplitudes, from the whole ensemble of per-source evidence, without naming any single source.

**Success metric.** Injected periodic fractions down to 10⁻⁴ are recovered in simulated populations of 10⁵
sources. The result is Stargazer's quasar census limit.

---

## 4. Experimental practice

### X1. A prediction ledger for astronomy

**Tier:** Required. **Horizon:** near. **Walls:** 1, 9. **Lineage:** the sister project's prediction ledger
(hash-chained commit, reveal and score, with confidence intervals on the hit rate), reused as it stands.

**Why.** Every periodic-quasar list failed when later data arrived. Forecasts committed before the data are the one
test that cannot be tuned.

**What to invent.** The practice around the tool: which forecasts are committed for which data releases, scoring
rules declared in advance and compared with the null model's forecast, and commitment hashes anchored in the public
Git history.

**Success metrics.** Forecasts are committed before each major release (Gaia DR4, the next ZTF release, the Rubin
previews) and every one is scored. The ledger's hit rate is published with its interval, whatever it is.

### X2. The rediscovery campaign

**Tier:** Required. **Horizon:** near. **Walls:** 5, 6.

**What.** A fixed list of known phenomena for each channel that the pipeline must find without being pointed at
them: pulsar periods, the periodic repeating burster, Boyajian's Star, a quasi-periodic eruption source, Voyager
1's carrier.

**Success metric.** Each is recovered at its expected significance and reaches its expected disposition.

### X3. Blind injection challenges

**Tier:** Enabling. **Horizon:** near to mid. **Walls:** 4.

**What.** Signals injected into real data by someone other than the analyst, under sealed ground truth: a second
person, or a separate automated session with no shared state.

**Success metric.** Blind detection efficiency agrees with the self-assessed completeness maps. A gap between them
measures the analyst's optimism.

### X4. Publishing null results as coverage

**Tier:** Required. **Horizon:** near to mid. **Walls:** 3, 4.

**What.** A standard, machine-readable limits file from every run (what was excluded, in how many sources, at what
completeness), and a first public note on modulation limits in quasars.

**Success metric.** Every run emits its limits file; the first note is released with code and data references.

### X5. Cross-calibration with the sister project

**Tier:** Enabling. **Horizon:** near. **Walls:** 6.

**Why.** The two projects ask the same question of different senders, and each has what the other lacks. The sister
project has communication systems with partly known structure. Stargazer will have natural non-communicative
sources in bulk.

**What.** Run the sister project's transcripts and synthetic systems through Stargazer's workbench, and
Stargazer's natural collection through the sister project's structure pipeline.

**Success metric.** The sister project's systems are placed as organised communication, with their known structure
recovered. Astronomical nulls are not reported as communication by either code base. Fixes flow back,
and the shared-library decision is made on this evidence.

### X6. A follow-up pathway

**Tier:** Enabling. **Horizon:** mid. **Walls:** 9.

**What.** A prepared route to independent observations for a candidate at R1 or above: amateur-observer campaigns,
proposal-based public telescopes, and contacts at the instrument teams, with request templates.

**Success metric.** One dry run completed on a benign target.

---

## 5. Grand challenges

Each bundles several inventions into one demonstrable result.

### G1. The honest quasar census (T2, T12, M3, M7, X1, X4)

Periodicity limits for 10⁵ to 10⁶ quasars, with completeness maps, a population-level periodic fraction, and
forecasts on the ledger for every surviving candidate. Deliverable: a public note and a data release.

### G2. Decode the solved, from the signal (T4 to T9, M1)

Every historical test message decoded end to end and blind, starting from the rawest public form: the Arecibo
bits, the Cosmic Call pages, Heller's challenge, the Busch and Reddick message, and "A Sign in Space" from its
public radio recordings. Deliverable: a scorecard, including time taken.

### G3. Organised is not artificial (T1, T2, T3, T9, M1)

A discrimination benchmark on which natural, human-made and designed signals are separated at a stated error
rate, and on which the known impossibility (language against organised non-language) is reported as impossible.
Deliverable: the benchmark collections, released for others.

### G4. What could we have seen (M2, M3, X4)

Programme-level upper limits on modulated signals across the archives Stargazer reaches, stated as the fraction of
the search space covered. Deliverable: a limits catalogue that grows with every run.

### G5. One structure toolkit for both projects (X5, M1, T8)

The shared library that both projects depend on, with the corrected estimators and null-calibration tests.
Deliverable: the library, with both projects' test suites passing on it.

---

## 6. Portfolio view

| ID | Invention | Type | Tier | Horizon | Plan phase | Main risk |
|---|---|---|---|---|---|---|
| T1 | Sampling-aware symboliser | tool | Required | near | 2 | survey cadence supports too few symbols |
| T2 | Natural-model simulator bank | tool | Required | near | 1, 2 | models too simple for real sources |
| T3 | Reference-generator bank | tool | Required | near–mid | 2, 5 | adds nothing beyond entropy and inventory size |
| T4 | Universal modulation search | tool | Enabling | mid | 5 | public products too coarse in time |
| T5 | Layout finder | tool | Required (decoding) | near | 5 | error tolerance |
| T6 | Blind code recogniser | tool | Enabling | mid | 5 | modern codes undetectable |
| T7 | Anchor search with a budget | tool | Required (decoding) | near–mid | 5 | no power left after correction |
| T8 | Decoding hypothesis manager | tool | Required | near–mid | 5 | one scale across stages |
| T9 | Benchmark collections | tool | Required | near | 0–2 | unrealistic generators |
| T10 | Vetting gates as code | tool | Required | near | 1, 4 | gates that stop real signals |
| T11 | Dossier generator | tool | Required | near | 1 | none major |
| T12 | Cross-survey stitcher | tool | Required (quasars) | near | 1 | offsets confounded with long periods |
| T13 | Coincidence engine | tool | Enabling | mid | 6 | schedule-driven false coincidences |
| T14 | Language-model proposer | tool | Visionary | mid | 5 | hallucinated readings |
| M1 | Error-aware estimators | mathematics | Required | near–mid | 2 | coverage at short lengths |
| M2 | Detectability bounds | mathematics | Enabling, decisive | mid | 5, 8 | bounds too loose |
| M3 | Open-ended error control | mathematics | Required | near | 0, 1 | correlated tests |
| M4 | Translatability bound | mathematics | Enabling | mid | 5 | too few calibration cases |
| M5 | Alignability to shared knowledge | mathematics | Visionary | mid | 5 | finite-sample behaviour |
| M6 | Calibrated origin odds | mathematics | Enabling | mid | 7 | prior sensitivity |
| M7 | Population inference | mathematics | Enabling | mid | 1, 8 | model misspecification |
| X1 | Prediction ledger practice | practice | Required | near | 1 | few data releases to test against |
| X2 | Rediscovery campaign | practice | Required | near | 1–4 | controls too easy |
| X3 | Blind injection challenges | practice | Enabling | near–mid | 2 | no second party |
| X4 | Null results as coverage | practice | Required | near–mid | 1, 8 | none major |
| X5 | Cross-calibration with the sister project | practice | Enabling | near | 2 | divergence of the two code bases |
| X6 | Follow-up pathway | practice | Enabling | mid | 7 | no access to instruments |

---

## 7. Sequencing

- **First three months (plan Phases 0 and 1):** T12, the red-noise part of T2, the optical gates of T10, T11, a
  first version of M3, X1, the optical part of X2, and the control and natural collections of T9. Everything here
  serves the quasar module and is cheap.
- **Months three to six (Phases 2 and 3):** T1, T3, the rest of T9, M1, X3, X5. This is where the sister project's code
  comes across and gets its fixes.
- **Months six to twelve (Phases 4 and 5):** T4 to T8, M2, M4, M7, X4, and challenge G2.
- **Beyond:** T13, T14, M5, M6, X6, and whatever G1 and G2 reveal.

Every invention that makes an empirical claim is tested on the benchmark collections before it is used on real
candidates, reports null results with their coverage, and is released under the project's licence. Ambition is
measured by what the ledger shows was predicted before it happened.
