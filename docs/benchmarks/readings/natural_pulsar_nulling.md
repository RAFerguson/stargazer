<!-- What the workbench was given: A 0.714 s pulsar with timing jitter, lognormal amplitudes and Markov nulling. It was told nothing else. -->

# Reading of natural/pulsar_nulling

**Content tier 1 of 5: discrete symbols recovered.** The signal takes discrete states, so it can be written down as a sequence of symbols. The sequence shows no order beyond what a short-memory or repeating natural process produces.

Input: EventList, 3,871 elements.

## Evidence

- (tier 1, clock, p = 0.005) The events fall on the ticks of a clock of period 0.714 s, scattering by 0.3% of a period; 77% of ticks carry an event. Each tick is therefore a slot that is either full or empty.

## Bits needed to write the stream down

- raw: 5,000
- short-memory model: 2,715

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| clock | 0.005 | yes | period 0.714 |
| structure battery | 0.67 | no | null: markov order 1 |
| shift register | 1 | no | linear complexity 2501 of 5000 bits |
| raster width | 1 | no | best width 2, gain -15 bits |
| frame length | 0.017 | no | best lag 1014 |
| parity checks | 2.8e-19 | no | (13,1) at offset 3, 1.6 sd above short-memory surrogates |
| anchor: run_length_sequence | 1 | no | library of 6; 12 settings inside this anchor |
| anchor: field_counting | 1 | no | library of 6; 810 settings inside this anchor |
| anchor: constants | 1 | no | library of 6; 327,700,000 settings inside this anchor |
| anchor: text_code | 1 | no | library of 6; 40 settings inside this anchor |
| anchor: functional_relations | 1 | no | library of 6; 362,880 settings inside this anchor |
| anchor: repeated_block | 1 | no | library of 6; 2,499 settings inside this anchor |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
