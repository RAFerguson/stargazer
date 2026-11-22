<!-- What the workbench was given: The 62 primes below 300 in unary. It was told nothing else. -->

# Reading of message/unary_primes_300

**Content tier 4 of 5: content anchored to shared knowledge.** Part of the content matches shared knowledge from the fixed library, with the probability of a chance match corrected for everything in that library.

Input: BitStream, 8,336 elements.

## Evidence

- (tier 4, anchor: run_length_sequence, p = 5.4e-110) 62 consecutive runs of ones, starting at bit 0, have lengths 2, 3, 5, 7, 11, ..., 293, which are consecutive prime numbers.

## Bits needed to write the stream down

- raw: 8,336
- short-memory model: 530

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| structure battery | 0.02 | no | null: markov order 0 |
| shift register | 1 | no | linear complexity 4169 of 8336 bits |
| raster width | 0.019 | no | best width 310, gain 6 bits |
| frame length | 4.8e-49 | no | frame 8, 0.1 sd above short-memory surrogates |
| parity checks | 1.7e-06 | no | (10,5) at offset 2, 4.8 sd above short-memory surrogates |
| anchor: run_length_sequence | 5.4e-110 | yes | library of 6; 12 settings inside this anchor |
| anchor: field_counting | 1 | no | library of 6; 810 settings inside this anchor |
| anchor: constants | 1 | no | library of 6; 546,341,440 settings inside this anchor |
| anchor: text_code | 1 | no | library of 6; 40 settings inside this anchor |
| anchor: functional_relations | 1 | no | library of 6; 38,616 settings inside this anchor |
| anchor: repeated_block | 1 | no | library of 6; 4,167 settings inside this anchor |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
