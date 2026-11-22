<!-- What the workbench was given: Telemetry frames with random payloads: the envelope without readable content. It was told nothing else. -->

# Reading of human/telemetry_random

**Content tier 3 of 5: layout recovered.** A layout was recovered: how the signal is organised is known. What it says is not. Nothing in it matched the fixed library of shared knowledge.

Input: BitStream, 19,840 elements.

## Evidence

- (tier 2, structure battery, p = 0.005) The stream has order that a markov order 2 chain fitted to it does not produce (probability 0.005 over 199 simulated chains; strongest sign: entropy drop).
- (tier 3, raster width, p = 1.2e-64) Written as rows of 496 bits the stream takes 212 fewer bits to describe than as a line with short memory, after paying for the search over widths: neighbouring rows resemble each other. That is a raster 496 wide and 40 deep. (19840 = 2^7 x 5 x 31: 30 exact rectangles; width 496 divides it (40 rows).)
- (tier 3, volume height, p = 4.8e-11) The rows in turn repeat in groups of 2: each group resembles the one before it, saving a further 34 bits. That is either a block 496 x 2 x 20 or a longer unit of 992 bits that recurs.
- (tier 3, frame length, p = 5.6e-240) The stream agrees with itself 496 bits later far more often than chance or short memory allows: it is built of units 496 bits long. 43 positions of every unit hold the same value each time (pattern 0001101011001111111111000001110100000000000, starting at bit 0 of the unit): a fixed marker.
- (tier 3, anchor: repeated_block, p = 1.3e-34) The stream agrees with itself shifted by 496 bits at 55% of positions, where a short-memory model predicts 50%: a unit of that length recurs with changing content.

## Bits needed to write the stream down

- raw: 19,840
- short-memory model: 19,827
- raster of width 496: 19,603

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| structure battery | 0.005 | yes | null: markov order 2 |
| shift register | 1 | no | linear complexity 9920 of 19840 bits |
| raster width | 1.2e-64 | yes | best width 496, gain 212 bits |
| volume height | 4.8e-11 | yes | best height 2 |
| frame length | 5.6e-240 | yes | frame 496, 13.3 sd above short-memory surrogates |
| anchor: repeated_block | 1.3e-34 | yes | library of 6; 9,919 settings inside this anchor |
| anchor: run_length_sequence | 1 | no | library of 6; 12 settings inside this anchor |
| anchor: field_counting | 1 | no | library of 6; 810 settings inside this anchor |
| anchor: constants | 1 | no | library of 6; 1,300,313,600 settings inside this anchor |
| anchor: text_code | 1 | no | library of 6; 40 settings inside this anchor |
| anchor: functional_relations | 1 | no | library of 6; 362,880 settings inside this anchor |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
