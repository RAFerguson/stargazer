<!-- What the workbench was given: A Lincos-style primer of true sums and products in 8-bit tokens. It was told nothing else. -->

# Reading of message/arithmetic_primer

**Content tier 5 of 5: a hypothesised reading.** A specific reading is proposed and has predicted held-out parts of the signal. It is a hypothesis, open to refutation by more of the same signal.

Input: BitStream, 19,200 elements.

## Evidence

- (tier 2, structure battery, p = 0.005) The stream has order that a markov order 3 chain fitted to it does not produce (probability 0.005 over 199 simulated chains; strongest sign: entropy drop).
- (tier 3, raster width, p = 0) Written as rows of 48 bits the stream takes 2,604 fewer bits to describe than as a line with short memory, after paying for the search over widths: neighbouring rows resemble each other. That is a raster 48 wide and 400 deep. (19200 = 2^8 x 3 x 5^2: 52 exact rectangles; width 48 divides it (400 rows).)
- (tier 3, volume height, p = 5.1e-138) The rows in turn repeat in groups of 17: each group resembles the one before it, saving a further 456 bits. That is either a block 48 x 17 x 23 or a longer unit of 816 bits that recurs.
- (tier 3, frame length, p = 0) The stream agrees with itself 16 bits later far more often than chance or short memory allows: it is built of units 16 bits long. 12 positions of every unit hold the same value each time (pattern 11111111000000, starting at bit 8 of the unit): a fixed marker.
- (tier 3, anchor: repeated_block, p = 6e-300) The stream agrees with itself shifted by 48 bits at 85% of positions, where a short-memory model predicts 51%: a unit of that length recurs with changing content.
- (tier 5, anchor: functional_relations, p = 4.7e-142) In records of six 8-bit fields, field 5 equals field 1 plus field 3 whenever field 2 is 252 (195 of 195 records) and field 1 times field 3 whenever field 2 is 253 (205 of 205 records); field 4 is always 254; field 6 is always 255 (from bit 0, most significant bit first). As a hypothesis: the stream states arithmetic, with 252 as '+', 253 as '*'. Worked out from the first half of the stream alone, this reading predicts the result field of 200 of 200 records in the second half.

## Bits needed to write the stream down

- raw: 19,200
- short-memory model: 10,820
- raster of width 48: 8,204

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| structure battery | 0.005 | yes | null: markov order 3 |
| shift register | 1 | no | linear complexity 9600 of 19200 bits |
| raster width | 0 | yes | best width 48, gain 2604 bits |
| volume height | 5.1e-138 | yes | best height 17 |
| frame length | 0 | yes | frame 16, 46.5 sd above short-memory surrogates |
| anchor: repeated_block | 6e-300 | yes | library of 6; 9,599 settings inside this anchor |
| anchor: functional_relations | 4.7e-142 | yes | library of 6; 439,848 settings inside this anchor |
| anchor: run_length_sequence | 1 | no | library of 6; 12 settings inside this anchor |
| anchor: field_counting | 1 | no | library of 6; 810 settings inside this anchor |
| anchor: constants | 1 | no | library of 6; 1,258,368,000 settings inside this anchor |
| anchor: text_code | 1 | no | library of 6; 40 settings inside this anchor |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
