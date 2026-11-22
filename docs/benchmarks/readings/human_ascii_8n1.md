<!-- What the workbench was given: ASCII text sent as 8N1 asynchronous serial frames. It was told nothing else. -->

# Reading of human/ascii_8n1

**Content tier 4 of 5: content anchored to shared knowledge.** Part of the content matches shared knowledge from the fixed library, with the probability of a chance match corrected for everything in that library.

Input: BitStream, 6,000 elements.

## Evidence

- (tier 2, structure battery, p = 0.005) The stream has order that a markov order 3 chain fitted to it does not produce (probability 0.005 over 199 simulated chains; strongest sign: entropy drop).
- (tier 3, frame length, p = 0) The stream agrees with itself 10 bits later far more often than chance or short memory allows: it is built of units 10 bits long. 7 positions of every unit hold the same value each time (pattern 01101000, starting at bit 5 of the unit): a fixed marker.
- (tier 3, parity checks, p = 3.6e-138) Cut into blocks of 5 bits from bit 0, the stream obeys 4 independent parity rules in nearly every block (the best in 89.6% of blocks). Each block therefore carries at most 1 free bits out of 5: fixed positions, or an error-detecting code.
- (tier 3, anchor: repeated_block, p = 6e-300) The stream agrees with itself shifted by 80 bits at 75% of positions, where a short-memory model predicts 51%: a unit of that length recurs with changing content.
- (tier 4, anchor: text_code, p = 6e-300) Read as asynchronous serial frames (a 0 start bit, eight data bits least significant first, a 1 stop bit) from bit 0, 600 of 600 frames are well formed and hold printable ASCII, against about 16% at an arbitrary alignment or by chance; the commonest character is the space. It begins: "About the word to is into the of was it what. And at in is and small cal".

## Bits needed to write the stream down

- raw: 6,000
- short-memory model: 3,271

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| structure battery | 0.005 | yes | null: markov order 3 |
| shift register | 1 | no | linear complexity 3001 of 6000 bits |
| raster width | 1 | no | best width 10, gain -532 bits |
| frame length | 0 | yes | frame 10, 30.7 sd above short-memory surrogates |
| parity checks | 3.6e-138 | yes | (5,1) at offset 0, 21.8 sd above short-memory surrogates |
| anchor: text_code | 6e-300 | yes | library of 6; 40 settings inside this anchor |
| anchor: repeated_block | 6e-300 | yes | library of 6; 2,999 settings inside this anchor |
| anchor: run_length_sequence | 1 | no | library of 6; 12 settings inside this anchor |
| anchor: field_counting | 1 | no | library of 6; 810 settings inside this anchor |
| anchor: constants | 1 | no | library of 6; 393,240,000 settings inside this anchor |
| anchor: functional_relations | 1 | no | library of 6; 360,000 settings inside this anchor |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
