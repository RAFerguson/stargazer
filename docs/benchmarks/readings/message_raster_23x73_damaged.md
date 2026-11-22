<!-- The same raster after damage: 17 random bits in front and 3% of all bits flipped. -->

# Reading of message/raster_23x73, 17 junk bits in front, 3% of bits flipped

**Content tier 3 of 5: layout recovered.** A layout was recovered: how the signal is organised is known. What it says is not. Nothing in it matched the fixed library of shared knowledge.

Input: BitStream, 1,696 elements.

## Evidence

- (tier 3, anchor: repeated_block, p = 7e-05) The stream agrees with itself shifted by 23 bits at 83% of positions, where a short-memory model predicts 77%: a unit of that length recurs with changing content.

## Bits needed to write the stream down

- raw: 1,696
- short-memory model: 865

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| structure battery | 0.12 | no | null: markov order 2 |
| shift register | 1 | no | linear complexity 849 of 1696 bits |
| raster width | 0.00034 | no | best width 23, gain 12 bits |
| frame length | 1 | no | best lag 811 |
| parity checks | 0.0029 | no | none |
| anchor: repeated_block | 7e-05 | yes | library of 6; 847 settings inside this anchor |
| anchor: run_length_sequence | 1 | no | library of 6; 12 settings inside this anchor |
| anchor: field_counting | 1 | no | library of 6; 810 settings inside this anchor |
| anchor: constants | 1 | no | library of 6; 111,155,840 settings inside this anchor |
| anchor: text_code | 1 | no | library of 6; 40 settings inside this anchor |
| anchor: functional_relations | 1 | no | library of 6; 410,952 settings inside this anchor |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
