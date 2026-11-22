<!-- What the workbench was given: English-like text XORed with a SHA-256 counter-mode keystream; must read as noise. It was told nothing else. -->

# Reading of control/encrypted_text

**Content tier 0 of 5: no content analysis.** No structure was found. By these tests the signal cannot be told from its null model. That is not proof that it carries nothing: well-compressed or encrypted content looks exactly like this.

Input: BitStream, 24,000 elements.

## Bits needed to write the stream down

- raw: 24,000
- short-memory model: 24,011

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| structure battery | 0.59 | no | null: markov order 0 |
| shift register | 1 | no | linear complexity 12000 of 24000 bits |
| raster width | 1 | no | best width 2, gain -27 bits |
| frame length | 0.38 | no | best lag 3053 |
| parity checks | 1 | no | none |
| anchor: text_code | 0.16 | no | library of 6; 40 settings inside this anchor |
| anchor: repeated_block | 1 | no | library of 6; 11,999 settings inside this anchor |
| anchor: run_length_sequence | 1 | no | library of 6; 12 settings inside this anchor |
| anchor: field_counting | 1 | no | library of 6; 810 settings inside this anchor |
| anchor: constants | 1 | no | library of 6; 1,572,960,000 settings inside this anchor |
| anchor: functional_relations | 1 | no | library of 6; 362,880 settings inside this anchor |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
