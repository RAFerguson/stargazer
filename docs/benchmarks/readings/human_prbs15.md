<!-- What the workbench was given: The bare x^15 + x^14 + 1 pseudo-random bit sequence. It was told nothing else. -->

# Reading of human/prbs15

**Content tier 3 of 5: layout recovered.** A layout was recovered: how the signal is organised is known. What it says is not. Nothing in it matched the fixed library of shared knowledge.

Input: BitStream, 20,000 elements.

## Evidence

- (tier 3, shift register, p = 0) Every one of the 20,000 bits follows from the first 15 by one fixed rule (a linear feedback shift register of length 15). Random bits would need a register about 10,000 long.

## Notes

- A shift register that reproduces the whole stream explains every other pattern in it, so no further tests were made. Such sequences are used as test patterns, spreading codes and scramblers; this one carries no content beyond its rule.

## Bits needed to write the stream down

- raw: 20,000
- shift register: 59

## Everything that was tried

| Test | Probability | Found | Note |
|---|---|---|---|
| structure battery | 0.02 | no | null: markov order 0 |
| shift register | 0 | yes | linear complexity 15 of 20000 bits |

## Limits

Tests made: discrete levels and clocks; a battery of structure statistics against a fitted short-memory chain; raster and volume layouts; frame lengths and fixed markers; linear parity checks up to the block length set; linear shift registers; and a fixed library of six anchors. Not seen by any of them: non-linear codes, variable-length frames, compression, encryption, and any convention outside the library.
