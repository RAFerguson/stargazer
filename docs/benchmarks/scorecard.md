| Signal | Kind | Justified tier | Tier given | Verdict | Evidence found |
|---|---|---|---|---|---|
| control/biased_bits_p30 | BitStream | 0 | 0 | reached | none |
| control/encrypted_text | BitStream | 0 | 0 | reached | none |
| control/fair_bits | BitStream | 0 | 0 | reached | none |
| control/iid_symbols_16 | SymbolSequence | 0 | 0 | reached | none |
| control/iid_symbols_4 | SymbolSequence | 0 | 0 | reached | none |
| control/markov1_bits | BitStream | 0 | 0 | reached | none |
| control/poisson_events | EventList | 0 | 0 | reached | none |
| control/white_noise_series | Series | 0 | 0 | reached | none |
| control/zlib_text | BitStream | 0 | 0 | reached | none |
| human/ascii_8n1 | BitStream | 4 | 4 | reached | structure battery (tier 2); frame length (tier 3); parity checks (tier 3); anchor: text_code (tier 4); anchor: repeated_block (tier 3) |
| human/ascii_manchester | BitStream | 3 | 3 | reached | structure battery (tier 2); frame length (tier 3); parity checks (tier 3); anchor: repeated_block (tier 3); anchor: text_code (tier 3) |
| human/hamming74_text | BitStream | 3 | 3 | reached | structure battery (tier 2); frame length (tier 3); parity checks (tier 3); anchor: repeated_block (tier 3); anchor: text_code (tier 3) |
| human/lfsr_scrambled_text | BitStream | 3 | 0 | short | none |
| human/ook_pulses | EventList | 4 | 4 | reached | clock (tier 1); structure battery (tier 2); frame length (tier 3); parity checks (tier 3); anchor: text_code (tier 4); anchor: repeated_block (tier 3) |
| human/parity87_text | BitStream | 3 | 3 | reached | structure battery (tier 2); frame length (tier 3); parity checks (tier 3); anchor: text_code (tier 3); anchor: repeated_block (tier 3) |
| human/prbs15 | BitStream | 3 | 3 | reached | shift register (tier 3) |
| human/telemetry_random | BitStream | 3 | 3 | reached | structure battery (tier 2); raster width (tier 3); volume height (tier 3); frame length (tier 3); anchor: repeated_block (tier 3) |
| human/telemetry_text | BitStream | 4 | 4 | reached | structure battery (tier 2); frame length (tier 3); anchor: repeated_block (tier 3); anchor: text_code (tier 4) |
| message/arithmetic_primer | BitStream | 5 | 5 | reached | structure battery (tier 2); raster width (tier 3); volume height (tier 3); frame length (tier 3); anchor: repeated_block (tier 3); anchor: functional_relations (tier 5) |
| message/binary_count_255 | BitStream | 4 | 4 | reached | structure battery (tier 2); raster width (tier 3); volume height (tier 3); frame length (tier 3); anchor: field_counting (tier 4); anchor: repeated_block (tier 3) |
| message/raster3d_11x13x17 | BitStream | 3 | 3 | reached | raster width (tier 3); volume height (tier 3); frame length (tier 3); anchor: repeated_block (tier 3) |
| message/raster_127x127 | BitStream | 3 | 3 | reached | structure battery (tier 2); raster width (tier 3); frame length (tier 3); anchor: repeated_block (tier 3) |
| message/raster_23x73 | BitStream | 3 | 3 | reached | raster width (tier 3); anchor: repeated_block (tier 3) |
| message/raster_23x73_framed | BitStream | 3 | 3 | reached | raster width (tier 3); anchor: repeated_block (tier 3) |
| message/unary_count_1_40 | BitStream | 4 | 4 | reached | anchor: run_length_sequence (tier 4) |
| message/unary_primes_300 | BitStream | 4 | 4 | reached | anchor: run_length_sequence (tier 4) |
| natural/drw_plus_weak_sine | Series | 0 | 0 | reached | none |
| natural/drw_quasar | Series | 0 | 0 | reached | none |
| natural/eclipsing_binary | Series | 1 | 1 | reached | discrete levels (tier 1) |
| natural/pulsar_nulling | EventList | 1 | 1 | reached | clock (tier 1) |
| natural/pulsar_nulling_bits | BitStream | 0 | 0 | reached | none |
| natural/repeating_burster | EventList | 0 | 0 | reached | none |
| natural/variable_star_p180 | Series | 0 | 0 | reached | none |

| Collection | Items | Reached | Short | Overclaim |
|---|---|---|---|---|
| control | 9 | 9 | 0 | 0 |
| human | 9 | 8 | 1 | 0 |
| message | 8 | 8 | 0 | 0 |
| natural | 7 | 7 | 0 | 0 |
