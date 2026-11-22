# Runs

Every search that touches real data has a directory here and an entry in [`../ledger/trials.jsonl`](../ledger/trials.jsonl)
that was written, committed and pushed before the search ran.

| Run | What it asked | Outcome |
|---|---|---|
| [`2026-11-22-quasar-pix2-100`](2026-11-22-quasar-pix2-100/) | Does any of 1,210 ZTF quasars in one patch of sky vary with a strict period? | One object promoted, and explained: an RR Lyrae star whose half-day pulsation aliases to 80 days. No quasar promoted; no decoy promoted |
| [`2026-11-22-controls`](2026-11-22-controls/) | What does the same search say about a famous periodic candidate and a famous strange star? | Neither promoted |
| [`2026-11-22-levels`](2026-11-22-levels/) | Does any of the same quasars sit on discrete brightness levels? | 3 real passes against 2 decoy passes: chance level |
| [`2026-11-22-frb-timing`](2026-11-22-frb-timing/) | Do the bursts of two repeating fast radio bursts fall on a clock or a spacing unit? | No: smallest probability 0.10 against a threshold of 0.0025 |

## What a blind run directory holds

| File | What it is |
|---|---|
| `config.json` | The full configuration and its hash, as registered |
| `sealed_commitment.txt` | A hash of the hidden label table, written before the analysis |
| `sealed_proof.json` | The label table and the salt, published after: anyone can recompute the hash |
| `sealed_results.json.gz` | The analysis as it stood before the labels were opened (runs made after the first one) |
| `results.csv` | One line per sealed subject: what it was, its probability, each gate's verdict |
| `results.json.gz` | Every subject's full record |
| `summary.json` | Counts, calibration, completeness and upper limits |
| `dossiers/` | A written record for every promoted object and for the five most significant real ones |

## Checking a run

```python
import json
from stargazer.trust.blind import SealedSet
from stargazer.ledger.trials import TrialsLedger

proof = json.load(open("runs/2026-11-22-quasar-pix2-100/sealed_proof.json"))
assert SealedSet.verify(proof)                       # the label table matches the hash committed beforehand
assert TrialsLedger("ledger/trials.jsonl").verify()  # no ledger entry has been altered or removed
```

The Git history dates each step independently: the registration commit comes before the results commit.
