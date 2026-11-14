# First light: can we reach real data, and does the method hold up?

This is throwaway exploratory code from the drafting stage, kept so the feasibility claims in
[`docs/project-draft.md`](../../docs/project-draft.md) can be rerun. It is not the start of the
Stargazer package.

It answers two questions:

1. **Which data sources can this machine reach?** `check_sources.py` probes each host in
   `sources.txt` once and reports reachable, blocked, or unreachable.
2. **Can we pull a real quasar light curve and test it for "pulsing" correctly?** The other three
   scripts fetch PG 1302−102 from the public ZTF DR24 cloud copy, test its strongest periodogram
   peak against a red-noise null, and draw the figure used in the draft.

## Run it

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt

python check_sources.py                                              # standard library only
python fetch_ztf_lightcurve.py pg1302 196.38754 -10.555389           # about 10 s, no account needed
python rednoise_test.py data/pg1302_ztf_dr24.json --claimed-period 1884   # about 20 s
python plot_first_light.py data/pg1302_rednoise.json ../../docs/figures/first-light-pg1302
```

## What it found on 14 November 2026

| Check | Result |
|---|---|
| Hosts reachable from the cloud development environment | 10 of 46: the five public AWS buckets, Google Cloud Storage, GitHub, PyPI and conda |
| Hosts blocked by the environment's network policy | 36 of 46: every archive website and API |
| Time to fetch PG 1302−102's light curves | 7.7 s, reading 2 of 6 row groups from one 409 MB Parquet file |
| Light curves returned | g: 254 epochs, r: 306 epochs, i: 69 epochs, spanning April 2018 to July 2025 |
| Strongest r-band periodogram peak | 2,638 days, which is the full length of the data (one cycle) |
| False-alarm probability assuming white noise | about 10⁻⁸⁰ |
| False-alarm probability against fitted red noise | 0.067 (r band), 0.087 (g band) |
| Red-noise false-alarm probability at the published 1,884-day period | 0.13 in both bands |

The last four rows are the point. A textbook periodogram test calls this quasar a
one-in-10⁸⁰ detection. Calibrated against the red noise that every quasar shows, the same peak
turns up by chance about once in fifteen tries. Seven years of ZTF data covers only 1.4 cycles
of the published period, so this data alone can neither confirm nor refute it.

## Limits of this test

- The red-noise model is a single damped random walk with log-uniform priors on amplitude
  (0.02 to 1 mag) and timescale (5 to 20,000 days). The timescale is poorly constrained when it
  approaches the length of the data, as it does here.
- Only one survey is used. A real test needs a longer baseline stitched from several surveys.
- Photometry is taken as delivered, keeping epochs with `catflags == 0`.

## Data credit

`data/pg1302_ztf_dr24.json` is an extract of the ZTF DR24 light-curve catalogue served by IRSA
(doi:[10.26131/IRSA598](https://www.ipac.caltech.edu/doi/10.26131/IRSA598)). Follow the ZTF and
IRSA acknowledgement guidelines if you reuse it.
