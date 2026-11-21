# Radio technosignature (radio SETI) data archives, access paths, formats and software (state as of 21 November 2026)

Tag legend: **[P]** read today (21 Nov 2026) in the primary source (page, paper, repository file, API response, or my own test run). **[S]** web-search summary only. **[R]** recalled from background knowledge, not re-checked today.

Method notes. Pages were fetched with curl through the sandbox proxy. Repository files come from shallow clones that an earlier, interrupted researcher made on 14 Nov 2026; I re-read them today and treat them as [P] for their contents as of that date. PyPI metadata was queried live today. Benchmarks were run on a 4-vCPU, 15 GB RAM, CPU-only VM. Hosts the proxy blocks (github.com, kaggle.com, bldata.berkeley.edu, blpd*.ssl.berkeley.edu, zooniverse.org, CASDA, CIRADA, NRAO, FAST data centre) are listed as gaps, not verified.

---

## 1. Breakthrough Listen (BL) Open Data Archive: what is downloadable, how, in what form, and programme status

### Takeaway
The BL Open Data Archive is online in November 2026. Its search front end (https://breakthroughinitiatives.org/opendatasearch) answers GET queries, and paging through results counts about 859,000 catalogued files from GBT, Parkes and APF. It has no documented JSON API. The "beta" JSON API at seti.berkeley.edu/opendata exists in source code, but the HTTPS server now redirects it to plain HTTP, so it could not be tested from here. MeerKAT/BLUSE, COSMIC/VLA and ATA data are not in the archive search (MeerKAT is not offered as a project); they are available only on request or through ad-hoc pages. Download links are spread across three host families, one of them HTTP-only. Breakthrough Listen is headquartered at Oxford, its original 10-year, $100M commitment (July 2015) has run out, and I found no public announcement of renewed funding.

### Cited findings

**Archive front end and query interface**
- The archive search page is live and lists projects "Breakthrough Listen at Green Bank", "Automated Planet Finder" and "Breakthrough Listen at Parkes". File types offered are "baseband data", "filterbank", "HDF5" and "data" (= APF FITS). Its text warns that "file sizes can be several gigabytes" and points to "a beta interface with additional search options … at seti.berkeley.edu/opendata" [P] — [BI Open Data Archive](https://breakthroughinitiatives.org/opendatasearch)
- The search form is a plain HTTP GET with parameters `project` (all|GBT|APF|Parkes), `file_type` (all|baseband+data|filterbank|HDF5|data), `ra`, `ra_range`, `decl`, `decl_range` (degrees), `mjd`, `mjd_range`, `freq`, `freq_range` (MHz), `target`, `perPage`, `page`, `sort` (utc|mjd|project|target|ra|decl|centerfreq|filetype|size) and `order` (asc|desc). Results come back as an HTML table with columns UTC, MJD, Project, Target, RA, Decl, Center Freq, File Type, Size (bytes), Md5sum, Download URL. I tested `?project=GBT&file_type=HDF5&perPage=10&search=Search` and `?target=VOYAGER1&perPage=100` today and both returned tables [P] — [BI Open Data Archive](https://breakthroughinitiatives.org/opendatasearch)
- `perPage=100` returns 100 rows. Unsupported values such as `perPage=1000` silently fall back to 20 rows. Pages past the end return an empty table, and the pagination widget only shows a window of about 50 pages [P] (my test) — [BI Open Data Archive](https://breakthroughinitiatives.org/opendatasearch)
- Inference for tool design: the only working machine interface today is HTML scraping of this GET form. Rows include MD5 sums, so downloads can be checked for integrity.
- File counts from paging through the archive on 21 Nov 2026 (exact last non-empty page found by binary search; 20-row pages) [P] (my scrape; script at `.../radio_seti/scripts/bl_count2.py`):

| Project | File type | Files | Example sizes (first page) | Download host(s) seen |
|---|---|---|---|---|
| GBT | baseband (GUPPI raw) | 42,297 | 1.6–16.9 GB per raw file | storage.googleapis.com/gbt_guppi (HTTPS), blpd1.ssl.berkeley.edu (HTTP) |
| GBT | filterbank (.fil) | 10,610 | 14.7 MB – 76.5 GB | blpd0/1/5.ssl.berkeley.edu (HTTP only) |
| GBT | HDF5 (.h5) | 78,132 | 0.24–15.6 GB | bldata.berkeley.edu (HTTP and HTTPS) |
| APF | data (FITS) | 4,767 | 19.2 MB | bldata.berkeley.edu/apf/ |
| Parkes | baseband | 8,942 | 7.6–11.8 GB | blpd8/12.ssl.berkeley.edu (HTTP) |
| Parkes | filterbank | 366,568 | 48 MB – 1.15 GB | blpd2/4.ssl.berkeley.edu (HTTP) |
| Parkes | HDF5 | 347,810 | 4–15 GB | blpd1/8.ssl.berkeley.edu (HTTP) |
| **Total** | | **≈859,126** | | |

- The newest entries in the listing: GBT HDF5 up to 2022-11-19 (e.g. `http://bldata.berkeley.edu/pipeline/AGBT22B_999_25/.../blc16_guppi_59902_23686_NGC3226_0056.rawspec.0000.h5`, 3.2 GB). GBT raw/filterbank up to 2020-09-11 (Kepler-160 DIAG data). Parkes up to 2021-05-03. APF up to 2019-03-25 [P] — [BI Open Data Archive](https://breakthroughinitiatives.org/opendatasearch)
- Inference: the public search index appears not to have been extended with GBT data after late 2022, or with Parkes data after mid-2021. Newer data (MeerKAT, COSMIC, 3I/ATLAS) are published, if at all, through ad-hoc directories (see below).
- BLDR 1.0 (2019) counted 72,925 GBT files (40,375 raw voltage, 1,456 filterbank, 31,094 HDF5) and 6,086 APF FITS. The PDF-extracted Parkes figures (7,358 raw, 8,165 filterbank, total 15,523) come out jumbled, so check them against the original table [P] — [Lebofsky et al. 2019, Table 8](https://arxiv.org/abs/1906.07391)
- Inference: the archive has grown about 10× in file count since BLDR 1.0, mostly from Parkes multibeam filterbank/HDF5. The APF count today (4,767) is lower than the 6,086 FITS files in 2019, which I cannot explain.
- The "beta" interface seti.berkeley.edu/opendata (Flask + MySQL) defines a JSON API: `/api/list-targets` (optional `simbad`), `/api/list-telescopes`, `/api/list-file-types`, and `/api/query-files` with parameters `target` (required; prefix `!` = exact, `/` = regex, otherwise SQL LIKE), `telescopes` (comma-separated), `file-types` (comma-separated), `pos-ra`, `pos-dec`, `pos-rad` (deg), `time-start`, `time-end` (MJD), `freq-start`, `freq-end` (MHz) and `limit`. It returns `{"result": "success", "data": [{target, telescope, utc, mjd, ra, decl, center_freq, file_type, size, md5sum, url}]}` [P] — source code in the repository sxyu/bl-opendata, `blopendata/core.py` (GitHub; read from an earlier clone)
- Live status of that API: today `https://seti.berkeley.edu/opendata/` and every `/opendata/api/*` path return HTTP 302 to `http://seti.berkeley.edu/...` (Apache/2.4.18 Ubuntu). The sandbox cannot follow the plain-HTTP redirect, so I could not confirm the API still answers [P] — [seti.berkeley.edu/opendata](https://seti.berkeley.edu/opendata)
- Download-host reachability today:
  - `https://storage.googleapis.com/gbt_guppi/blc3_guppi_57386_VOYAGER1_0004.0000.raw` answers HTTP 200, content-length 16,929,030,144, last-modified 2016-01-05. Listing the `gbt_guppi` and `gbt_fil` buckets returns 403.
  - `https://storage.googleapis.com/gbt_fil/voyager_f1032192_t300_v2.fil` downloaded in full (528,482,692 bytes in 19.4 s, about 27 MB/s).
  - bldata.berkeley.edu is blocked by the sandbox proxy (403), and blpd*.ssl.berkeley.edu does not resolve here.

  [P] (my tests) — [Google Cloud bucket object](https://storage.googleapis.com/gbt_guppi/blc3_guppi_57386_VOYAGER1_0004.0000.raw)
- Inference: most archive links outside the two GCS buckets use plain-HTTP Berkeley hosts (blpd*.ssl.berkeley.edu). A tool should expect mixed HTTP/HTTPS, no checksums from the server itself (use the MD5 in the listing), and large single files.

**File formats and standard data products**
- GBT reduced products are three "spliced" filterbank files per observation [P] — [UCBerkeleySETI/breakthrough GBT/README.md](https://raw.githubusercontent.com/UCBerkeleySETI/breakthrough/master/GBT/README.md):
  - high frequency resolution (HSR; about 3 Hz, about 18 s), file suffix `.0000`
  - high time resolution (HTR; about 366 kHz, about 349 µs), suffix `.0001`
  - medium resolution (MR; about 3 kHz, about 1 s), suffix `.0002`
- Exact channelisation, Revision 1A [P] — [Lebofsky et al. 2019, Table 4](https://arxiv.org/abs/1906.07391):
  - HSR: FFT length 1,033,216, 51 integrations, 2.84 Hz, 17.98 s
  - HTR: FFT length 8, 128 integrations, 366.21 kHz, 349.53 µs
  - MR: FFT length 1024, 3027 integrations, 2.86 kHz, 1.06 s
- Revision 1B uses 999,424 channels per coarse channel (2.93 Hz, 17.40 s). The newer rawspec products use 1,048,576 fine channels per 2.93 MHz coarse channel, about 2.79 Hz and 18.25 s [P] — [Lebofsky et al. 2019](https://arxiv.org/abs/1906.07391); [rawspec README](https://raw.githubusercontent.com/UCBerkeleySETI/rawspec/master/README.md)
- Formats [P] — [Lebofsky et al. 2019](https://arxiv.org/abs/1906.07391); [rawspec README](https://raw.githubusercontent.com/UCBerkeleySETI/rawspec/master/README.md):
  - SIGPROC filterbank (.fil): header plus a float32 (time, IF, channel) array.
  - FBH5 (.h5): HDF5 with dataset `data` of shape (time, feed_id, frequency), header attributes, and bitshuffle/LZ4 compression, which needs the HDF5 bitshuffle plugin (hdf5plugin).
  - GUPPI raw (.raw): repeated header units plus 2- or 8-bit complex voltage blocks, documented in `breakthrough/doc/RAW-File-Format.md`.
  - BL converts filterbank to HDF5 before public release.
- GBT Gen-1 recorder: 750 MB/s raw per compute node, 64 nodes, 187.5 MHz per node, 2.92 MHz coarse channels [P] — [GBT README](https://raw.githubusercontent.com/UCBerkeleySETI/breakthrough/master/GBT/README.md)
- APF products: raw 2D echelle image (18.3 MB each) and extracted 1D spectrum (1.4 MB each); 789 targets in BLDR 1.0 [P] — [Lebofsky et al. 2019](https://arxiv.org/abs/1906.07391)
- Sample file I downloaded today, `voyager_f1032192_t300_v2.fil` (528,482,692 bytes), has header nchans = 66,060,288 (63 coarse channels × 1,048,576), foff = −2.838 Hz, fch1 = 8588.965 MHz, tsamp = 105.70 s, and only **2** time integrations; source VOYAGER1, telescope_id 6 (GBT), MJD 57386.86 [P] (blimpy header read).

**ON/OFF (ABACAD) cadence**
- The GBT strategy is 5 min on primary target A, 5 min on a nearby B, A, C, A, D (an "ABACAD" cadence). A signal is interesting only if it appears in all three A scans and none of B/C/D. The primary list is about 1,700 stars [P] — [GBT README](https://raw.githubusercontent.com/UCBerkeleySETI/breakthrough/master/GBT/README.md)
- Archive target names use suffixes such as `_OFF`, `_R`/`_S` (Parkes reference and source) and `_B1…B13` (Parkes multibeam beams). Targets ending `_OFF` are pointed "typically less than 2 degrees Dec" away [P] — [BI Open Data Archive](https://breakthroughinitiatives.org/opendatasearch)
- Parkes BLC1/Proxima data used an ABABAB (6-pointing) and longer cadences [P] — [Smith et al. 2021](https://arxiv.org/abs/2111.08007)

**Volumes**
- BLDR 1.0 is about 1 PB of public raw and reduced data, including the 400 TB FRB 121102 raw dataset [P] — [Lebofsky et al. 2019](https://arxiv.org/abs/1906.07391)
- The Price et al. 2020 data are about 219 TB (142 TB GBT, 77 TB Parkes) [P] — [Price et al. 2020](https://arxiv.org/abs/1906.07750)
- The L-band 2017 survey: 692 stars, 400 h, 4,798 observations, 8.0 PB raw, 180 TB reduced filterbank, 28,858,710 hits [P] — [lband2017 statistics](https://seti.berkeley.edu/lband2017/statistics.html)
- BL-operated storage on 25 March 2026: 41.8 PB capacity, 31.4 PB used. By system (capacity/used, PB): GBT 9.7/5.8; Parkes 3.5/3.1; SRT 0.2/0.2; Berkeley data centre 16.0/15.0; COSMIC 3.6/3.1; MeerKAT 7.0/2.7; ATA 1.8/1.5. The main archive is "hosted 'off prem' at the Berkeley data center" [P] — [MacMahon & … 2026, "The First Decade of Breakthrough Listen", IAC-26](https://arxiv.org/abs/2610.07482)
- BL policy: preserve all reduced data plus a curated raw-voltage subset (about 1 PB/yr), "with both types of data publicly available" [P] — [MacMahon 2026](https://arxiv.org/abs/2610.07482)
- Inference: only a small fraction of the 31 PB is indexed in the public search (the indexed files are a few PB at most). The exact public volume is not published.

**Other BL public pages and datasets (checked today)**
- L-band 2017 (Enriquez et al.): the page is live (last-modified 2025-04-01) and lists top-11 events with HDF5 zips. The hit CSVs ("1.1 million 3Hz hits from the 11 top ranking stars (~420 MB)" and "All 28.9 million 3Hz hits from all 692 stars (~10 GB)") are licensed CC BY 4.0. The page text also says "12 million hits", which contradicts the link label. The files live at `http://blpd0.ssl.berkeley.edu/lband2017/AAA_candidates.v4_1492476400.csv` and `.../All_hits_turbo_seti.csv`, which are HTTP and unverifiable from here [P] — [lband2017 downloads](https://seti.berkeley.edu/lband2017/downloads.html)
- 1327-star DR1 (Price et al. 2020): the event CSVs are live. `lband2019_events.csv` is 908,918 B with about 5,480 event rows; `sband2019_events.csv` is 357,347 B; `pks10cm2019_events.csv` is 10,754 B. Columns are Source, DriftRateMax, Nevent, FreqMid, DriftBW, DriftRates, Freqs, FileID, SNR [P] — [listen2019 open data](https://seti.berkeley.edu/listen2019/opendata.html)
- Zenodo record 7600403 (CC BY 4.0) holds the turboSETI `.dat` outputs for the Price et al. S-band data (`sband2019.tar.gz`, 2.67 GB) and `static.tar.gz` (1.62 GB: event database as HDF5 plus event-viewer images) [P] — [Zenodo 7600403](https://zenodo.org/records/7600403)
- BLC1 / Proxima Centauri: datasets for Smith 2021, Sheikh 2021a and 2021b are linked via `http://seti.berkeley.edu/opendata/?onLoad=1&paperName=Smith21a` (and `Sheikh21a`, `Sheikh21b`) [P] — [seti.berkeley.edu/blc1](https://seti.berkeley.edu/blc1/)
- 3I/ATLAS (2025): data directories at `https://bldata.berkeley.edu/ATLAS/` (`GB_ATLAS/`, `ATA_ATLAS`, `MeerKAT/`). The paper says "All data collected by BL will be publicly available … can be downloaded from the Breakthrough Listen Data Portal" [P] — [seti.berkeley.edu/atlas](https://seti.berkeley.edu/atlas); [Jacobson-Bell et al. 2025](https://arxiv.org/abs/2512.19763). A search summary says raw voltages were still "in the process of being transferred from Green Bank" [S].

**MeerKAT / BLUSE**
- BLUSE has run autonomously since mid-2022. It forms 64 coherent beams plus 1 incoherent beam at about 1 Hz resolution and searches with seticore (Taylor tree, ±10 Hz/s, initial SNR threshold 6). Beams processed: about 73,500 (mid-2022 to mid-2023), 386,200 (2023–24) and 436,700 (2024–25); about 1.2 million viable pointings on about 360,000 unique objects, now about 29,000 unique objects per month. The target list is about 32 million Gaia-derived stars [P] — [Czech et al. 2026](https://arxiv.org/abs/2607.23651)
- BLUSE outputs are hit files plus "stamp" files (upchannelised voltage cutouts around detections, HDF5). Data availability: "shared on request to the corresponding author." Science results are "forthcoming" [P] — [Czech et al. 2026](https://arxiv.org/abs/2607.23651)
- MeerKAT is not a project option in the BL archive search [P] — [BI Open Data Archive](https://breakthroughinitiatives.org/opendatasearch)

**Programme status (2025–2026)**
- Oxford became BL's international headquarters (announced 18 Oct 2023), with PI Andrew Siemion [P] — [BI news 48](https://breakthroughinitiatives.org/news/48)
- The 2026 review restates "a 10-year, $100M commitment" from July 2015. It lists GBT (2016–present), Parkes (2016–present), APF (2017–present), MeerKAT (2018–present), SRT (2019–present), LOFAR international stations (2019–present), VLA/COSMIC (2023–present) and ATA (2023–present). It describes the "fully commensal era" since 2023 [P] — [MacMahon 2026, Table 2](https://arxiv.org/abs/2610.07482)
- New in 2024–2025: SRT partnership (Oct 2024); NVIDIA Holoscan real-time pilot at the ATA (Dec 2024); MeerKAT archive-mining (RATT/TRON, Feb 2025); a Westerbork all-sky PAF monitor with ASTRON and Manchester (May 2025); an ATA real-time ML FRB system (Nov 2025, the latest item on the BI news page today) [P] — [BI news](https://breakthroughinitiatives.org/news)
- Funding: a search summary reports a March 2026 article saying "the original funding cycle ends in 2026" with infrastructure expected to keep operating. I found no announcement of renewed or successor funding [S] — [newspaceeconomy.ca, Oct 2026](https://newspaceeconomy.ca/2026/10/09/what-does-breakthrough-listens-first-decade-show-about-scientific-data-infrastructure/)
- **Flag:** whether GBT and Parkes BL observing continues routinely in 2026 is unconfirmed. The 2026 review says "present" for both. The latest GBT BL observation I could document is 3I/ATLAS on 18 Dec 2025 [P] — [Jacobson-Bell 2025](https://arxiv.org/abs/2512.19763)
- The seti.berkeley.edu home page was last modified 2024-05-16 and still says "Green Bank, Parkes, and Automated Planet Finder" [P] — [seti.berkeley.edu](https://seti.berkeley.edu/)

### Inferences
- An independent tool can reliably get (a) GUPPI raw files from the HTTPS GCS bucket, (b) small hit/event tables from seti.berkeley.edu and Zenodo, and (c) any listed file from Berkeley hosts over plain HTTP. HDF5 products need hdf5plugin for bitshuffle.
- A downloader should handle HTTP-only links, files of 10 GB or more, MD5 checks, and HTML-table scraping. The JSON beta API may work from the open internet but should be treated as unstable.
- The ABACAD cadence is not stored as a single object. It has to be rebuilt from target names (A plus `_OFF` or neighbouring HIP targets) and MJD ordering.

### Gaps
- Could not test the JSON beta API (HTTP redirect) or confirm bldata.berkeley.edu / blpd*.ssl.berkeley.edu downloads (blocked here).
- Total public byte volume is not published. My counts give file numbers only, plus first-page sizes.
- No official 2026 statement on post-2025 BL funding, or on which telescopes still observe for BL every week.
- "Breakthrough Listen Data Portal" (Jacobson-Bell 2025) may be a new name or a new portal. Its URL was not found.

---

## 2. Other public radio SETI data (ATA, COSMIC/VLA, FAST, LOFAR, MWA, UCLA SETI, SETI@home, Big Ear/Wow!, Phoenix)

### Takeaway
Outside BL's GBT/Parkes archive, very little raw radio-SETI data is openly downloadable. Most new programmes (COSMIC, BLUSE, ATA-ASP) release software and papers but provide data only on request. The best non-BL downloads are UCLA SETI's candidate-detection tables on Zenodo (CC0) and the Arecibo Wow! project's archival Ohio SETI material. FAST data follow the observatory's 12-month proprietary policy. SETI@home released code and a web interface, not a bulk candidate database.

### Cited findings
- **COSMIC (VLA):**
  - Operational since March 2023, recording all 27 antennas commensally during VLASS. Coherent beams cover about 2,000 sources per hour. Search uses seticore at SNR > 10 (8 in the post-processing paper) and ±50 Hz/s. Hits are logged to an SQL database with "postage stamp" voltage cutouts.
  - By 1 Sep 2024: more than 950,000 unique pointings, and 359,353,639 hits in the first six months. EIRP limits are 10^11–10^16 W.
  - Data availability: "can be accessed by contacting the lead author." Software is on GitHub (seticore, BLADE, CRICKETS, targets-minimal).

  [P] — [Tremblay et al. 2025](https://arxiv.org/abs/2501.17997); [Tremblay et al. 2023 (system paper)](https://arxiv.org/abs/2310.09414)
  - Inference: the COSMIC hit/candidate database is **not public**.
- COSMIC high-frequency SETI extension paper (Aug 2026) [P, abstract only] — [Myburgh et al. 2026](https://arxiv.org/abs/2608.18275)
- **ATA (SETI Institute):**
  - ATA-ASP anti-solar-point survey (since 2023): 28 dishes; 2,356–3,700 MHz (two 672 MHz tunings); 1.9 Hz and 16.8 s resolution; 5-min scans of 125 GB each.
  - Planned 87,660 pointings (≥ 11.25 PB). Search uses bliss (±4 Hz/s, SNR 15), then NBeamAnalysis; about 27 million hits were processed, about 11,000 plotted.
  - Software ATSAT is MIT-licensed (GitHub SETIatHCRO/ATSAT, Zenodo DOI 10.5281/zenodo.22751140). No public data release is stated.

  [P] — [Gerrard et al. 2026](https://arxiv.org/abs/2610.07065)
- ATA 3I/ATLAS search: 1–9 GHz, about 74 million hits in 7.25 h, AJ 2026 [P abstract via search] — [arXiv 2512.18142](https://arxiv.org/abs/2512.18142). Data are listed under `https://bldata.berkeley.edu/ATLAS/ATA_ATLAS` [P] — [seti.berkeley.edu/atlas](https://seti.berkeley.edu/atlas)
- The SETI Institute earlier released some ATA data in SigMF format; no date was found [S] — [SETI Institute news](https://seti.org/breakthrough-listen-and-uc-berkeley-seti-research-center-announce-open-format-data-release)
- The older ATA/IBM "setiQuest" data (ibmseti package, v2.0.0, PyPI last upload 2018-03-24) is legacy [P] — [PyPI ibmseti](https://pypi.org/project/ibmseti/)
- **FAST:**
  - FAST data policy: free-application projects have a 12-month protection period; priority and key projects are released in phases with up to 12 months proprietary each.
  - Quarterly public releases continue: DR23 (1 Nov 2025), DR24 (1 Feb 2026), DR25 (1 May 2026, covering Feb–Apr 2025 data). Use is restricted to scientific research, with no commercial use. Whether SETI data are included is unconfirmed.

  [S] — [FAST Data Release 25](https://lssf.cas.cn/en/facilities/sa/fast/notice/202605/t20260502_5108627.html); [FAST Data Policy](https://fast.bao.ac.cn/static/uploadfiles/数据中心相关文件/FAST-Data -Policy.pdf)
- FAST SETI papers in 2025–2026 (all [P, abstract]): no candidate table was confirmed as public.
  - M dwarfs multibeam — [Luan et al. 2025](https://arxiv.org/abs/2502.20419)
  - TRAPPIST-1 deep search — [Song et al. 2025](https://arxiv.org/abs/2509.06310)
  - Dual-backend multibeam position switching — [Li et al. 2025](https://arxiv.org/abs/2509.09654)
  - 3I/ATLAS — [Li et al. 2026](https://arxiv.org/abs/2603.19023)
  - Wavelet pipeline on 33 exoplanet systems — [Li et al. 2026](https://arxiv.org/abs/2605.23739)
- **LOFAR:** simultaneous dual-site (I-LOFAR Birr plus LOFAR-SE Onsala) search at 109.96–190.04 MHz; 2.98 Hz and 0.671 s resolution; 27,000,832 channels; 15-min observations; turboSETI ±4 Hz/s, S/N 10. Targets are TESS and Gaia stars over 232 deg² in 11 h; zero mutual hits. Full hit tables are "available with paper's supplementary material" [P] — [Johnson et al. 2023](https://arxiv.org/abs/2310.15704)
- LoFreq-SETI pipeline for the LWA Swarm (Oct 2026) [P, abstract] — [Taylor et al. 2026](https://arxiv.org/abs/2610.04126)
- **MWA:** opportunistic image-cube SETI searches (Galactic Centre, 400 deg², 115 min; Galactic anticentre; 'Oumuamua at 72–102 MHz; Vela region). These use MWA archive data; MWA data are public through the ASVO after a proprietary period [P for papers, abstracts] — [Tingay et al. 2016](https://arxiv.org/abs/1607.08818); [Tremblay & Tingay 2020](https://arxiv.org/abs/2009.03267). ASVO access terms [R].
- **UCLA SETI:**
  - 2020–2023 GBT L-band (1.15–1.73 GHz) search of 11,680 stars: ±8.86 Hz/s over 1,023 drift trials, S/N ≥ 10, 41.2 million detections. Limits: fewer than 6.6% of stars within 100 pc host a transmitter with EIRP > 10^13 W.
  - UCLA data feed the Zooniverse project "Are we alone in the universe?" (http://arewealone.earth).

  [P] — [Margot et al. 2023](https://arxiv.org/abs/2308.02712)
- UCLA candidate-detection tables on Zenodo (CC0) [P] (Zenodo API):
  - TRAPPIST-1/LHS 1140/Kepler field 2017 data: `GBT20170504.csv`, 1.07 GB — [Zenodo 3969036](https://zenodo.org/records/3969036)
  - 14 Kepler-field systems: 0.64 GB — [Zenodo 3897303](https://zenodo.org/records/3897303)
  - 31 Sun-like stars: 2.92 GB — [Zenodo 4023234](https://zenodo.org/records/4023234)
- **Zooniverse "Are we alone in the universe?":**
  - Launched 14 Feb 2023 with Planetary Society and NASA citizen-science funding. More than 20,000 volunteers and 700,000 classifications (by Aug 2024).
  - The project page reportedly says "Finished! … out of data" and anticipates new uploads in early 2026. Not confirmed (zooniverse.org blocked).

  [S] — [Zooniverse project](https://www.zooniverse.org/projects/ucla-seti-group/are-we-alone-in-the-universe)
- **SETI@home (final analysis, 2025):**
  - Data are Arecibo recordings over 14 years (the project ran 1998/1999–2020).
  - About 12 billion detections (spikes 5.03 billion, Gaussians 0.318 billion, pulses 3.08 billion, triplets 2.56 billion, autocorrelations 1.12 billion) were reduced by RFI removal to ranked "multiplets". A few hundred were selected; 80 were re-observed with FAST from 24 Sep 2022.
  - Source code is at https://sourceforge.net/projects/setiscience/ and the web interface ("Nebula") at https://setiathome.berkeley.edu/nebula/. No bulk downloadable candidate database is described.

  [P] — [Anderson et al. 2025](https://arxiv.org/abs/2506.14737); [Korpela et al. 2025](https://arxiv.org/abs/2506.14718)
- **Big Ear / Wow! signal:**
  - The Arecibo Wow! project retrieved "decades of previously unpublished Ohio SETI observations". It revised the Wow! signal to 1420.726 ± 0.005 MHz, peak flux above 250 Jy, at RA 19h25m02s ± 3s or 19h27m55s ± 3s, Dec −26°57′ ± 20′ (J2000).

  [P] — [Méndez et al. 2025, Arecibo Wow! II](https://arxiv.org/abs/2508.10657); [Méndez et al. 2024, Arecibo Wow! I](https://arxiv.org/abs/2408.08513)
  - The project states a goal to "digitize, catalog, and safeguard" Big Ear records. A data page exists (https://phl.upr.edu/wow/data, "last update September 2026"). The original printout is held by the Ohio History Connection (MSS 1151). Downloadability was not verified [S] — [Arecibo Wow! data](https://phl.upr.edu/wow/data)
  - The BL modern reobservation of the Wow! field is Perez et al. 2022 [P, as cited] — [MacMahon 2026](https://arxiv.org/abs/2610.07482)
- **Project Phoenix (SETI Institute, 1995–2004):** I found no public raw-data archive [R]. This is a gap.

### Inferences
- For a Stargazer prototype, the downloadable non-BL radio-SETI assets are hit/candidate tables (UCLA Zenodo CSVs; LOFAR supplementary tables) and a few per-campaign directories (3I/ATLAS). Raw voltages from COSMIC, BLUSE and ATA require contacting the teams.
- MWA (ASVO) and FAST (after the proprietary period) are general-purpose archives where SETI-useful spectral data could in principle be found. They need account registration, and FAST releases are for scientific use only.

### Gaps
- Kaggle, Zooniverse, FAST, MWA ASVO and Arecibo Wow! pages could not be fetched (blocked). Their terms are search-summary or recalled only.
- No public Project Phoenix or Big Ear raw-data archive was confirmed.
- SETI@home's per-candidate list (the few hundred selected multiplets) was not found as a downloadable table.

---

## 3. Peer-reviewed searches with machine-readable hit/candidate tables: parameters and where the tables live

### Takeaway
Most BL narrowband searches use turboSETI with about 3 Hz resolution, ±2 to ±4 Hz/s drift and S/N 10–25 on 5-min ABACAD cadences. Hit tables are public for Enriquez 2017 (CSV, CC BY), Price 2020 (event CSVs plus Zenodo `.dat` files) and Choza 2024 (AJ machine-readable table). Commensal-era searches (COSMIC, BLUSE, ATA-ASP) widen drift ranges to ±10–50 Hz/s but release data only on request.

### Cited findings (search parameters → table location)

| Search | Telescope / band | Targets | Drift range / S/N | Key results / EIRP | Machine-readable data |
|---|---|---|---|---|---|
| Enriquez et al. 2017 [P] — [arXiv 1709.03491](https://arxiv.org/abs/1709.03491) | GBT 1.1–1.9 GHz, about 3 Hz | 692 stars, about 4,800 5-min obs | ±2 Hz/s (400 steps, 0.01 Hz/s), S/N 25 for analysis (20 in the abstract) | about 29M hits; 11 events, all RFI; no transmitter > 10^13 W EIRP; Smin 17 Jy | Hit CSVs (CC BY 4.0) at [lband2017 downloads](https://seti.berkeley.edu/lband2017/downloads.html) |
| Price et al. 2020 [P] — [arXiv 1906.07750](https://arxiv.org/abs/1906.07750) | GBT L (1.10–1.90), S (1.80–2.70), Parkes 10-cm (2.60–3.45 GHz) | 1,327 stars | ±4 Hz/s (about 800 trials), S/N 10 | GBT L 37.14M hits, 15,998 events, 4,522 groups; S 10.12M hits; Parkes 4.45M hits; 0 final | [event CSVs](https://seti.berkeley.edu/listen2019/opendata.html); [Zenodo 7600403](https://zenodo.org/records/7600403) |
| Smith et al. 2021 [P] — [arXiv 2111.08007](https://arxiv.org/abs/2111.08007) | Parkes UWL 0.704–4.032 GHz | Proxima Cen | ±3 Hz/s (±4–5 Hz/s in sub-bands), S/N 10 | 4,172,702 hits, 5,160 events, 1 signal of interest (BLC1); EIRPmin 1.9 GW | "available for download via https://seti.berkeley.edu/blc1" |
| Sheikh et al. 2021 (BLC1) [P] — [arXiv 2111.06350](https://arxiv.org/abs/2111.06350) | Parkes | Proxima (982.002 MHz signal of interest) | — | BLC1 attributed to RFI (intermodulation) | data via the BLC1 page |
| Sheikh et al. 2019/2020 (drift-rate choice) [P] — [arXiv 1910.01148](https://arxiv.org/abs/1910.01148) | theory | — | recommends up to ±200 nHz | — | — |
| Gajjar et al. 2021 [P] — [arXiv 2104.14148](https://arxiv.org/abs/2104.14148) | GBT C (3.9–8 GHz) plus Parkes UWL (1–4 GHz), Galactic Centre | GC fields | ±4 Hz/s (860 / 1,720 steps), S/N 20 | EIRP ≥ 4×10^18 W (60M stars) and ≥ 5×10^17 W (0.5M stars); 304,658 non-zero-drift hits; SPANDAK dispersed-pulse search | "will be publicly available"; via the archive |
| Ma et al. 2023 (Nature Astronomy) [P] — [arXiv 2301.12670](https://arxiv.org/abs/2301.12670) | GBT 1.023–1.926 GHz | 1,004 cadences, 820 stars | β-VAE + random forest; snippets 4096 × 16, ±10 Hz/s; trained with setigen injections at S/N 10–50 | 20,515 signals of interest; 8 top candidates (MLc1–8) in Table 1 | input data from the open data archive; no full candidate table confirmed |
| Margot et al. 2023 (UCLA) [P] — [arXiv 2308.02712](https://arxiv.org/abs/2308.02712) | GBT 1.15–1.73 GHz | 11,680 stars | ±8.86 Hz/s, S/N 10 | 41.2M detections; < 6.6% of stars within 100 pc with EIRP > 10^13 W | Zooniverse; earlier UCLA tables on Zenodo (CC0) |
| Johnson et al. 2023 (LOFAR) [P] — [arXiv 2310.15704](https://arxiv.org/abs/2310.15704) | I-LOFAR plus LOFAR-SE, 110–190 MHz | TESS and Gaia stars in 232 deg² | ±4 Hz/s, S/N 10 | 0 mutual hits | supplementary tables (journal) |
| Choza et al. 2024 [P] — [arXiv 2312.03943](https://arxiv.org/abs/2312.03943) | GBT L, S, C, X (1.1–11.2 GHz) | 97 galaxy centres, 459 cadences | turboSETI ±4 Hz/s, input S/N 10 (effective threshold 33) | 6,034,313 hits, 1,519 events; EIRPmin about 2.3–3.1×10^24 W; documents turboSETI v2.3.2 S/N bias and drift-rate dechirping loss | "table is available in machine-readable form" (AJ) |
| Tremblay et al. 2025 (COSMIC, VLASS) [P] — [arXiv 2501.17997](https://arxiv.org/abs/2501.17997) | VLA 2–4 GHz | about 2,000 sources/h; 511 sources analysed in detail | ±50 Hz/s, SNR 8–10 | 359M hits; EIRPmin 2.3×10^11 – 2.1×10^16 W | on request |
| Painter et al. 2025 [P] — [arXiv 2412.05786](https://arxiv.org/abs/2412.05786) | GBT L, S, C, X archive | BL archive cadences | kurtosis-based anomaly search, ±4 Hz/s | EIRPmin 2–3×10^12 W (≤ 50 pc) | archive |
| Barrett et al. 2025 [P] — [arXiv 2506.13459](https://arxiv.org/abs/2506.13459) | Parkes UWL | 27 TESS eclipsing TOIs | ±0.1 to ±4 Hz/s, S/N 10 | 1,954,880 hits, 14,639 events; EIRPmin 9.5×10^11 – 1.3×10^15 W | via seti.berkeley.edu/opendata |
| Jacobson-Bell et al. 2025 [P] — [arXiv 2512.19763](https://arxiv.org/abs/2512.19763) | GBT 1–12 GHz | 3I/ATLAS (18 Dec 2025) | ±4 Hz/s, about 16σ effective | nondetection to about 100 mW–1 W EIRP | bldata.berkeley.edu/ATLAS/ |
| Czech et al. 2026 (BLUSE) [P] — [arXiv 2607.23651](https://arxiv.org/abs/2607.23651) | MeerKAT UHF/L/S (544 MHz–3.5 GHz) | about 1.2M pointings | ±10 Hz/s, SNR 6 | system paper; results forthcoming | on request |
| Gerrard et al. 2026 (ATA-ASP) [P] — [arXiv 2610.07065](https://arxiv.org/abs/2610.07065) | ATA 2.36–3.70 GHz | anti-solar point | ±4 Hz/s, SNR 15 | about 27M hits, about 11k plotted | software only |
| Anderson et al. 2025 / Korpela et al. 2025 (SETI@home) [P] — [2506.14737](https://arxiv.org/abs/2506.14737), [2506.14718](https://arxiv.org/abs/2506.14718) | Arecibo, about 1420 MHz ± 1.25 MHz [R for band] | nearly the whole Arecibo sky | barycentric and non-barycentric multiplets | about 12 billion detections → a few hundred candidates → 80 re-observed with FAST | code (SourceForge) and Nebula web pages only |

- Other recent BL-related analyses (abstract-level [P]):
  - Uno et al. 2023, extragalactic transmitter-rate limits — [arXiv 2304.02756](https://arxiv.org/abs/2304.02756)
  - Mason et al. 2025, "stellar bycatch" — [arXiv 2511.20231](https://arxiv.org/abs/2511.20231)
  - Suresh et al. 2023, periodic technosignatures in the Galactic Centre at 4–8 GHz — [arXiv 2305.18527](https://arxiv.org/abs/2305.18527)
- Choza et al. found that turboSETI v2.3.2 does not search exactly zero drift (zero-drift signals are reported at the smallest searched rate). Its S/N calculation (time-sum, discard top and bottom 5%, median/stdev) inflates S/N for setigen injections. Dechirping efficiency falls outside ±0.16 Hz/s [P] — [Choza et al. 2024](https://arxiv.org/abs/2312.03943)
- With about 3 Hz and 17 s data, a "BL-like" (turboSETI) process recovered only 570 of 10,000 injected signals over ±8.88 Hz/s, and 99.1% of recoveries were within ±1 Hz/s. The UCLA pipeline recovered 98.7% [P] — [Margot et al. 2023](https://arxiv.org/abs/2308.02712)

### Inferences
- Stargazer can use published hit tables (Enriquez CSV, Price CSV/Zenodo, Choza MRT, UCLA Zenodo) as ready-made "hit catalogues" for meta-analysis (cross-survey frequency/drift clustering, recurrence at the same sky position) without downloading raw spectra.
- Re-searching archive data at high drift rates is a known gap: turboSETI-style searches lose sensitivity above about 1 Hz/s on 18 s integrations.

### Gaps
- I could not open journal supplementary files (AJ/Nature/PASA sites blocked) to confirm column formats for Choza 2024, Johnson 2023 and Ma 2023.
- No 2025–2026 FAST, COSMIC or BLUSE hit catalogue is public as far as I found.

---

## 4. Labelled / benchmark datasets for training and testing detectors

### Takeaway
The main labelled set is the 2021 Kaggle "SETI Breakthrough Listen – E.T. Signal Search" competition data: GBT cadence snippets of 6 × 273 × 256 with injected signals, relaunched after a test-set leak. For controlled benchmarks, setigen (MIT) is the standard injection tool, and turboSETI recovers its injections on a laptop-class CPU (my test).

### Cited findings
- Kaggle 2021 competition (https://www.kaggle.com/competitions/seti-breakthrough-listen) [P] — [arXiv 2610.04449 review](https://arxiv.org/abs/2610.04449); [arXiv 2610.22050](https://arxiv.org/abs/2610.22050):
  - Data are GBT "cadence snippets (small sections of frequency across all six scans in a cadence)" with simulated signals injected into a fraction.
  - More than 750 teams entered; winners used ensembles of image CNNs. "A leak in the original test set forced a reset partway."
  - Each sample has "6 panels … (273 × 256 each)"; panels 0, 2, 4 are ON and 1, 3, 5 are OFF; labels mark injected-signal presence. Use the relaunched release.
- Further Kaggle details from search summary only [S] — [HeroX listing](https://www.herox.com/SETIBreakthroughListen):
  - $15,000 prize pool, 979 competitors in 768 teams.
  - Snippets are normalised after injection.
  - File count, size (I recall about 50k train samples as `.npy`, around 150 GB) and licence (Kaggle competition rules, likely non-commercial or competition-use only) are [R] and unconfirmed. kaggle.com is blocked here.
- Ma et al. 2023 built training data from 14,711 real GBT background snippets (3 cadences) with setigen injections (S/N 10–50, ±8 Hz/s) [P] — [Ma et al. 2023](https://arxiv.org/abs/2301.12670)
- setigen [P] — [PyPI setigen](https://pypi.org/project/setigen/); repository files `LICENSE`, `setup.py`; [Brzycki et al. 2022, AJ 163, 222](https://doi.org/10.3847/1538-3881/ac5e3d):
  - Latest PyPI release 2.7.0 (2024-08-02), MIT licence ("Copyright (c) 2018-2024 Bryan Brzycki").
  - Depends on numpy, scipy, astropy, blimpy ≥ 2.1.4 and matplotlib.
  - The GitHub repository is active: last commit 2026-04-26 (voltage-level polyphase filterbank/channelisation work).
  - It synthesises spectrogram frames and voltages and writes .fil/.h5.
- BL GBT README suggests Voyager observations and "birdie" injections as training material. Archive targets include `BIRDIE_001/002` and `BIRDIE_OFF_001/002` (injected test signals) and VOYAGER1 [P] — [GBT README](https://raw.githubusercontent.com/UCBerkeleySETI/breakthrough/master/GBT/README.md); [BI archive](https://breakthroughinitiatives.org/opendatasearch)
- My benchmark of setigen to turboSETI recovery (synthetic chi² noise; BL HSR geometry 2.79 Hz × 18.25 s × 16 integrations; 3 signals per coarse channel at S/N 15/30/100, drift ±3 Hz/s) [P] (my runs):
  - In the 8-coarse-channel run, turboSETI reported the injected signals with S/N and drift close to injected values (for example S/N about 117–120 for S/N 100, and about 16–19 for S/N 15).
  - In a synthetic ABACAD cadence (ET signal at +0.5 Hz/s in ON scans only, plus two RFI tones in all six), `find_event_pipeline(filter_threshold=3)` returned exactly one event: the injected ON-only signal at drift 0.500009 Hz/s, S/N 61.5. The zero-drift tone (S/N 196) and the −0.8 Hz/s RFI were rejected.
  - The event's `in_n_ons` field was 2, not the expected 3. Possible causes: drift or frequency window matching across scans, or a 0-based count. Worth checking before relying on the field.

### Inferences
- For Stargazer, setigen-injected frames into real archive backgrounds (as in Ma 2023 and Choza 2024) are a reproducible, licence-clean benchmark. The Kaggle set is useful but its licence terms must be checked before redistribution.

### Gaps
- Kaggle dataset size, file layout beyond the 6 × 273 × 256 shape, and licence text were not verified.
- No other public labelled radio-SETI benchmark with RFI-class labels was found. The ATA-ASP authors say about 75% of about 11,000 spectrograms are labelled, but the labels were not released.

---

## 5. Open-source software that reads and searches radio SETI data

### Takeaway
The Python stack (blimpy + turboSETI + setigen) still works in 2026 on CPU. turboSETI has had no release since May 2022. blimpy's repository got a 2.1.5 fix (Oct 2026) that is not yet on PyPI. GPU successors are hyperseti (CuPy), seticore (C++/CUDA, used at MeerKAT, VLA and Berkeley) and bliss (C++/CUDA with Python bindings, used by the SETI Institute's ATSAT). None of the three is on PyPI under its own name, so a tool would build them from source or need a GPU.

### Cited findings

| Package | Repository | Licence | Status (2026) | Requirements | What it does / limitations |
|---|---|---|---|---|---|
| blimpy | github.com/UCBerkeleySETI/blimpy | BSD (PyPI) | PyPI latest 2.1.4 (2022-11-16), 47 releases. Repo VERSION-HISTORY has 2.1.5 dated 2026-10-27 ("Replace deprecated pkg_resources with importlib.metadata"); last commit 2026-10-30; **2.1.5 not on PyPI** [P] | astropy, numpy ≥ 1.19, h5py, hdf5plugin, scipy, pandas, matplotlib; **pins pyparsing==2.4.7** [P] — [PyPI](https://pypi.org/project/blimpy/) | Reads and writes .fil, .h5 and GUPPI raw; `Waterfall` class, plotting, `fil2h5`/`h5fil`, `dsamp`, `bl_scrunch`. 2.1.4 imports `pkg_resources` (warns on setuptools ≥ 81, slated for removal). It rejects HDF5 with fewer than 3 time integrations (see §7) |
| turboSETI (`turbo-seti`) | github.com/UCBerkeleySETI/turbo_seti | MIT | PyPI 2.3.2 (2022-05-25), 43 releases; last commit on default branch 2022-05-25 [P] | astropy, numpy, blimpy ≥ 2.0.34, pandas, numba, dask; optional GPU via CuPy (`requirements_gpu.txt`) [P] — [PyPI](https://pypi.org/project/turbo-seti/) | Taylor-tree dedoppler per coarse channel → `.dat` hit tables; `find_event_pipeline` for ON/OFF cadence filtering; `plot_event_pipeline`. It still installs and runs with numpy 2.5 and Python 3.13 (my venv). Limits: no zero-drift search, S/N bias, dechirping loss above about 0.16 Hz/s ([Choza 2024](https://arxiv.org/abs/2312.03943)) |
| setigen | github.com/bbrzycki/setigen | MIT | PyPI 2.7.0 (2024-08-02); repo commits to 2026-04-26 [P] | numpy, scipy, astropy, blimpy, matplotlib | Signal synthesis and injection into spectrograms and raw voltages |
| hyperseti | github.com/UCBerkeleySETI/hyperseti | MIT (pyproject) | **Not on PyPI**; version 1.2.0 in pyproject; last commit 2025-12-15 ("Improved speed of _group_hits") [P] | Python ≥ 3.8, **CuPy ≥ 10 (CUDA GPU required)**, pandas ≥ 2, h5py, hdf5plugin, setigen | Brute-force GPU dedoppler with a "ddsk" kernel (dedoppler plus spectral kurtosis), boxcar trials, hit browser; CLI `findET`, `tsdat`, `hyperfrost`. Described as "spiritual successor to turboSETI" |
| seticore | github.com/lacker/seticore (also UCBerkeleySETI/seticore) | MIT ("Copyright (c) 2022 Kevin Lacker") | Not on PyPI; lacker clone last commit 2023-08-03. The UCBerkeleySETI copy could not be fetched (newer state unknown) [P] | C++ (meson/ninja, Boost, HDF5 with the bitshuffle plugin), CUDA GPU | High-performance turboSETI-compatible dedoppler (.h5 → .dat) plus beamforming (.raw + .bfr5 → filterbank). Used at MeerKAT (BLUSE), VLA (COSMIC) and Berkeley. Writes hits and "stamps" |
| bliss (`blissdedrift`) | github.com/UCBerkeleySETI/bliss (mirror n-west/bliss) | BSD-3-Clause (LICENSE "Copyright (c) 2026, Nathan West, Berkeley SETI Research Center") | Last commit 2026-03-20 (CITATION.cff added); Zenodo release "2025-0522"; pyproject name `blissdedrift`, "Development Status :: 3 - Alpha". **Not on PyPI** (the PyPI package `bliss` is an unrelated ESRF project, LGPL) [P] | C++17, CMake, libhdf5, FFTW, bitshuffle; optional CUDA; Python ≥ 3.8 bindings; Linux | Cadence-aware flagging (filter rolloff, spectral kurtosis), noise estimation, drift integration (example ±500 Hz/s with desmear), hit search. Used by SETI Institute ATSAT ([Gerrard 2026](https://arxiv.org/abs/2610.07065)) |
| rawspec | github.com/UCBerkeleySETI/rawspec | — (not checked) | Last commit 2024-05-21 [P] | CUDA/cuFFT GPU | GUPPI RAW → SIGPROC .fil or FBH5 .h5 (`-j`), multiple products at once (default FFT lengths 1048576, 8, 1024 and integrations 51, 128, 3072, i.e. the standard HSR/HTR/MR products); UDP output |
| BLADE | github.com/luigifcruz/blade | — | last commit 2026-09-13; Zenodo DOI 10.5281/zenodo.21454459 cited [P] | CUDA | BL Accelerated DSP Engine (beamforming and upchannelisation at ATA and COSMIC) |
| ATSAT, NBeamAnalysis | github.com/SETIatHCRO/ATSAT; isabelgerrard/NBeamAnalysis | MIT (ATSAT) | ATSAT 2.0.0 (2026, Zenodo 10.5281/zenodo.22751140) [P] | bliss | Multibeam sky-localisation filtering of hits |
| Blio.jl | github.com/david-macmahon/Blio.jl | — | last commit 2026-07-26 [P] | Julia | Julia readers for BL formats |
| sigpyproc (sigpyproc3) | github.com/FRBs/sigpyproc3 | MIT | PyPI 2.0.0 (2026-10-17), Python ≥ 3.12 [P] — [PyPI](https://pypi.org/project/sigpyproc/) | numpy ≥ 2, numba | General SIGPROC filterbank toolbox (pulsar/FRB-oriented; useful generic .fil reader) |
| your | github.com/thepetabyteproject/your | — | PyPI 0.6.7 (2022-01-22) [P] | — | Unified filterbank/PSRFITS reader |
| baseband | github.com/mhvk/baseband | GPL v3+ | PyPI 4.3.0 (2025-03-10) [P] | astropy | Baseband readers (VDIF, GUPPI, …) |
| blscint | github.com/bbrzycki/blscint | — | PyPI 0.0.1 (2022) [P] | setigen | Scintillation utilities for technosignature filtering |
| hdf5plugin / h5py | — | — | hdf5plugin 7.1.0 (2026-09-01); h5py 3.16.0 (2026-03-06) [P] | — | Required for bitshuffle-compressed BL .h5. blimpy 2.1.4 triggers hdf5plugin deprecation warnings ("Bitshuffle's lz4 argument is deprecated") [P] (my run) |
| SETI@home science code | sourceforge.net/projects/setiscience/ | — | [P] cited in [Anderson 2025](https://arxiv.org/abs/2506.14737) | — | RFI removal and multiplet finding |

- seticore implements "turboSETI's dedoppler algorithm" and "should already be installed on the BL machines on MeerKAT, the VLA, and Berkeley". HDF5 plugin set-up requires copying hdf5plugin `.so` files into the HDF5 plugin directory [P] — [seticore README](https://raw.githubusercontent.com/lacker/seticore/master/README.md)
- hyperseti "uses the GPU heavily, so a working CUDA environment is needed" [P] — [hyperseti README](https://raw.githubusercontent.com/UCBerkeleySETI/hyperseti/master/README.md)
- Not found on PyPI under these names: hyperseti, seticore, rawspec, bliss-seti and blio [P] (PyPI JSON API 404).

### Inferences
- For a CPU-only open-source tool: blimpy (reader) plus turboSETI (baseline search) plus setigen (tests) is the only pip-installable chain. Pin numpy and setuptools carefully, since blimpy pins pyparsing 2.4.7 and 2.1.4 uses `pkg_resources`. Consider writing a thin FBH5/SIGPROC reader directly on h5py + hdf5plugin to avoid blimpy's quirks.
- For scale or high-drift searches, bliss (BSD, CPU or CUDA) is the most actively maintained, permissively licensed engine. It needs a C++ build.

### Gaps
- The current UCBerkeleySETI/seticore and UCBerkeleySETI/turboSETI default-branch state could not be fetched (the earlier fetch failed). The rawspec licence was not checked.
- The GitHub web/API was not reachable, so commit dates come from 14 Nov 2026 clones.

---

## 6. General radio sky surveys with public multi-epoch data that could be mined for anomalous or variable sources

### Takeaway
The main public multi-epoch radio image/catalogue surveys are VLASS (2–4 GHz, three or more epochs, NRAO/CIRADA), ASKAP RACS and VAST (CASDA), LoTSS DR2 (LOFAR 120–168 MHz, single-epoch but large) and MeerKAT archive/survey products (SARAO). These are continuum imaging surveys at arcsecond resolution, not Hz-resolution spectra. They suit searches for variable or anomalous compact sources, not narrowband carriers. Access needs accounts or TAP/VO services I could not test from this sandbox.

### Cited findings
- VLASS: 2–4 GHz, about 2.5″ resolution, 70 µJy/beam coadded goal, all sky above Dec −40° (33,885 deg²) in three epochs, about 5,500 h; started Sept 2017 [P] — [Lacy et al. 2020](https://arxiv.org/abs/1907.01981)
- Epoch 3 Quick Look interim component lists were released before Epoch 4.1 (NRAO memo, June 2025). CIRADA catalogues "transitioned to NRAO". VLASS 2.1/3.1 Quick Look images are preferred over 1.1 because of pointing-calibration problems [S] — [VLASS memo 022](https://library.nrao.edu/public/memos/vla/vlass/VLASS_022.pdf); [NRAO VLASS basic data products](https://science.nrao.edu/vlass/vlass-data/basic-data-products)
- COSMIC runs commensally during VLASS, so VLASS pointings are also the COSMIC SETI footprint [P] — [Tremblay 2025](https://arxiv.org/abs/2501.17997)
- RACS: full 36-antenna ASKAP, whole sky south of about Dec +41°, 700–1800 MHz; "All RACS survey products will be public" (about 15″ resolution) [P] — [McConnell et al. 2020](https://arxiv.org/abs/2012.00747)
- VAST: ASKAP variables and slow transients, 5 s to about 5 yr timescales. The Pilot covered about 162 h at 888 MHz, rms 0.24 mJy/beam, 113 fields, 12-min epochs [P] — [Murphy et al. 2021](https://arxiv.org/abs/2108.06039)
- VAST Extragalactic Survey Data Release 1 came out in Feb 2026 [P, abstract] — [de Ruiter et al. 2026](https://arxiv.org/abs/2602.22739)
- vast-tools (Python, GitHub askap-vast/vast-tools, last commit 2026-07-20) queries VAST pilot data, light curves and forced photometry. "Some features require access to the VAST Pilot Survey data" [P] — vast-tools README (repository clone)
- ASKAP images and catalogues come from CASDA [S] — [arXiv 2503.09108](https://arxiv.org/pdf/2503.09108)
- LoTSS DR2: 120–168 MHz, 27% of the northern sky (4,178 + 1,457 deg²), from 3,451 h (7.6 PB) of HBA data; 4,396,228 sources [P] — [Shimwell et al. 2022](https://arxiv.org/abs/2202.11733)
- MeerKAT: BL and RATT are mining the MeerKAT archive for transients (TRON). katdal (BSD-3; PyPI 0.24, 2026-11-13) is the MeerKAT visibility data-access library [P] — [BI news 52](https://breakthroughinitiatives.org/news/52); [PyPI katdal](https://pypi.org/project/katdal/)
- MWA data via ASVO; giant-squid (Rust CLI for the MWA ASVO, last commit 2026-10-06) [P for the repo date] — MWATelescope/giant-squid clone. ASVO registration [R].

### Inferences
- For "anomalous source" mining, VLASS epochs 1–3+ (northern sky) and VAST/RACS (southern sky) give the widest public multi-epoch coverage. Cross-matching variable compact sources against nearby-star catalogues (Gaia) is feasible on catalogues alone (tens of GB), without images.

### Gaps
- VizieR (vizier.cds.unistra.fr) returned "Page unavailable" today, so I could not confirm catalogue IDs or row counts. CASDA, CIRADA, LOFAR Surveys and SARAO endpoints are blocked here. Access methods (TAP URLs, registration) are not verified.

---

## 7. Practical feasibility: smallest useful slice and turboSETI run time on modest hardware

### Takeaway
One coarse channel of BL HSR data (2.93 MHz × 16 integrations, about 112 MB as .h5) is searched by turboSETI in about 8 s on a 4-vCPU VM. That scales linearly to about 8 s per coarse channel, so a full ABACAD cadence of one compute node (6 × 64 coarse channels) takes about 50 min, and a full spliced L-band cadence about 7 h of CPU time. The smallest useful real-data slice is one ABACAD cadence of single-node HSR .h5 files (6 × about 1–2 GB). The public Voyager .fil sample is not usable with turboSETI 2.3.2 / blimpy 2.1.4 because it has only 2 integrations.

### Cited findings (my benchmarks on 4 vCPU, 15 GB RAM, no GPU; turboSETI 2.3.2, blimpy 2.1.4, h5py 3.16.0, hdf5plugin 7.1.0, Python 3.13, numpy 2.5.4) [P]
- 1 coarse channel (1,048,576 × 16, 2.79 Hz, 18.25 s; file 112 MB .h5), max_drift ±4 Hz/s, S/N 10: search time 0.13 min (about 8 s). Drift resolution 0.0102 Hz/s.
- 8 coarse channels (8,388,608 channels; 895 MB .h5): 1.09 min, i.e. about 8.2 s per coarse channel. 24 injected signals were reported with correct drift and S/N.
- Synthetic ABACAD cadence (6 × 112 MB): FindDoppler 0.15 min per scan, then `find_event_pipeline` in seconds. It isolated the ON-only injected signal.
- The real public sample `voyager_f1032192_t300_v2.fil` (528 MB, 66,060,288 channels, 2 integrations of 105.7 s) failed:
  - turboSETI first converts .fil to .h5 (370 MB written).
  - blimpy `examine_h5` then reads `h5["data"][2]`, which fails ("Index (2) out of range"), and the tool exits with code 86 ("HDF5 Time-dimension data corruption").
  - Inference: turboSETI/blimpy need at least 3 time integrations, and this sample is a header/plotting demo, not a drift-search test file.
- Scale reference:
  - A GBT 5-min HSR L-band spliced file in the archive is about 15.3 GB compressed .h5 (e.g. HIP56445 `.gpuspec.0000.h5`, 15,271,502,711 B); the matching MR file is about 243 MB [P] — [BI archive](https://breakthroughinitiatives.org/opendatasearch)
  - One node covers 187.5 MHz = 64 coarse channels [P] — [GBT README](https://raw.githubusercontent.com/UCBerkeleySETI/breakthrough/master/GBT/README.md)
  - Price 2020 used about 800 drift trials at ±4 Hz/s [P] — [Price 2020](https://arxiv.org/abs/1906.07750)
- The newer GBT rawspec HDF5 in the archive is split per compute node (e.g. `blc16_guppi_59902_..._NGC3226_0056.rawspec.0000.h5`, 3.2 GB), which gives a natural per-node download unit [P] — [BI archive](https://breakthroughinitiatives.org/opendatasearch)

### Inferences
- Extrapolated from about 8 s per coarse channel (not measured on real data):
  - Per-node file (64 coarse channels): about 9 min per scan, about 50 min per ABACAD cadence.
  - Full spliced L-band file (about 512 coarse channels over about 1.5 GHz recorded, about 800 MHz usable): about 70 min per scan, about 7 h per cadence.
  - RFI-rich real data produce more hits and slower event processing than synthetic noise.
- Download is the bottleneck: at the 27 MB/s I measured from GCS, a 15 GB file takes about 10 min, and a full L-band cadence (6 × 15 GB = 90 GB) about 1 h.
- Recommended minimal slice: pick one cadence via the archive GET search (`project=GBT&file_type=HDF5&target=<HIP>`), download only the `.0002` MR files (about 0.24 GB each) for energy detection, or per-node `rawspec.0000.h5` HSR files (about 1–3 GB) for drift search. Or use Price/Enriquez hit tables and skip spectra entirely.
- Memory: peak RSS was not captured for the 8-channel run. The 15 GB VM handled 895 MB files easily, and turboSETI processes one coarse channel at a time.

### Gaps
- No measurement on a real multi-integration BL .h5 cadence, because the Berkeley hosts holding them are blocked from this sandbox. The timings are synthetic-data timings.
- No GPU benchmarks (hyperseti, seticore, bliss) were possible here.
