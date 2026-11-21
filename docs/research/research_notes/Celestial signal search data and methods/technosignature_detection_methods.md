# Technosignature detection and encoded-signal testing: methods reference (state of the art to November 2026)

> **Verification status — read before using these notes.** Sections 1–3 were drafted on 14 November 2026 when almost nothing could be fetched; on 21 November 2026 they were re-checked against arXiv abstracts and, for the load-bearing numbers, the paper PDFs. Sections 4–10 were written on 21 November 2026. Tags:
> - **[P]** read on 21 November 2026 in the primary source: the arXiv abstract or full text, a Crossref/Semantic Scholar record, or the text of a protocol document.
> - **[P†]** read in a primary source (package source code or metadata downloaded from PyPI) in the earlier session on 14 November 2026 and not re-read today. I kept these because package code is a primary source, but versions may have moved.
> - **[S]** from a web-search result summary only (mainly 2025–2026 news and protocol updates). The underlying page was not read.
> - **[M]** recalled from background knowledge and **not re-checked**. Treat these as leads.
> - Where a numeric detail in the old notes turned out to be wrong, the correction is marked **(corrected 21 Nov)**.
> - Links of the form `arxiv.org/search/?query=...` are search links, given where the exact identifier was not confirmed. They are not verified identifiers.
> - Arithmetic derived here from cited values is placed under "Inferences" and marked "(computed)".
> - Hosts that could not be reached from the sandbox (for example ADS, journal websites, iaaspace.org, st-andrews.ac.uk) are recorded as gaps. Content from those hosts appears only as [S].

## 1. Frameworks: how the field classifies technosignatures and search strategies, and which classes can be searched in public archival data

### Takeaway
The field classifies searches by signature class (deliberate or leaked communication signals; energy-use or "Dysonian" artefacts such as waste heat and megastructures; atmospheric and surface signatures; Solar-System artefacts), by figures of merit (Sheikh's nine axes; the Wright–Kanodia–Lubar haystack fraction, about 6×10⁻¹⁸ for radio SETI to 2018) and by target breadth (the Exotica catalog). Public archives support these searches: narrowband, pulsed and periodic radio; optical and infrared time-domain anomalies; infrared excess; laser lines in archival spectra; vanishing sources; and synchronisation ("SETI ellipsoid") targeting. They do not support nanosecond optical-pulse searches, phase-level modulation analysis of most radio data, or atmospheric technosignatures at useful sensitivity.

### Cited Findings

#### Classification frameworks and figures of merit
- [P] The NASA Technosignatures Workshop (Lunar and Planetary Institute, Houston, September 2018) produced the report "NASA and the Search for Technosignatures", convened so that NASA could learn about the state of the art and what role it might play. [M] It reviews limits by class: radio and optical/near-infrared communication, megastructures and waste heat, atmospheric and surface signatures, and Solar-System artefacts. — [NASA Technosignatures Workshop Participants 2018, arXiv:1812.08681](https://arxiv.org/abs/1812.08681)
- [M] Tarter's review framed SETI as a search of a nine-dimensional "cosmic haystack" (three spatial dimensions, time, two polarisations, frequency, sensitivity and modulation) and is the standard reference for classical radio and optical strategies (targeted, sky survey and commensal). — [Tarter 2001, ARA&A 39, 511](https://doi.org/10.1146/annurev.astro.39.1.511)
- [P] Sheikh proposes a framework for comparing technosignature searches on nine "axes of merit". [M] The axes are observing capability, cost, ancillary benefits, detectability, duration, ambiguity, extrapolation, inevitability and information. — [Sheikh 2019/2020, IJA; arXiv:1908.02683](https://arxiv.org/abs/1908.02683)
- [P] Wright, Kanodia and Lubar turned the haystack metaphor into a quantitative eight-dimensional model, evaluated the integral analytically, and released a Python script for computing haystack volumes. The total model volume is 6.4×10¹¹⁶ m⁵ Hz² s W⁻¹. The 692-star survey of Enriquez et al. covered a fraction of 3.8×10⁻¹⁹. **(corrected 21 Nov)** The sum over the large radio programmes is **6.0×10⁻¹⁸** in the arXiv text (v3), about 8,000 litres of seawater out of Earth's oceans, "between the volumes of a large hot tub and a small swimming pool". The 5.8×10⁻¹⁸ figure in the old notes was not found in the arXiv version and may come from the journal version. — [Wright, Kanodia & Lubar 2018, AJ 156, 260; arXiv:1809.07252](https://arxiv.org/abs/1809.07252)
- [P] The Breakthrough Listen Exotica Catalog introduces "survey breadth". It is "a 963 entry collection of 816 distinct targets" in four samples: Prototype, Superlative, Anomaly and Control. It comes with a new classification system for anomalies. — [Lacki et al. 2021, ApJS 257, 42; arXiv:2006.11304](https://arxiv.org/abs/2006.11304)
- [P] Socas-Navarro et al. (TechnoClimes 2020) set out the community's concepts for future technosignature missions. [M] The paper introduces the "ichnoscale" (the size of a signature in units of present-day Earth's). — [Socas-Navarro et al. 2021, Acta Astronautica 182, 446; arXiv:2103.01536](https://arxiv.org/abs/2103.01536)
- [M] Wright et al. argue that technosignatures may be more abundant, longer-lived, more detectable and less ambiguous than biosignatures. — [Wright et al. 2022, ApJL 927, L30 (search link)](https://arxiv.org/search/?query=Wright+technosignatures+abundant+long-lived+detectable+unambiguous&searchtype=all)
- [P] The Ad Hoc Committee on SETI Nomenclature (convened in 2018 at Frank Drake's suggestion) recommended standardised definitions of terms used inconsistently in the literature. These are useful as labels in software. — [Wright et al. 2018, arXiv:1809.06857](https://arxiv.org/abs/1809.06857)
- [M] Standard radio figures of merit for null results are: EIRP_min = 4πd²·F_min, with F_min = (S/N)_min·SEFD·√(δν/(n_pol·τ_obs)) for a signal narrower than the channel; the Drake figure of merit; and the continuous-waveform transmitter figure of merit (CWTFM). — [Enriquez et al. 2017, arXiv:1709.03491](https://arxiv.org/abs/1709.03491)
- [P] Wright's 2021 guide to SETI strategy focuses "not on the variety of technosignatures… but on" the philosophies and strategies behind searches. — [Wright 2021, arXiv:2107.07283](https://arxiv.org/abs/2107.07283)
- [P] **2026 development:** Vidal, Fields, Sowinski, Bartlett et al. (posted 20 May 2026, v2) review proposed technosignatures, organised by spatial scale "from the Earth and the solar system to the edge of the Universe", deliberately not tied to specific observatories. Preprint; peer-review status not checked. — [Vidal et al. 2026, arXiv:2605.21093](https://arxiv.org/abs/2605.21093)
- [P] **2026 development:** Haqq-Misra (18 March 2026) defines a two-parameter matrix: domain accessibility (can the signature be reached within today's spatial domain?) and recognisability. The four cells are four search strategies. — [Haqq-Misra 2026, arXiv:2603.17741](https://arxiv.org/abs/2603.17741)

#### Technosignature classes that have been searched with archival or public data
- [P] **Narrowband drifting radio carriers.** Breakthrough Listen's first data release (about 219 TB, for 1327 stars) is public. — [Price et al. 2020, arXiv:1906.07750](https://arxiv.org/abs/1906.07750); [Lebofsky et al. 2019, arXiv:1906.07391](https://arxiv.org/abs/1906.07391)
- [P] **Broadband pulsed, artificially (negatively) dispersed beacons:** 1883 stars, 233 h, 4–8 GHz. — [Gajjar et al. 2022, arXiv:2205.02964](https://arxiv.org/abs/2205.02964)
- [P] **Periodic pulse trains:** a fast folding algorithm (FFA) search of the Galactic Centre. — [Suresh et al. 2023, arXiv:2305.18527](https://arxiv.org/abs/2305.18527)
- [P] **Anomalous optical light curves:** Kepler, TESS and ZTF. **(corrected 21 Nov)** The Giles & Walkowicz identifier in the old notes (1812.01150) points to an unrelated mathematics paper; the correct one is arXiv:1812.07156. — [Giles & Walkowicz 2018, arXiv:1812.07156](https://arxiv.org/abs/1812.07156); [Malanchev et al. 2021, arXiv:2012.01419](https://arxiv.org/abs/2012.01419); [Martínez-Galarza et al. 2021, arXiv:2009.06760](https://arxiv.org/abs/2009.06760)
- [P] **Infrared excess:** Gaia + 2MASS + WISE for partial Dyson spheres, and WISE for galaxy-scale waste heat. — [Suazo et al. 2024, arXiv:2405.02927](https://arxiv.org/abs/2405.02927); [Griffith et al. 2015, ApJS 217, 25; arXiv:1504.03418](https://arxiv.org/abs/1504.03418)
- [P] **Continuous laser lines in archival high-resolution spectra:** Keck HIRES (5600 stars) and APF. — [Tellis & Marcy 2017, arXiv:1704.02535](https://arxiv.org/abs/1704.02535); [Zuckerman et al. 2023, arXiv:2301.06971](https://arxiv.org/abs/2301.06971)
- [P] **Vanishing or appearing sources:** USNO-B1.0 against Pan-STARRS, and POSS-I against Pan-STARRS and Gaia. — [Villarroel et al. 2020, arXiv:1911.05068](https://arxiv.org/abs/1911.05068); [Solano et al. 2022, arXiv:2206.00907](https://arxiv.org/abs/2206.00907)
- [P] **SETI ellipsoid:** Gaia target selection, then TESS light curves. — [Davenport et al. 2022, arXiv:2206.04092](https://arxiv.org/abs/2206.04092); [Cabrales et al. 2024, arXiv:2402.11037](https://arxiv.org/abs/2402.11037)
- [P] **Spatio-temporal survey SETI** (Astro2020 white paper). — [Davenport 2019, arXiv:1907.04443](https://arxiv.org/abs/1907.04443)
- [M] **High-energy and gravitational-wave channels** have been examined only at proof-of-concept level. — [Hippke & Forgan, X-ray (search link)](https://arxiv.org/search/?query=Hippke+Forgan+X-ray+SETI&searchtype=all); [Sellers et al. 2022, gravitational waves (search link)](https://arxiv.org/search/?query=Sellers+gravitational+wave+technosignature&searchtype=all)
- [P] **Nanosecond optical pulses** need dedicated fast photon-counting hardware. VERITAS analysed 30 h of dedicated observations of 136 targets and 249 archival observations of 140 targets, with no detection. — [Acharyya et al. 2023, AJ 166, 84; arXiv:2306.17680](https://arxiv.org/abs/2306.17680)
- [M] **Atmospheric technosignatures** (NO₂, CFCs) are modelled targets for future spectroscopy, not archive searches. — [Kopparapu et al. 2021 (search link)](https://arxiv.org/search/?query=Kopparapu+nitrogen+dioxide+technosignature&searchtype=all)

### Inferences
- An archive-only tool works on five data shapes: (1) radio dynamic spectra; (2) evenly sampled space photometry; (3) irregular ground-based light curves; (4) catalogue spectral energy distributions and astrometry; (5) spectra and digitised images. Each needs its own detector and null model.
- Narrowband radio scores well on Sheikh's "ambiguity" axis, because its confounder is human, not astrophysical. Anomaly, infrared-excess and vanishing-source searches score poorly on ambiguity. Every published candidate in those classes has needed, and mostly received, a natural or instrumental explanation (Sections 4 and 7).
- Log the covered haystack (targets, frequency or time span, sensitivity, drift or period range, and pipeline efficiency) for every run. Null results can then be stated as upper limits in the same currency as Wright et al. 2018 and Margot et al. 2023.

### Gaps
- The NASA workshop report's category list and its archive-mining recommendations were not re-read in full; only the abstract was.
- The Vidal et al. 2026 review and the Haqq-Misra matrix were read only at abstract level.
- The origin of the 5.8×10⁻¹⁸ haystack figure (versus 6.0×10⁻¹⁸ in the arXiv text) is unresolved. The journal version could not be reached.

## 2. Radio narrowband searches: de-Doppler search, drift rates, thresholds, ON/OFF cadences, coincidence rejection, scintillation, machine learning, false-positive rates and BLC1

### Takeaway
The established method is an incoherent tree de-Doppler search of fine-resolution spectrograms (about 3 Hz and about 18 s), with S/N thresholds of 10–25 over drift rates of ±2 to ±4 Hz/s, followed by ON/OFF cadence filtering. Published funnels run from 10⁶–10⁸ hits to tens of events to zero confirmed candidates, so false positives are set by human interference, not by noise. BLC1 showed that a signal can pass cadence filtering for hours and still be local intermodulation interference. Injection tests show that the standard fine-resolution product recovers only about 6% of injected signals.

### Cited Findings

#### Why narrowband, and what counts as artificial
- [M] The narrowest known natural radio emitters are masers, with line widths of a few hundred Hz, so a signal of a few Hz is treated as artificial; the confounder is interference. — [Tarter 2001](https://doi.org/10.1146/annurev.astro.39.1.511); [Enriquez et al. 2017](https://arxiv.org/abs/1709.03491)
- [P] **New natural confounder (2024–2025).** Cold H I clouds produce narrowband (about 10 kHz) emission near 1420 MHz resembling the Wow! signal, and a maser-like flare of such a cloud is proposed as "a new source of false positives in technosignature searches". These are kHz-wide, not Hz-wide (Section 10). — [Méndez et al. 2024, arXiv:2408.08513](https://arxiv.org/abs/2408.08513)

#### The de-Doppler (tree) algorithm as implemented in turboSETI
- [P†] turboSETI searches for narrowband drifting signals in filterbank data, using a Taylor tree in Numba (with optional CuPy). It reads `.fil`/`.h5` and writes `.dat` hit tables. Latest PyPI release 2.3.2 (25 May 2022). — [turbo-seti on PyPI](https://pypi.org/project/turbo-seti/2.3.2/)
- [P†] Code defaults: `max_drift=10.0` Hz/s, `min_drift=0.00001`, `snr=25.0`, `blank_dc=True`. The README example uses `max_drift=4, snr=25`. — [turbo-seti 2.3.2 source](https://pypi.org/project/turbo-seti/2.3.2/)
- [P†] The tree sums along drift paths in log₂(N_t) stages; N_t must be a power of two; drift resolution is |δν|/((N_t−1)·δt). — [turbo-seti 2.3.2 source](https://pypi.org/project/turbo-seti/2.3.2/)
- [P†] The example Breakthrough Listen file has δt = 18.253611 s, δν = 2.793968 Hz and 16 integrations (about 292 s), with 1,048,576 fine channels per coarse channel. — [turbo-seti README](https://pypi.org/project/turbo-seti/)
- [P†] Robust noise statistics: the median and standard deviation are computed after clipping to the 5th–95th percentile. S/N = (sum − median)/σ. The maximum over drifts is kept per channel. A "top hit" must be the local maximum within ±2·ceil(max_drift/drift resolution) channels. — [turbo-seti 2.3.2 source](https://pypi.org/project/turbo-seti/2.3.2/)
- [M] The tree algorithm comes from pulsar dedispersion (Taylor 1974): O(N_t log N_t) per channel instead of O(N_t²). — [Taylor 1974, A&AS 15, 367 (search link)](https://arxiv.org/search/?query=Taylor+1974+dispersed+radio+emission&searchtype=all)

#### Thresholds, drift ranges and outcomes of the main surveys
- [P] Enriquez et al. 2017: 692 stars, Green Bank 1.1–1.9 GHz, three 5-minute ON scans interspersed with 5-minute OFF scans. Eleven events passed thresholds, all consistent with interference. No high-duty-cycle transmitters at about 10¹³ W EIRP; fewer than about 0.1% of systems within 50 pc host such transmitters. [M] S/N ≥ 25, drift ±2 Hz/s. — [arXiv:1709.03491](https://arxiv.org/abs/1709.03491)
- [P] Price et al. 2020: 1138 stars at Green Bank (1.10–1.90 and 1.80–2.80 GHz) and 189 at Parkes (2.60–3.45 GHz), drift −4 to +4 Hz/s, zero candidates. Limits are 2×10¹² W (GBT) and 9×10¹² W (Parkes); DR1 is about 219 TB. **(corrected 21 Nov)** Funnel from the full text: **51.71 million hits → 21,117 events detected only in ON observations → 6,154 event groups**. The old "about 6,150 → 20 events" was partly wrong: 6,154 is the number of clustered event groups, and the "20 events in all ON scans" figure was not found. [M] S/N ≥ 10. — [arXiv:1906.07750](https://arxiv.org/abs/1906.07750)
- [P] Smith et al. 2021 (Proxima Centauri, Parkes UWL): **4,172,702 hits** in total, of which **5,160** were present in multiple on-source pointings but not in off-source pointings. Only one event passed all filtering and visual inspection: BLC1 at **982.002571 MHz**, drift **0.038 Hz/s**, detected over about 2.5 h. — [Smith et al. 2021, Nature Astronomy; arXiv:2111.08007](https://arxiv.org/abs/2111.08007)
- [P] UCLA pipeline (Margot et al. 2021; 31 Sun-like stars; Green Bank L-band 1.15–1.73 GHz): detection by topographic prominence; direction-of-origin filters that uniquely link signals across scans. Injection recovery was 93% (98% outside dense interference); 99.73% of recovered injections were correctly classified as candidates; over 99.84% of about 26 million signals were classified as interference; 4,539 candidates outside known interference regions were inspected and all were anthropogenic. — [arXiv:2011.05265](https://arxiv.org/abs/2011.05265)
- [P] Margot et al. 2023 (62 TESS Objects of Interest, about 11,680 stars in the beam): the UCLA pipeline recovers 94.0% of injections (98.7% excluding dense interference), while a pipeline using incoherent sums of 51 consecutive spectra recovers about 15 times fewer, about 6%. They introduce an improved Drake figure of merit including pipeline efficiency and duty cycle. At 95% confidence, fewer than 6.6% of stars within 100 pc host a transmitter with EIRP > 10¹³ W, and at most 3×10⁻⁴ of stars within 20,000 ly host one above 5×10¹⁶ W. The pipeline natively detects the Ma et al. 2023 signals. — [arXiv:2308.02712](https://arxiv.org/abs/2308.02712)
- [P] **2026 development (imaging-domain search):** OVRO-LWA all-sky narrowband imaging at 50–86 MHz with about 10 Hz channels (more than 3×10⁶ images per 30 s epoch) used multi-kernel matched filtering, empirical noise standardisation **and false-discovery-rate control**. Three candidates above 10σ were all resolved as non-celestial. Sensitivity about 100 Jy per channel; EIRP limits about 10¹⁴ W at 10 pc. — [Kosogorov et al. 2026, arXiv:2606.04304](https://arxiv.org/abs/2606.04304)
- [P] **2025 development:** the SETI@home final analysis paper describes stage two: interference removal, then grouping detections spread over 14 years of Arecibo data into candidates. Sensitivity was calibrated with injected artificial "birdies" that model persistent ET signals over a range of power, bandwidth and planetary-motion parameters. A few hundred top candidates were selected for re-observation with FAST. — [Anderson et al. 2025, arXiv:2506.14737](https://arxiv.org/abs/2506.14737)
- [P] **2025 development:** a Breakthrough Listen observation of the interstellar object 3I/ATLAS on 18 December 2025 (GBT, 1–12 GHz) found no candidates down to 100 mW. — [Jacobson-Bell et al. 2025, arXiv:2512.19763](https://arxiv.org/abs/2512.19763)
- [M] Other programmes using the same recipe: Galactic Centre 1–8 GHz; TESS hosts; 97 galaxies; FAST multibeam; LOFAR dual-site. — [Gajjar et al. 2021, arXiv:2104.14148](https://arxiv.org/abs/2104.14148) ([P] confirms the GC strategy paper exists); [Choza et al. 2024 (search link)](https://arxiv.org/search/?query=Choza+97+nearby+galaxies+technosignature&searchtype=all)

#### Drift-rate ranges
- [P] Sheikh et al. 2019 recommend a normalised maximum drift of **200 nHz** (for example 200 Hz/s at 1 GHz) as "a generous, physically motivated guideline". — [arXiv:1910.01148](https://arxiv.org/abs/1910.01148)
- [M] Li et al. 2023: 99% of simulated drifts for about 5,300 known exoplanets lie within ±53 nHz; the bias-corrected population is within about ±0.44 nHz. — [Li et al. 2023 (search link)](https://arxiv.org/search/?query=Li+drift+rate+distribution+technosignature+exoplanets&searchtype=all)
- [P] Fine-resolution products lose sensitivity at higher drift. Ma et al. note that their 4,096-bin snippets limit the maximum drift and that signals can "drift out" of snippets. — [Ma et al. 2023, arXiv:2301.12670](https://arxiv.org/abs/2301.12670)

#### ON/OFF cadence filtering
- [P] ABACAD cadence: three 5-minute ON scans interspersed with OFF scans. — [Enriquez et al. 2017](https://arxiv.org/abs/1709.03491)
- [P†] turboSETI `find_events` filter levels: level 1 applies thresholds only; level 2 keeps hits in at least one ON scan and no OFF scan; level 3 (the default) requires all ON scans and no OFF scan. Zero-drift hits are excluded by default. The pipeline example uses `SNR_cut=10`, `filter_threshold=3`, `number_in_cadence=6`. — [turbo-seti 2.3.2 source](https://pypi.org/project/turbo-seti/2.3.2/)
- [P†] Matching windows: OFF-scan veto within ±|max drift|×(Δt + obs length); ON-scan follow within ±|drift|×2×(obs length) of the extrapolated frequency. — [turbo-seti 2.3.2 source](https://pypi.org/project/turbo-seti/2.3.2/)

#### Coincidence, multibeam and multi-site rejection
- [P] Perytons, millisecond swept-frequency signals mimicking dispersion at Parkes, were traced to a terrestrial source. [M] The source was microwave ovens opened mid-cycle, and the signals appeared in all 13 beams. — [Petroff et al. 2015, arXiv:1504.02165](https://arxiv.org/abs/1504.02165)
- [M] FAST multibeam coincidence; LOFAR dual-site (Ireland and Sweden); Project Phoenix pseudo-interferometry. — [Tao et al. 2022 (search link)](https://arxiv.org/search/?query=Tao+multi-beam+SETI+FAST+exoplanet&searchtype=all); [Johnson et al. 2023 (search link)](https://arxiv.org/search/?query=Johnson+dual-site+technosignature+LOFAR&searchtype=all)
- [P] Sheikh et al. 2021 recommend that "simultaneous observations of potential technosignatures should be conducted at two different observing sites". For Proxima Centauri, MeerKAT and Murriyang share about 1 GHz of frequency overlap and a window of about 4.5 h. — [arXiv:2111.06350](https://arxiv.org/abs/2111.06350)

#### Scintillation-based discrimination
- [P] Brzycki et al. 2023 investigate detecting intensity scintillation from the ionised interstellar medium as a discriminant that does not depend on sky position. [M] The statistics are the standard deviation, the minimum, a KS test against an exponential distribution, and an ACF timescale, compared with NE2001 predictions. — [arXiv:2307.08793](https://arxiv.org/abs/2307.08793)
- [P†] `blscint` 0.0.1 (January 2022) is on PyPI. — [PyPI](https://pypi.org/project/blscint/)

#### Machine-learning approaches
- [P] Harp et al. 2019 used classical machine vision and deep learning on archived and simulated spectrograms; Brzycki et al. 2020 localised narrowband signals in synthetic spectrograms. — [arXiv:1902.02426](https://arxiv.org/abs/1902.02426); [arXiv:2006.04362](https://arxiv.org/abs/2006.04362)
- [M] Pinchuk & Margot 2022: a CNN direction-of-origin filter. — [search link](https://arxiv.org/search/?query=Pinchuk+Margot+machine+learning+direction-of-origin&searchtype=all)
- [P] Ma et al. 2023: a β-VAE plus random forest on 820 targets and more than 480 h of GBT data. **(confirmed from full text)** There were 115 million snippet windows (57 million unique). The model returned 2,917,789 signals of interest, which were reduced to 20,515 by flagging RFI-dense bins; visual inspection then gave 8 signals of interest from five stars (HIP 13402, 54677, 62207, 56802, 118212; 30–90 ly). Training used 120,000-sample sets with injected signals. — [arXiv:2301.12670](https://arxiv.org/abs/2301.12670)
- [P] Jacobson-Bell et al. 2024 (GLOBULAR) cluster turboSETI hits with HDBSCAN. Benchmarked on the 97-galaxy L-band search, they report false-positive reductions of 93.1% (hits) and 99.3% (events). — [arXiv:2411.16556](https://arxiv.org/abs/2411.16556)
- [P†] `setigen` 2.7.0 (August 2024) injects synthetic drifting signals and cadences, including raw voltages. — [PyPI setigen](https://pypi.org/project/setigen/); [Brzycki et al. 2022, AJ 163, 222](https://doi.org/10.3847/1538-3881/ac5e3d)

#### BLC1
- [P] Sheikh et al. 2021: BLC1 "is not an extraterrestrial technosignature, but rather an electronically-drifting intermodulation product of local, time-varying interferers aligned with the observing cadence"; dozens of look-alikes were found at frequencies harmonically related to common clock oscillators. The ten-step framework is listed in Section 10. — [arXiv:2111.06350](https://arxiv.org/abs/2111.06350)

#### Software (PyPI, checked 14 November 2026)
- [P†] On PyPI: `turbo-seti` 2.3.2, `blimpy` 2.1.4, `setigen` 2.7.0, `blscint` 0.0.1, `sigpyproc` 2.0.0 (October 2026), `your` 0.6.7, `baseband` 4.3.0, `iqrm` 0.2.0. Not on PyPI under these names: `hyperseti`, `seticore`, `bliss-seti`. — [PyPI turbo-seti](https://pypi.org/project/turbo-seti/); [PyPI setigen](https://pypi.org/project/setigen/)

### Inferences
- (computed) For δν = 2.794 Hz, δt = 18.25 s and 16 integrations: drift resolution ≈ 0.0102 Hz/s; ±4 Hz/s is about 784 trials. The one-channel-per-sample rate δν/δt ≈ 0.153 Hz/s. Above it, only a fraction of about 0.153/|drift| of the power lands on the path (about 15% at 1 Hz/s, about 4% at 4 Hz/s). This is consistent with Margot's "about 6%" recovery for the 51-spectrum product. Report sensitivity as a function of drift.
- (computed) About 2.9×10⁸ channels × 784 drifts ≈ 2.2×10¹¹ trials per 800 MHz scan. The Gaussian tail at 10σ is 7.6×10⁻²⁴, so the expected number of thermal false alarms is about 2×10⁻¹². Thresholds of 10–25σ are therefore set by the human capacity to vet interference, not by noise.
- (computed) Earth's rotation contributes up to about 0.11 nHz. A ±4 Hz/s window is about ±3.6 nHz at 1.1 GHz and about ±0.5 nHz at 8 GHz, far short of 53 nHz or 200 nHz.
- Built-in failure modes of the cadence filter: (a) a signal just below threshold in OFF scans passes, so search OFF scans at a lower threshold (Sheikh step 2); (b) level 3 rejects transmitters with less than 100% duty cycle over 30 minutes; (c) interference whose on/off pattern matches the cadence passes (BLC1); (d) the zero-drift exclusion removes acceleration-compensated transmitters.
- Recommended null test: swap the ON and OFF labels and count the events that pass. This measures chance cadence coincidence directly.

### Gaps
- No paper gives a calibrated false-alarm probability for the cadence filter itself; only funnels exist.
- COSMIC (VLA) and MeerKAT BLUSE thresholds and results, and Brzycki's 2024 Galactic Centre scintillation search, were not checked.
- The current status of `hyperseti`, `seticore` and `bliss` was not checked; they are not on PyPI under those names.

## 3. Other signal morphologies: wideband pulses, artificial dispersion, periodic pulse trains, modulated carriers, spread spectrum; cyclostationary and modulation-recognition techniques

### Takeaway
Beyond drifting carriers, published searches cover dispersed broadband pulses (including negative dispersion), periodic pulse trains (FFA) and autocorrelation tests for repeated wideband waveforms. All returned null results. Cyclostationary analysis and automatic modulation classification need complex voltages, and in practice they classify human interference. They are best used to reject interference and to characterise a candidate after detection.

### Cited Findings

#### Broadband pulses and artificial dispersion
- [P] Gajjar et al. 2022: broadband pulsed beacons are energetically efficient over long operating times. The survey covered 1883 stars in 233 h at 4–8 GHz and searched three classes of artificially (negatively) dispersed signals, using a CNN classifier on GPUs. Result: at most about 1 in 1000 stars have transmitter power densities of at least about 10⁵ W/Hz repeating at intervals of 500 s or less. — [arXiv:2205.02964](https://arxiv.org/abs/2205.02964)
- [P] The FDMT algorithm gives optimal detection of pulses with unknown dispersion measure (DM). FETCH classifies transient candidates with deep learning. [M] The cold-plasma delay is Δt ≈ 4.15 ms·DM·(ν_lo⁻² − ν_hi⁻²), with ν in GHz. — [Zackay & Ofek 2017, arXiv:1411.5373](https://arxiv.org/abs/1411.5373); [Agarwal et al. 2020, arXiv:1902.06343](https://arxiv.org/abs/1902.06343)
- [P] A CNN found 72 new FRB 121102 pulses in Breakthrough Listen C-band data from 26 August 2017, where 21 had been found previously. — [Zhang et al. 2018, arXiv:1809.03043](https://arxiv.org/abs/1809.03043)
- [P] **Speculative:** FRBs could be beams powering extragalactic light sails. — [Lingam & Loeb 2017, arXiv:1701.01109](https://arxiv.org/abs/1701.01109)

#### Periodic pulse trains
- [P] BLIPSS (open-source, CPU, per-channel FFA): 4.5 h of 4–8 GHz GBT data on the central 6′ of the Galaxy, searching kHz-wide signals with periods of 11–100 s and duty cycles of 10–50%. Non-detection; fewer than 1 in about 600,000 Galactic Centre stars host a transmitter above a 7σ EIRP of about 2×10¹⁸ W at δ ≈ 10%. — [arXiv:2305.18527](https://arxiv.org/abs/2305.18527)
- [P] An FFA search "closely approximates" the optimal periodic-signal search, unlike FFT with incoherent harmonic summing. [P†] `riptide-ffa` 0.2.7 (September 2026) is on PyPI. — [Morello et al. 2020, arXiv:2004.03701](https://arxiv.org/abs/2004.03701); [PyPI](https://pypi.org/project/riptide-ffa/)
- [P] **Speculative:** civilisations modulating pulsar signals; pulsar positioning as a "quest for evidence". — [Chennamangalam et al. 2015, arXiv:1311.4608](https://arxiv.org/abs/1311.4608); [Vidal 2019, arXiv:1704.03316](https://arxiv.org/abs/1704.03316)

#### Wideband, spread-spectrum and information-bearing signals
- [P] Messerschmitt advocates spread spectrum, which is noise-like and has much higher bandwidth than necessary, for its immunity to interference. — [Messerschmitt 2012, arXiv:1111.0547](https://arxiv.org/abs/1111.0547)
- [P] Harp et al. propose information-bearing beacons that remain discoverable, motivated by the efficiency of narrowband searches. [M] The mechanism is a delayed repeat of a wideband signal, detected as an autocorrelation peak. — [Harp et al. 2012, arXiv:1211.6470](https://arxiv.org/abs/1211.6470)
- [M] Related proposals: antipodal-signalling autocorrelation (Morrison 2012) and the KLT (Maccone 2010); identifiers not verified.
- [M] The generalised spectral kurtosis estimator has expectation 1 for Gaussian noise, falls below 1 for continuous-wave signals and rises above 1 for pulsed signals. — [Nita & Gary 2010 (search link)](https://arxiv.org/search/?query=Nita+Gary+generalized+spectral+kurtosis&searchtype=all)

#### Cyclostationary analysis and modulation classification
- [M] Spectral correlation S(f; α): stationary noise contributes only at α = 0, so features at α ≠ 0 reveal symbol, chip and carrier rates. Cycle-frequency resolution ≈ 1/T_obs. Estimators are the FFT accumulation method and the strip spectral correlation analyser. — [Gardner 1991](https://doi.org/10.1109/79.81007)
- [M] Fourth-order cumulants separate BPSK, PSK and QAM. — [Swami & Sadler 2000](https://doi.org/10.1109/26.837045)
- [P] Over-the-air deep-learning classification was compared against a higher-order-moment baseline with boosted trees across channel impairments. — [O'Shea, Roy & Clancy 2018, arXiv:1712.04578](https://arxiv.org/abs/1712.04578)
- [P†] `torchsig` 2.2.0 (August 2026), `scikit-commpy` 0.8.0 and `sionna` 2.2.0 are on PyPI; there is no `cyclostationary` package. — [PyPI torchsig](https://pypi.org/project/torchsig/)

### Inferences
- Public fine-resolution Breakthrough Listen products are total power at about 18 s sampling, so only slow amplitude modulation can be tested on them. Symbol and carrier analysis needs raw voltages.
- Order of blind characterisation steps for a voltage snippet, cheapest first: spectral kurtosis; spectra of |x|², x² and x⁴; voltage autocorrelation; spectral coherence on the strongest hits; a trained classifier. Positive outputs are also interference fingerprints (for example the 1.023 MHz GPS C/A chip rate).
- Low-duty-cycle beacons fail the "all ON scans" rule, so each morphology needs its own ON/OFF logic (consistent period and phase across scans).

### Gaps
- No published search was found that applies cyclostationary or modulation-recognition methods to a real technosignature candidate.
- Sensitivity of Allen Telescope Array autocorrelation searches and 2025 real-time ML pulse searches at the ATA were not checked.

## 4. Optical and infrared approaches: nanosecond pulses, laser lines, anomalous transits and megastructures, infrared excess / Dyson spheres, vanishing sources and the SETI ellipsoid

### Takeaway
Optical and infrared technosignature searches have produced several high-profile "anomalies". Each was later given, or is being given, a natural or instrumental explanation:
- Boyajian's Star dips were attributed to non-grey extinction by ordinary dust.
- The Project Hephaistos Dyson-sphere candidates turned out to be blends with background dust-obscured galaxies; JWST confirmed this for two candidates in July 2026.
- The VASCO plate transients remain actively contested (emulsion or copying defects and schedule artefacts versus real sub-second flashes), with a burst of claims and rebuttals in 2025–2026.

Laser-line and nanosecond-pulse searches are clean null results with injection-calibrated thresholds. Only laser-line searches in archival spectra can be run on public archives.

### Cited Findings

#### Nanosecond optical pulse searches (established; not archive-searchable)
- [P] Air-Cherenkov telescopes can monitor 10–10² stars simultaneously for nanosecond pulses of about 6th magnitude or brighter, collecting tens of photons from a pulse of isotropic-equivalent power of tens of solar luminosities at 100 pc. — [Eichler & Beskin 2001, astro-ph/0111081](https://arxiv.org/abs/astro-ph/0111081)
- [P] Whipple 10 m: cosmic-ray (Cherenkov) background was discriminated by image shape, since a point-source optical pulse produces a different camera image from an air shower. — [Holder et al. 2005, astro-ph/0506758](https://arxiv.org/abs/astro-ph/0506758)
- [P] VERITAS/Breakthrough Listen: 30 h of dedicated observations of 136 targets plus 249 archival observations of 140 targets (since 2012), with no technosignature. Analysis efficiency was validated by detecting the CALIOP space-borne pulsed laser on CALIPSO. — [Acharyya et al. 2023, arXiv:2306.17680](https://arxiv.org/abs/2306.17680)
- [P] PANOSETI aims at nanosecond-to-second optical transients over about 2,500 deg², using two assemblies of telescopes for coincidence rejection, with clocks synchronised by White Rabbit at the nanosecond level. **2026:** the first light of "Dark100" (0.5 m Fresnel-lens telescopes with 1024-pixel SiPM cameras at Palomar) is reported as doubling as an ultra-high-energy gamma-ray Cherenkov array. — [Wright et al. 2018, arXiv:1808.05772](https://arxiv.org/abs/1808.05772); [Maire et al. 2022, arXiv:2210.01356](https://arxiv.org/abs/2210.01356); [Korzoun et al. 2026, arXiv:2610.41775](https://arxiv.org/abs/2610.41775)
- [M] The Harvard optical SETI programmes (Howard et al. 2004, 2007) used pairs of hybrid avalanche photodiodes in coincidence, with nanosecond thresholds of a few photons per detector; they reported occasional unexplained single coincidences that did not repeat. Identifiers not verified.

#### Continuous laser-line searches in high-resolution spectra (established; archive-searchable)
- [P] Tellis & Marcy 2017 searched Keck HIRES spectra of 5600 FGKM stars (3640–7890 Å, within 0.5″ of the star) for emission lines that are both inconsistent with a natural origin and spatially unresolved. Thresholds were set by injecting synthetic laser lines and analysing them blind. Laser power limits are 3 kW–13 MW (for a diffraction-limited 10 m transmitter). No detections; models in which more than 0.1% of warm Earth-size planets beam lasers at us are excluded. — [arXiv:1704.02535](https://arxiv.org/abs/1704.02535)
- [P] Zuckerman et al. 2023 built an APF/Levy pipeline: SpecMatch-Emp residuals against best-matching library spectra, rejection of events that do not persist between consecutive same-night spectra, and analysis of 1983 spectra of 388 stars. Sensitivity is an 84 kW laser at the median distance of 78.5 ly. No convincing candidates. — [arXiv:2301.06971](https://arxiv.org/abs/2301.06971)
- [P] Hippke argues that the "magic colours" for optical SETI are Nd:YAG 1064 nm and its second harmonic at 532.1 nm, then 393.8, 656.5 and 589.1 nm. Several of these coincide with strong stellar absorption features (Ca K, Hα, Na D2) and so need special handling. — [Hippke 2018, arXiv:1804.01249](https://arxiv.org/abs/1804.01249)
- [P] **2025 development (proposal):** a spectroscopic laser search with ESO's future highly multiplexed spectrographs in the 2040s. — [Ivanov 2025, arXiv:2512.18903](https://arxiv.org/abs/2512.18903)
- [M] The main false positives in laser-line searches are cosmic-ray hits (single-pixel, not reproduced in a second exposure), sky emission lines (OH, Na, Hg street lamps), fibre and echelle ghosts, and the residuals of telluric correction. Persistence across exposures and the spatial profile along the slit are the key discriminants.

#### Anomalous transit shapes and megastructures (Arnold 2005; Wright 2016; Boyajian's Star)
- [P] Arnold 2005 shows that triangular, two-screen and louvred six-screen objects give transit light curves distinguishable from those of spherical or oblate planets (with possible ambiguity against ringed planets), and that multiple artificial objects would be "easily distinguishable". A companion paper argues that artificial transits could be used as attention-getting signals. — [Arnold 2005, ApJ; astro-ph/0503580](https://arxiv.org/abs/astro-ph/0503580); [Arnold 2005, astro-ph/0509431](https://arxiv.org/abs/astro-ph/0509431)
- [P] Wright et al. 2016 enumerate **ten** ways a transiting megastructure's silhouette, orbit and transmission would differ from a planet's, and the natural mimics of each. They note that KIC 12557548 and CoRoT-29 show such anomalies but have natural explanations. They introduce a **normalised information-content statistic M** for a discrete series of bounded measurements (for example transit depths), intended to separate constant sources, beacons and information-rich signals. — [Wright et al. 2016, ApJ 816, 17; arXiv:1510.04606](https://arxiv.org/abs/1510.04606)
- [P] Boyajian et al. 2016 (Planet Hunters X): KIC 8462852 showed aperiodic dips of up to about 20% lasting 5–80 days. The star is a typical F3 V with no significant infrared excess. Instrumental artefacts were excluded. The favoured explanation was a family of exocomet fragments. — [arXiv:1509.03622](https://arxiv.org/abs/1509.03622)
- [P] Long-term fading: 0.341 ± 0.041%/yr over the first about 1000 days of Kepler full-frame images, then a drop of more than 2% (Montet & Simon). ASAS-SN measured 6.3 ± 1.4 mmag/yr from 2015, significant at more than 99.4% against about 1000 comparison stars (Simon et al.). Warm Spitzer showed no 3.6 µm excess and only 0.43 ± 0.18 mJy at 4.5 µm. — [Montet & Simon 2016, arXiv:1608.01316](https://arxiv.org/abs/1608.01316); [Simon et al. 2018, arXiv:1708.07822](https://arxiv.org/abs/1708.07822); [Marengo et al. 2015, arXiv:1511.07908](https://arxiv.org/abs/1511.07908)
- [P] **Resolution:** the 2017 dips ("Elsie", "Celeste", "Skara Brae", "Angkor"; 1–2.5%) showed no spectral or polarisation changes and showed non-grey (chromatic) extinction, which is "inconsistent with dip models that invoke optically thick material" (an opaque megastructure would be grey) and in line with an occulter of ordinary dust. Swift/Spitzer fading was likewise wavelength-dependent (22.1 ± 9.7 mmag/yr in the UV versus 5.0 ± 1.2 mmag/yr at 3.6/4.5 µm). — [Boyajian et al. 2018, arXiv:1801.00732](https://arxiv.org/abs/1801.00732); [Meng et al. 2017, arXiv:1708.07556](https://arxiv.org/abs/1708.07556); [Bodman et al. 2019, arXiv:1806.08842](https://arxiv.org/abs/1806.08842)
- [P] A faint star 2″ away is not a bound companion: it has a relative proper motion of 23.9 ± 2.6 mas/yr. — [Clemens et al. 2018, arXiv:1803.03299](https://arxiv.org/abs/1803.03299)
- [P] Nonparametric detectors recover Boyajian's Star. The "Weird Detector" (phase-dispersion minimisation on 161,786 Kepler sources) recovered it along with a comet crossing the field, and found 18 new heartbeat or planet-like signals plus single and double transits. Density-based clustering by Giles & Walkowicz also flags it. — [Wheeler & Kipping 2019, arXiv:1903.09213](https://arxiv.org/abs/1903.09213); [Giles & Walkowicz 2018, arXiv:1812.07156](https://arxiv.org/abs/1812.07156)
- [P] **Speculative proposal:** planetary mean-motion resonances could encode sequences (primes 2–11, Fibonacci 1–13, the Lazy Caterer sequence 1–16). N-body tests show these remain stable for 10 Gyr plus post-main-sequence evolution, whereas the consecutive-integer system destabilised. — [Clement et al. 2022, arXiv:2204.14259](https://arxiv.org/abs/2204.14259)

#### Infrared excess and Dyson spheres (Ĝ; Project Hephaistos I–IV; 2026 resolution)
- [P] Ĝ: a formalism for turning mid-infrared photometry into upper limits on energy supply, with dust named as the dominant false positive. A WISE survey of about 1×10⁵ galaxies found 563 red extended sources (|b| ≥ 10°, visually vetted for artefacts) and **no** galaxy reprocessing more than 85% of its starlight into the mid-infrared. — [Wright et al. 2014, arXiv:1408.1134](https://arxiv.org/abs/1408.1134); [Griffith et al. 2015, arXiv:1504.03418](https://arxiv.org/abs/1504.03418)
- [P] Hephaistos I (Gaia DR2 + AllWISE; Dyson spheres absorbing 1–90% of the light, with waste heat at 100–1000 K): for about 2.7×10⁵ stars within 100 pc, fewer than about 2×10⁻⁵ could host about 300 K spheres at 90% completion. For about 2.9×10⁸ stars within 5 kpc, the limit is ≤ 8×10⁻⁴. Limits weaken for partial spheres because of confusion with natural mid-infrared sources. — [Suazo et al. 2022, arXiv:2201.11123](https://arxiv.org/abs/2201.11123)
- [P] Hephaistos II (Gaia DR3 + 2MASS + WISE): a multi-filter pipeline on about 5 million objects, including a CNN to flag WISE confusion. It left **7 candidates, all M dwarfs**. — [Suazo et al. 2024, MNRAS; arXiv:2405.02927](https://arxiv.org/abs/2405.02927)
- [P] Ren, Garrett & Siemion (May 2024) found radio sources offset by about 4.9″, 0.4″ and 5.0″ from candidates A, B and G. They argued that dust-obscured galaxies along the line of sight dominate the W3/W4 fluxes, and that with a Hot DOG sky density of 9×10⁻⁶ arcsec⁻², chance blends can probably explain all 7 of 5 million. — [arXiv:2405.14921](https://arxiv.org/abs/2405.14921)
- [P] Contardo & Hogg 2024: a data-driven contextual anomaly search on 4,898,812 main-sequence FGK stars (G < 16; Gaia DR3 + 2MASS + unWISE). An ensemble of ML models predicts W1/W2 from optical and near-infrared inputs; after cuts, 53 objects remain (a rate of 1.1×10⁻⁵), with fractional luminosities of 0.005–0.1, consistent with extreme debris disks. — [arXiv:2403.18941](https://arxiv.org/abs/2403.18941)
- [P] **2026 developments, July:**
  - **Hephaistos III** finds the 7 + 3 candidates' stellar parameters normal and sees no signs of youth; it says the excess is "still" unexplained and needs JWST or ALMA.
  - **Hephaistos IV** (JWST/MIRI imaging and spectroscopy of candidates D and E) attributes the excess to **background galaxies within about 1″**: a Hot DOG at z ≈ 0.9 for D and a dusty starburst at z ≈ 0.4 for E.
  - **Ren et al. 2026** (archival diagnostics of nine candidates) compared Gaia positions propagated to the AllWISE epoch with the mid-infrared centroids. Significant contamination is found in B (a radio counterpart with α = 0.63 ± 0.11) and C (a near-infrared companion at 3.75″); A is suggestive (α = 0.40 ± 0.35); E, F, H and J are marginal; D and I showed no obvious contamination archivally (but D is now resolved by JWST).

  — [Korn et al. 2026, arXiv:2607.25701](https://arxiv.org/abs/2607.25701); [Zackrisson et al. 2026, arXiv:2607.09460](https://arxiv.org/abs/2607.09460); [Ren et al. 2026, arXiv:2607.03619](https://arxiv.org/abs/2607.03619)

#### Vanishing-source and plate-transient searches (VASCO) — contested
- [P] VASCO I: 600 million USNO-B1.0 objects were cross-matched with Pan-STARRS DR1, giving about 150,000 objects with no counterpart within 30″. These are redder and have higher proper motions than typical objects. Visual inspection of about 24,000 left about 100 single-epoch red point sources (flares, distant supernovae and similar). — [Villarroel et al. 2020, AJ 159, 8; arXiv:1911.05068](https://arxiv.org/abs/1911.05068)
- [P] A Virtual Observatory search of POSS-I red plates against Pan-STARRS DR2 and Gaia EDR3 found 298,165 POSS-I-only sources. Of these, 288,770 were matched within 5″ in other archives (mostly infrared), 189 were asteroids, 35 variables, 3,592 artefacts (by comparison with a second, SuperCOSMOS, digitisation) and 180 high-proper-motion stars, leaving **5,399 unidentified transients**. — [Solano et al. 2022, arXiv:2206.00907](https://arxiv.org/abs/2206.00907)
- [P] Method: in an image with nine "simultaneous transients", four or five aligned within a 10′×10′ box would indicate reflective objects near geosynchronous orbit at 2.5–3.9σ. The proposal is to search pre-1957 plates. — [Villarroel et al. 2022, arXiv:2110.15217](https://arxiv.org/abs/2110.15217); [Villarroel et al. 2022, arXiv:2204.06091](https://arxiv.org/abs/2204.06091)
- [P] A triple transient (about 15th mag, within 10″, POSS-I, 19 July 1952) is absent in GTC imaging to r ≈ 25.5, implying a fade of more than 10 mag within an hour. — [Solano et al. 2023, arXiv:2310.09035](https://arxiv.org/abs/2310.09035)
- [P] **Critique:** Hambly & Blair 2024 built a random-forest morphology classifier on SuperCOSMOS catalogues and found the putative transients "likely to be spurious artefacts of the photographic emulsion", possibly arising in the glass-copy process. — [arXiv:2402.00497](https://arxiv.org/abs/2402.00497)
- [P] **Response (2025):** Villarroel, Solano & Marcy argue that the narrower, rounder profiles are expected for sub-second flashes on long exposures, where stars are blurred by seeing and tracking. — [arXiv:2507.15896](https://arxiv.org/abs/2507.15896)
- [P] **Critique (January 2026):** Watters et al. find the following.
  - The uniform-random background assumed in the Earth-shadow test is false for the datasets they examined; after deriving the null distribution, there is no significant shadow deficit.
  - The nuclear-test correlation becomes insignificant once normalised by observation days, and it "is almost completely determined by the observation schedule of the Palomar telescope".
  - The datasets contain catalogue stars, scan artefacts and plate defects; feature density rises toward plate edges and corners; there are north–south empty strips spanning plates and geometric clusters.
  - The argument used is circular.

  — [Watters et al. 2026, arXiv:2601.21946](https://arxiv.org/abs/2601.21946)
- [P] **Pro-transient papers (2026):**
  - Doherty (March 2026) reproduces relative risk 1.45 (p = 0.011), incidence rate ratio (IRR) 1.80 and a permutation p of 0.006 for the nuclear correlation, and a shadow fraction of 0.46% against about 1.4% expected, using the original authors' dataset.
  - Bruehl et al. (April 2026) trained an ML real-versus-defect classifier on 250 expert-labelled pairs (AUC 0.81; sensitivity and specificity 0.71), applied it to 107,875 transients, and report a nuclear-window p of 0.024 and a shadow deficit with p < 0.0001.
  - Doherty (May 2026) reports 7 plates with 5–8-source alignments narrower than 3″ (p < 0.03, 10,000 Monte Carlo runs) among 36,215 candidates with probability ≥ 0.5.
  - A Hamburg-plate (APPLAUSE) search reports similar narrow-FWHM transients.
  - A microscopy study (August 2026) of the 1950 nine-transient and 1952 triple-transient events concludes that microscopy "cannot by itself separate a short-duration flash from an intrinsic emulsion defect".

  — [arXiv:2604.00056](https://arxiv.org/abs/2604.00056); [arXiv:2604.18799](https://arxiv.org/abs/2604.18799); [arXiv:2605.01190](https://arxiv.org/abs/2605.01190); [Busko 2026, arXiv:2603.20407](https://arxiv.org/abs/2603.20407); [Yang & Villarroel 2026, arXiv:2608.26371](https://arxiv.org/abs/2608.26371)
- [P] **Fringe-leaning (flag as such):** single-author preprints claim that geomagnetic-storm suppression of transient rates "rul[es] out emulsion defects" (Z = −3.391, p = 0.0007) and propose trapped icy dust at L ≈ 6.6. Further preprints model geosynchronous reflectors (September 2026). — [Cann 2026, arXiv:2604.04950](https://arxiv.org/abs/2604.04950); [Cann 2026, arXiv:2604.06234](https://arxiv.org/abs/2604.06234); [Villarroel et al. 2026, arXiv:2609.05105](https://arxiv.org/abs/2609.05105); [Bruehl et al. 2026, arXiv:2609.09461](https://arxiv.org/abs/2609.09461)
- [P] Modern analogue confounder: W-FAST (25 Hz) shows that 0.1–0.3 s flares of mag 9–11 are mostly glints off geosynchronous and graveyard-orbit satellites. The rate is about 30–40 events per day per deg² brighter than about 11 mag; diluted into a 30 s exposure they look like 14–16 mag point sources. — [Nir et al. 2021, arXiv:2011.03497](https://arxiv.org/abs/2011.03497)

#### The SETI ellipsoid (synchronisation)
- [P] Seto proposed a coordinated signalling scheme keyed to a conspicuous event such as a neutron-star merger, so that the receiver need examine only a tiny fraction of the sky as a function of time. **2025:** a hybrid scheme uses the Galactic Centre as the spatial reference and an extragalactic burst as the time marker, reducing the area to survey by more than two orders of magnitude. — [Seto 2019, arXiv:1904.00536](https://arxiv.org/abs/1904.00536); [Seto 2025, arXiv:2509.20718](https://arxiv.org/abs/2509.20718)
- [P] Davenport et al. 2022 used the Gaia Catalogue of Nearby Stars with SN 1987A and 278 classical novae as reference events. Fewer than 8% of stars within 100 pc are inside the SN 1987A ellipsoid; on average 734 stars per year cross it, about 10% of them with distance uncertainties better than 0.1 ly. — [arXiv:2206.04092](https://arxiv.org/abs/2206.04092)
- [P] Cabrales et al. 2024: 32 SN 1987A ellipsoid targets in the TESS continuous viewing zone (5% of TESS data, first three years) with Gaia EDR3 uncertainties better than 0.5 ly. No anomalous signatures during crossing. — [arXiv:2402.11037](https://arxiv.org/abs/2402.11037)
- [P] SN 2023ixf: more than 100 stars within 100 pc are near its ellipsoid, within about 3° of the supernova; radio monitoring with the ATA and GBT began. — [Davenport et al. 2023, arXiv:2306.03118](https://arxiv.org/abs/2306.03118)

### Inferences
- Common pattern: every IR-excess and anomaly candidate failed on **blending and confusion**. A tool should automatically check for background-source blending before flagging a source: the offset between the Gaia position (propagated to the IR epoch) and the IR centroid; radio counterparts and their spectral index; high-resolution imaging; and the expected chance-blend count, density × search area × N_stars. (computed) 9×10⁻⁶ arcsec⁻² × π(6″)² ≈ 1×10⁻³ blends per star, roughly 5,000 blends among 5×10⁶ stars before colour cuts, so finding 7 candidates is unsurprising.
- Chromaticity is a cheap megastructure discriminant: an opaque structure is grey, dust is not. Multi-band dip depths should be part of transit-anomaly vetting.
- Treat VASCO-type plate transients as **contested**. Their significance depends on the null model for spatial and temporal non-uniformity (plate edges, observing schedule). This is exactly the Section 8 issue: test statistics need nulls built from the real observation schedule and plate geometry, by permutation over actual observing dates and positions, not uniform-random assumptions.
- Synchronisation targeting is cheap to implement from Gaia: compute ellipsoid crossing times for each star, with distance uncertainties, and attach them as a prior or time window to any time-domain search.

### Gaps
- Harvard OSETI and NIROSETI results and thresholds were not verified (not found in arXiv title searches).
- No 2025–2026 Boyajian's Star papers were found in arXiv searches.
- Whether Hephaistos candidates other than D and E have been resolved by JWST or ALMA after July 2026 is unknown.
- The peer-reviewed version of Bruehl & Villarroel 2025 (the original nuclear-test correlation) was not located; it is known only through 2026 citations.

## 5. Proposals that natural variable sources could be modulated to carry information, and follow-up

### Takeaway
Learned et al. (2008) proposed that civilisations might "tickle" Cepheids to phase-modulate their pulsation as galactic beacons. Hippke, Learned et al. searched Kepler RR Lyrae and Cepheid photometry: they reported period ratios near √2 (2014) and a prime-number autocorrelation oddity in KIC 5520878. They themselves concluded the latter is natural (two simultaneous pulsation modes). Asteroseismology now explains the common unusual period ratios in RRc stars as non-radial modes. The lasting lesson is methodological: O–C (phase) time series of regular variables are a legitimate, automatable search space, and pattern hunting in them is prone to pareidolia.

### Cited Findings
- [P] Learned, Kudritzki, Pakvasa & Zee propose "the Cepheid Galactic Internet": advanced civilisations could use Cepheids as beacons, and Cepheids and other regular variables should be searched "for signs of phase modulation (in the regime of short pulse duration) and patterns". The rationale is that Cepheids are seen out to Virgo and are monitored by any society that uses them as distance markers. — [Learned et al. 2008, arXiv:0809.0339](https://arxiv.org/abs/0809.0339)
- [P] Hippke, Learned & Zee report a few Cepheids and RR Lyrae stars with main/secondary period ratios P1/P2 very close to √2, with deviations in discrete steps of about 0.000390(4) (multiples of about 0.06%), "indicating the possible presence of a sort of fine structure". — [Hippke, Learned & Zee 2014/2015, arXiv:1407.3442](https://arxiv.org/abs/1407.3442)
- [P] KIC 5520878 (an RRc in Kepler): two regimes of long and short pulse durations. The period-length sequence is strongly autocorrelated, with correlation coefficients at prime-number lags "significantly higher (p = 99.8%)". However, "the prime number oddity originates from two simultaneous pulsation periods and is likely of natural origin". The authors encourage automated testing of archival photometry for modulated stars as a by-product of standard analyses. — [Hippke et al. 2015, ApJ 798, 42; arXiv:1409.1265](https://arxiv.org/abs/1409.1265)
- [P] Natural explanation of "odd" period ratios: low-amplitude signals at period ratios of about 0.61 (range 0.60–0.64) relative to the first overtone are common in RRc and RRd stars and in first-overtone Cepheids. They are modelled as harmonics of non-radial modes of degree ℓ = 8 or 9 (7–9 for Cepheids). An OGLE census covered 11,415 RRc and 148 RRd stars, and there is a second group at about 0.68. — [Netzel & Smolec 2019, arXiv:1908.07251](https://arxiv.org/abs/1908.07251); [Netzel & Smolec 2022, arXiv:2206.12853](https://arxiv.org/abs/2206.12853); [Netzel et al. 2024, arXiv:2404.14373](https://arxiv.org/abs/2404.14373)
- [P] Related engineered-natural-system proposals: resonant planetary systems encoding integer sequences, and Chennamangalam et al. on modulating pulsars. Both are speculative. — [Clement et al. 2022, arXiv:2204.14259](https://arxiv.org/abs/2204.14259); [Chennamangalam et al. 2015, arXiv:1311.4608](https://arxiv.org/abs/1311.4608)
- [P] The information-content statistic M of Wright et al. 2016 (Section 4) applies to any bounded discrete series, including pulse-to-pulse periods or transit depths. — [arXiv:1510.04606](https://arxiv.org/abs/1510.04606)

### Inferences
- Pipeline recipe for "modulated variable" tests:
  1. Fit a linear (or quadratic) ephemeris.
  2. Build the cycle-by-cycle O–C or period series.
  3. Remove known physics (multi-mode beating, the Blazhko effect, light-travel time in binaries) by fitting a multi-frequency model first.
  4. Test the residual series against surrogates that preserve its power spectrum (Section 8).
  5. Only then compute pattern statistics (prime lags, ratios to irrational constants).

  The KIC 5520878 case shows that two beating modes alone create apparently "structured" lag correlations.
- (computed) The look-elsewhere effect is acute here. There are about 25 primes below 100, so testing "prime lags" against "non-prime lags" is one of many partitions a searcher could have chosen; ratios close to √2, φ, π or e among hundreds of stars are expected by chance at the 0.1% level. Pre-register the constants and statistics (Section 8).

### Gaps
- I found no independent follow-up paper that specifically re-examined the √2 ratios of Hippke et al. 2014. The non-radial-mode literature explains the 0.61 and 0.68 groups but I did not confirm that it covers the √2 stars.
- No 2025–2026 papers on modulated Cepheids or RR Lyrae as beacons were found.

## 6. Information-theoretic and "encoded data" tests

### Takeaway
No single statistic separates "message" from "nature". The tests that discriminate are those comparing a sequence against null models that match its marginal distribution and power spectrum: entropy-rate estimates (block or conditional entropies, Lempel–Ziv, compression), ordinal and permutation entropy, statistical complexity, the shape of mutual-information decay, and recurrence structure. Zipf slopes and higher-order entropy slopes (McCowan, Hanser and Doyle; Doyle et al. 2011) and Elliott's language-structure work are proposed SETI tools. They are not validated on astrophysical data and need long symbol sequences.

There is a hard limit: an optimally compressed or encrypted message is indistinguishable from noise, and an information-optimal electromagnetic transmission looks thermal (Lachmann, Newman and Moore). So these tests can only catch deliberately redundant, beacon-like signals. The beacon literature (Benford, Messerschmitt, Hippke) predicts sparse, low-duty-cycle, narrowband or pulsed signals, not dense coded streams.

### Cited Findings

#### The fundamental limit
- [P] Lachmann, Newman & Moore: "a message transmitted with optimal efficiency over a channel of limited bandwidth is indistinguishable from random noise to a receiver who is unfamiliar with the language". For electromagnetic transmission, "the most information-efficient format for a given message is indistinguishable from black-body radiation", at a temperature set by the transmitted energy. If information is not encoded in direction, the optimal format has a **one-dimensional** black-body spectrum, which is distinguishable from the three-dimensional case. — [Lachmann, Newman & Moore, arXiv:cond-mat/9907500](https://arxiv.org/abs/cond-mat/9907500); [M] published in Am. J. Phys. 72, 1290 (2004)
- [P] Messerschmitt (end-to-end design for power efficiency): the fundamental limit is set by noise alone and is unaffected by carrier frequency, dispersion, scattering or motion. Power-efficient signals have "wide bandwidth relative to the information rate and a sparse distribution of energy in both time and frequency". Minimum-power beacons need not be wideband but concentrate energy sparsely in time. Discovery "most closely resembles approaches that have been employed in optical SETI". — [Messerschmitt 2013, arXiv:1305.4684](https://arxiv.org/abs/1305.4684)

#### What a deliberately detectable message or beacon looks like
- [P] Benford, Benford & Benford: cost-optimal beacons split capital cost equally between antenna gain and radiated power. Cost, spectral lines near 1 GHz and scintillation favour frequencies well above the water hole, and the transmission strategy is "a rapid scan of the galactic plane". To the receiver this means intermittent, brief, possibly wideband pulses; searches need long dwell times or repeated revisits. — [arXiv:0810.3964](https://arxiv.org/abs/0810.3964); [arXiv:0810.3966](https://arxiv.org/abs/0810.3966)
- [P] Benford 2010 asks how beacons could be distinguished from transient pulsars, pointing to bandwidth, pulse width and frequency. — [arXiv:1003.5938](https://arxiv.org/abs/1003.5938)
- [P] Hippke's "Interstellar communication" series compares information carriers and concludes that photons are superior by orders of magnitude for fast point-to-point links. It covers solar-gravitational-lens relays and optical laser colours, and a network series asks where relay nodes in our Solar System could be. — [IV, arXiv:1711.07962](https://arxiv.org/abs/1711.07962); [II, arXiv:1706.05570](https://arxiv.org/abs/1706.05570); [X, arXiv:1804.01249](https://arxiv.org/abs/1804.01249); [Network I, arXiv:1912.02616](https://arxiv.org/abs/1912.02616); [Network III, arXiv:2104.09564](https://arxiv.org/abs/2104.09564)
- [P] The Arecibo message (November 1974) was documented in Icarus 26, 462–466. [M] It consists of 1679 bits, the product of the primes 23 and 73, so the only non-trivial raster is 23×73 (or 73×23). It is the canonical example of a self-describing raster in which a semiprime length is the framing clue. — [Staff at NAIC 1975, Icarus 26, 462](https://doi.org/10.1016/0019-1035(75)90116-5)

#### Entropy and block-entropy tests
- [M] For a symbol sequence with alphabet size k, compute the block entropies H_n (over n-grams), the conditional entropies h_n = H_{n+1} − H_n, and the entropy rate h = lim h_n. Comparisons:
  - i.i.d. noise: h_n ≈ log k for all n.
  - Periodic sequence: h_n → 0 at n ≈ the period.
  - Markov sources: h_n plateaus at n = the order.
  - Natural language: h_n decreases slowly over many n.

  The excess entropy E = Σ(h_n − h) measures stored structure.
- [M] Shannon's estimate of English entropy was about 1 bit per character; it is the classic benchmark.
- [P] Grassberger estimates English at about 1.82 bit/character from 135 GB of Project Gutenberg text, extrapolating to about 0.8 bit/character for infinitely long input. Long-range structure therefore keeps lowering h_n out to very long blocks, unlike Markov or periodic processes. — [Grassberger 2002, arXiv:physics/0207023](https://arxiv.org/abs/physics/0207023)

#### Zipf and higher-order entropy as applied to animal communication and proposed for SETI
- [P] McCowan, Hanser & Doyle introduced information-theoretic tools (Zipf's statistic, first- and higher-order entropies) for comparing animal repertoires, applied to bottlenose dolphin whistles. — [McCowan, Hanser & Doyle 1999, Animal Behaviour 57, 409–419](https://doi.org/10.1006/anbe.1998.1000)
- [P] Doyle, McCowan, Johnston & Hanser apply these measures to SETI. They conclude that the study of non-human communication can provide "quantitative general measures of communicative complexity". (Semantic Scholar TL;DR; abstract not available.) — [Doyle et al. 2011, Acta Astronautica 68, 406–417](https://doi.org/10.1016/j.actaastro.2009.11.018)
- [M] Their reported numbers: the log rank–log frequency (Zipf) slope is about −1 for human languages; adult dolphin whistles have a slope of about −0.95, while infant dolphins are flatter (about −0.82) and babbling-like. The slope of entropy against order (Shannon entropic order) is proposed as a measure of rule structure. These values and the claim that the Zipf slope indicates communicative complexity are contested: random-typing ("monkey") texts also produce Zipf-like slopes near −1, so a Zipf slope alone is not diagnostic (Li 1992, [M]).
- [P] Elliott published SETI-oriented work on detecting the signature of intelligent authorship in an unknown signal, including "Detecting the signature of intelligent life" (2010) and "Constructing the matrix" (2012), both in Acta Astronautica. [M] The approach looks for language-like structure (a word-length distribution, function-word-like high-frequency tokens, phrase-length chunking around 6–9 words, and recurring separators) in symbolised data. — [Elliott 2010, Acta Astronautica 67, 1419](https://doi.org/10.1016/j.actaastro.2010.02.018); [Elliott 2012, Acta Astronautica 78, 26](https://doi.org/10.1016/j.actaastro.2011.03.006)

#### Compression and Lempel–Ziv complexity
- [P] Lempel & Ziv 1976 defined the complexity of finite sequences (the LZ76 parse). — [IEEE Trans. Inf. Theory 22, 75–81](https://doi.org/10.1109/tit.1976.1055501)
- [M] Normalised LZ76 complexity c(n)·log_k(n)/n converges to the entropy rate for stationary ergodic sources. It is biased upward for short n, so compare against surrogates of the same length. Practical compressors (gzip, bzip2, xz/LZMA) give a compression ratio usable as an upper bound on the entropy rate, but header overhead and finite windows (32 KB for gzip) distort results for sequences shorter than about 10⁴ symbols or with structure longer than the window.

#### Sample, approximate and permutation entropy
- [P] Approximate entropy (ApEn) "can classify complex systems, given at least 1000 data values". — [Pincus 1991, PNAS 88, 2297](https://doi.org/10.1073/pnas.88.6.2297)
- [P] Sample entropy (SampEn) was introduced because ApEn is ill-suited to short, noisy data. [M] ApEn counts self-matches and is therefore biased, and SampEn removes them. Typical parameters are m = 2 and r = 0.1–0.25 × SD; SampEn is undefined when no template matches occur. — [Richman & Moorman 2000, Am J Physiol 278, H2039](https://doi.org/10.1152/ajpheart.2000.278.6.h2039)
- [P] Permutation entropy (Bandt & Pompe) compares neighbouring values (ordinal patterns of order d). It is fast, robust to observational noise, invariant under nonlinear monotonic transformations, and for chaotic systems behaves like Lyapunov exponents. [M] Use d = 3–7 with N ≫ d! (a rule of thumb is N ≥ 5·d!). Normalise by log d!. Forbidden ordinal patterns indicate determinism. Ties (from quantised data) must be broken randomly or handled explicitly. — [Bandt & Pompe 2002, PRL 88, 174102](https://doi.org/10.1103/physrevlett.88.174102)

#### Statistical complexity and ε-machines
- [P] Crutchfield & Young reconstruct "minimal equations of motion from the recursive structure of measurement sequences" (ε-machines). Statistical complexity C_μ is the entropy of the causal-state distribution, which separates structure from randomness: both i.i.d. noise and periodic sequences have low C_μ. — [Crutchfield & Young 1989, PRL 63, 105](https://doi.org/10.1103/physrevlett.63.105)
- [M] The practical inference algorithm is CSSR (Shalizi & Klinkner 2004); it needs a maximum history length L with roughly k^L ≪ N. The complexity–entropy plane (Rosso et al. 2007) pairs normalised permutation entropy with Jensen–Shannon statistical complexity, which separates chaos from coloured noise.

#### Mutual-information decay
- [P] Lin & Tegmark: mutual information between symbols separated by distance d decays **exponentially** for any probabilistic regular grammar (Markov or hidden Markov process), but can decay as a **power law** for context-free grammars, as in natural languages. — [Lin & Tegmark 2017, Entropy; arXiv:1606.06737](https://arxiv.org/abs/1606.06737)
- [M] In astrophysical time series, power-law-like decay of mutual information also arises from 1/f (red) noise and non-stationarity. Test against surrogates matched in power spectrum (IAAFT, Section 8) before interpreting a power-law decay.

#### Recurrence analysis
- [P] Recurrence plots and recurrence quantification analysis are reviewed comprehensively by Marwan et al. — [Marwan et al. 2007, Phys. Rep. 438, 237](https://doi.org/10.1016/j.physrep.2006.11.001)
- [M] Quantities: recurrence rate RR (fix it at about 1–5% by choosing ε), determinism DET (the fraction of recurrent points on diagonal lines of length ≥ l_min = 2), laminarity LAM and the diagonal-line entropy. Embedding is chosen by the first minimum of mutual information (delay) and false nearest neighbours (dimension). Python packages: `pyunicorn` and `PyRQA` [M].

#### Small-sample bias of entropy and mutual-information estimators
- [P] Paninski proved consistency and central limit theorems for common discretised estimators and showed their bias for undersampled distributions. — [Paninski 2003, Neural Computation 15, 1191](https://doi.org/10.1162/089976603321780272)
- [M] The plug-in (maximum-likelihood) entropy is biased low by about (m − 1)/(2N) nats, with m the number of occupied bins. **Miller–Madow** adds this term back (Miller 1955).
- [P] **Grassberger:** analytic estimators for N points in M boxes with N/M finite, with explicitly known biases that are much smaller than most alternatives. A 2021 generalisation of Schürmann's estimator is "completely free of bias" for special parameters. [M] The 2003 form is H = log N − (1/N) Σ n_i G(n_i), with G(n) = ψ(n) + ½(−1)ⁿ[ψ((n+1)/2) − ψ(n/2)]. — [Grassberger 2003, arXiv:physics/0307138](https://arxiv.org/abs/physics/0307138); [Grassberger 2021, arXiv:2111.11175](https://arxiv.org/abs/2111.11175)
- [P] **NSB (Nemenman–Shafee–Bialek):** near-uniform Dirichlet priors give "disastrous results" for undersampled distributions; an infinite mixture of Dirichlet priors, approximately flat in the entropy, resolves this. A follow-up handles unknown cardinality. — [Nemenman et al. 2002, arXiv:physics/0108025](https://arxiv.org/abs/physics/0108025); [Nemenman 2002, arXiv:physics/0207009](https://arxiv.org/abs/physics/0207009)
- [P] **Chao–Shen:** combines the Horvitz–Thompson adjustment with sample coverage to correct for unseen categories. [M] Coverage C = 1 − f₁/N, where f₁ is the number of singletons; then p̃ = C·p̂ and Ĥ = −Σ p̃ log p̃ / (1 − (1 − p̃)^N). — [Chao & Shen 2003, Environ. Ecol. Stat. 10, 429](https://doi.org/10.1023/a:1026096204727)
- [M] Plug-in mutual information is biased upward by about (m_xy − m_x − m_y + 1)/(2N ln 2) bits (the Panzeri–Treves correction). Shuffle-based bias estimation, which subtracts the mean mutual information of shuffled pairs, is standard. The KSG k-nearest-neighbour estimator (k = 3–10) avoids binning for continuous data.

#### Symbolising a continuous, noisy, irregularly sampled series
- [M] Options:
  - Equal-probability (quantile) binning into k = 2–8 symbols. This is robust to outliers and maximises the alphabet entropy.
  - SAX (piecewise aggregate approximation followed by Gaussian breakpoints).
  - Ordinal patterns (Bandt–Pompe).
  - Threshold crossing or event times (inter-pulse intervals).
  - Derivative sign.

  For irregular sampling, either resample on a uniform grid only where the gap is less than the coherence time (for example with Gaussian-process interpolation, which injects model smoothness), or work with event-based or inter-arrival sequences. A badly placed partition can create spurious entropy and complexity structure (Bollt et al. 2001), so repeat every test across several k and partition schemes and require robustness.

### Inferences
- Minimum viable "encoded structure" battery for Stargazer:
  1. Symbolise in several ways (quantile k = 2, 4, 8; ordinal d = 4–6).
  2. Compute bias-corrected block entropies h_n (with NSB or Grassberger for n beyond log_k N − 2), normalised LZ76, permutation entropy and C_μ.
  3. Compute mutual information I(d) for d = 1 … N/10 with shuffle-bias subtraction.
  4. Run RQA DET and LAM at fixed RR.
  5. Compare every statistic with ≥ 1,000 IAAFT surrogates plus a periodic-plus-red-noise model fitted to the data.

  Report the statistic as a surrogate rank or p-value, never as a raw number.
- (computed) Sample size: block entropy at order n with alphabet k needs N ≫ kⁿ. For k = 4 and N = 10⁴ (a long TESS light curve at 2 min cadence over 2 weeks), reliable orders stop near n ≈ 5–6 (4⁶ = 4096). Most archival light curves (10²–10⁴ points) can therefore test only short-range structure; language-like long-range tests (power-law mutual information, h_n decreasing out to n ≈ 20) need ≳ 10⁶ symbols, achievable only for radio voltage or high-rate photon data.
- Because of the Lachmann–Newman–Moore limit, a negative result from these tests says nothing about encrypted or compressed messages. A positive result most often indicates a deterministic physical process (pulsation, rotation, eclipses) or an instrument (periodic systematics, quantisation). Compare against astrophysical nulls before invoking design.
- Pareidolia guard: mathematical-marker hunts (primes, π, √2, Fibonacci) must fix their target list in advance, count every pattern family tried as a trial, and test against surrogates. KIC 5520878 (Section 5) is the worked cautionary example.

### Gaps
- Abstracts for Doyle et al. 2011 and Elliott 2010 were not available through Crossref or Semantic Scholar, and the papers were not read; the specific Zipf-slope numbers are [M].
- No published application of ε-machines, permutation entropy or RQA to a technosignature candidate was found.
- Messerschmitt's book-length design work was not retrieved; only the arXiv papers were read.

## 7. Unsupervised anomaly detection in time-domain surveys

### Takeaway
The working pattern in astronomy is: feature extraction (hand-crafted light-curve statistics, Gaussian-process or periodogram features, or learned embeddings), then an outlier scorer (isolation forest, density or kNN distance, unsupervised random forest, variational autoencoder), then human review, increasingly with active learning. Contamination by bogus or instrumental data is high: **68%** of SNAD's 277 ZTF DR3 outliers were bogus. Real finds are almost always rare but known astrophysics. Active learning (AAD/PineForest, Astronomaly) roughly doubles the yield of "interesting" anomalies per object inspected, compared with static scoring.

### Cited Findings
- [P] Giles & Walkowicz 2018: density-based clustering (DBSCAN-like) of Kepler long-cadence feature vectors, using Boyajian's Star as ground truth. It flagged fewer than 4% of each quarter as outliers; these included unusual phenomena, objects rare in the field and data artefacts. 2020 follow-up: kNN-distance outlier scores for all 2.8 million long-cadence light curves of about 200k objects (Q1–17), testing PCA, feature removal and the choice of k. Rare classes (fewer than 1000 members) consistently score higher. — [arXiv:1812.07156](https://arxiv.org/abs/1812.07156); [arXiv:2003.00109](https://arxiv.org/abs/2003.00109)
- [P] Wheeler & Kipping (the "Weird Detector"): nonparametric phase-dispersion minimisation for short-duty-cycle periodic signals of arbitrary shape, run on 161,786 Kepler sources. It recovers most known KOIs, Boyajian's Star and a comet. — [arXiv:1903.09213](https://arxiv.org/abs/1903.09213)
- [P] Martínez-Galarza et al. combine tree-based anomaly scores (unsupervised random forest) with manifold learning and clustering, so that anomalies are grouped with their "analogs"; they used Kepler data augmented by Gaia colour and luminosity. Crake & Martínez-Galarza applied the method to TESS: anomalies include rapid pulsators, long-period deep eclipsing binaries and young stellar objects with dimming from obscuration, and the anomalous types differ between Kepler and TESS (anomaly is relative to the parent population). — [arXiv:2009.06760](https://arxiv.org/abs/2009.06760); [arXiv:2301.10264](https://arxiv.org/abs/2301.10264)
- [P] SNAD on the Open Supernova Catalog: Gaussian-process fits, then dimensionality reduction, then isolation forest, then expert review. Of 81 candidates (1.4% of about 2000 objects), **27 (33%)** were astrophysically peculiar; 16 objects labelled as supernovae were likely quasars or stars. — [Pruzhinskaya et al. 2019, arXiv:1905.11516](https://arxiv.org/abs/1905.11516)
- [P] SNAD on ZTF DR3 (3 fields, more than 2.25 million objects; 4 algorithms): **277 outliers → 188 (68%) bogus** (image-subtraction artefacts, star–asteroid overlaps), 66 (24%) known sources and 23 (8%) uncatalogued objects. A simple two-dimensional relation helps filter bogus light curves. Code: `zwad`. — [Malanchev et al. 2021, arXiv:2012.01419](https://arxiv.org/abs/2012.01419)
- [P] Active Anomaly Discovery (AAD) modifies isolation-forest weights after each expert label. On the Open Supernova Catalog it found about 80% more true anomalies than a static isolation forest among the top about 2% of scores. — [Ishida et al. 2021, arXiv:1909.13260](https://arxiv.org/abs/1909.13260)
- [P] `coniferest` (open source) implements isolation forest, AAD and PineForest. Real–bogus scores from a separate classifier can be folded into AAD to suppress plane and satellite tracks, bad CCD columns and ghosts. A supernova-score hybrid followed in 2026. — [Kornilov et al. 2024, arXiv:2410.17142](https://arxiv.org/abs/2410.17142); [Semenikhin et al. 2024, arXiv:2409.10256](https://arxiv.org/abs/2409.10256); [Semenikhin et al. 2026, arXiv:2603.09762](https://arxiv.org/abs/2603.09762)
- [P] **2025:** SNAD ran PineForest on ZTF DR23 fields also observed by LSSTComCam. Of 400 visually inspected candidates, they found 6 new variable stars and refined 6 known ones. An M-dwarf flare catalogue came from more than 35 million ZTF DR8 light curves. — [Pruzhinskaya et al. 2025, arXiv:2507.06217](https://arxiv.org/abs/2507.06217); [Voloshina et al. 2024, arXiv:2404.07812](https://arxiv.org/abs/2404.07812)
- [P] Astronomaly: active learning that gives personalised recommendations and separates "interesting" anomalies from artefacts and from rare-but-uninteresting sources. It roughly **doubles** the interesting anomalies found in the first 100 objects viewed (Galaxy Zoo plus simulations). At scale, it was applied to about 4 million DECaLS galaxies with deep-learned representations, an isolation forest and active learning. "Protege" (2024) adds optimised human labelling. — [Lochner & Bassett 2021, arXiv:2010.11202](https://arxiv.org/abs/2010.11202); [Lochner & Bassett 2022 guide, arXiv:2201.10189](https://arxiv.org/abs/2201.10189); [Etsebeth et al. 2023, arXiv:2309.08660](https://arxiv.org/abs/2309.08660); [Lochner & Rudnick 2024, arXiv:2411.04188](https://arxiv.org/abs/2411.04188)
- [P] **November 2026:** enhanced isolation-forest variants (pooled gain, kurtosis-based, hypervolume) were evaluated on 1,031 MeerKAT extended sources with 52-dimensional BYOL features. — [Thavhana et al. 2026, arXiv:2611.21674](https://arxiv.org/abs/2611.21674)
- [P] Radio analogues: a β-VAE plus random forest for cadences (Ma et al. 2023) and HDBSCAN clustering of hits (Jacobson-Bell et al. 2024); see Section 2. Contardo & Hogg's mid-infrared "contextual anomaly" search predicts expected fluxes with ML and flags residual outliers; see Section 4. — [arXiv:2301.12670](https://arxiv.org/abs/2301.12670); [arXiv:2411.16556](https://arxiv.org/abs/2411.16556); [arXiv:2403.18941](https://arxiv.org/abs/2403.18941)
- [M] Typical light-curve features: amplitude, standard deviation, skewness, kurtosis, MAD, beyond-1σ fraction, Stetson K, von Neumann η, linear trend, Lomb–Scargle peak period, power and false-alarm probability, periodogram percentiles, flux percentile ratios, and Gaussian-process length scales. The `light-curve` package (SNAD) and `feets` implement them. Isolation forest defaults are 100 trees and a subsample of 256 (Liu, Ting & Zhou 2008).

### Inferences
- Plan for most flagged anomalies to be instrumental or bogus: 68% in ground-based difference imaging. Build the real–bogus, known-catalogue cross-match and asteroid or satellite filters **before** the anomaly scorer, not after.
- "Anomalous" depends on the reference population (Crake & Martínez-Galarza), so scores are not portable across surveys. Store the population definition and model version with each score.
- Use active learning with explicit label classes such as "artefact", "known-rare" and "unexplained", and keep a held-out injected-anomaly set (for example synthetic Arnold-type transits or step functions) to measure recall per anomaly type (Section 8).
- No time-domain anomaly search has produced a technosignature candidate that survived vetting. The realistic output of Stargazer's anomaly module is a ranked list of rare astrophysics and systematics, which is still valuable.

### Gaps
- No calibrated instrumental fraction was found for Kepler/TESS anomaly lists (only qualitative "includes data artefacts").
- The Wheeler & Kipping candidate-to-artefact ratio and the 2025–2026 Rubin-era SETI anomaly programmes were not checked.

## 8. Statistics and validation

### Takeaway
A defensible pipeline needs five things:
1. **Calibrated false-alarm probabilities** for each detector, from analytic formulas only where their assumptions hold and otherwise from empirical nulls (OFF scans, inverted or scrambled data, surrogates, time slides).
2. **Explicit trials accounting** (the look-elsewhere effect).
3. **Injection–recovery** to measure completeness.
4. **Pre-registered analysis choices**, because post hoc partitions are how pattern claims fail.
5. **Error control that stays valid as the search grows**: online FDR (LORD/SAFFRON) or alpha-investing for streams of tests, and e-values or e-processes for continuous monitoring and optional stopping.

### Cited Findings

#### False-alarm probability and look-elsewhere
- [P] Lomb–Scargle periodograms: conceptual assumptions and practical pitfalls (frequency grids, aliasing from the window function, false-alarm-probability estimation). — [VanderPlas 2018, ApJS 236, 16; arXiv:1703.09824](https://arxiv.org/abs/1703.09824)
- [M] Baluev (2008) gives an analytic upper bound on the maximum-peak false-alarm probability, FAP ≈ 1 − (1 − e^{−z})·exp(−W·√z·e^{−z}) for the normalised power z with effective bandwidth W. It is reliable for small FAP; bootstrap is preferred when the number of points is small or the noise is red. [M] Implemented in `astropy.timeseries.LombScargle.false_alarm_probability(method='baluev'|'bootstrap'|'davies'|'naive')`.
- [P] Gross & Vitells: the trials factor for a search over a continuous parameter can be estimated from the number of up-crossings of a low threshold in a few background simulations, using Davies' result; asymptotically it grows linearly with the local significance. — [Gross & Vitells 2010, arXiv:1005.1891](https://arxiv.org/abs/1005.1891)

#### Surrogate data
- [P] The method of surrogate data tests a null hypothesis such as "linear Gaussian process" by computing a statistic on surrogates generated under that null. — [Theiler et al. 1992, Physica D 58, 77](https://doi.org/10.1016/0167-2789(92)90102-s)
- [P] Schreiber & Schmitz review surrogates and their caveats, and introduce **iterative amplitude-adjusted Fourier-transform (IAAFT)** surrogates. These handle nonlinear rescalings of Gaussian linear processes by preserving both the autocorrelation (power spectrum) and the amplitude distribution, which simple amplitude adjustment fails to do (it causes spurious nonlinearity detections). Non-stationary nulls need separate treatment. — [Schreiber & Schmitz 2000, arXiv:chao-dyn/9909037](https://arxiv.org/abs/chao-dyn/9909037); [arXiv:chao-dyn/9909041](https://arxiv.org/abs/chao-dyn/9909041); [arXiv:chao-dyn/9904023](https://arxiv.org/abs/chao-dyn/9904023)
- [M] Phase-randomised surrogates keep the Fourier amplitudes and draw uniform random phases. End-point mismatch creates spurious high frequencies, so trim the series so that the first and last values (and slopes) match. A one-sided rank test with M surrogates gives a minimum p of 1/(M+1): use M = 19 for α = 0.05 or M = 999 for α = 0.001. For unevenly sampled data, use Gaussian-process or CARMA model surrogates, or time-shuffled residuals, instead of FFT surrogates.

#### Red-noise null models
- [M] Astrophysical variability is red (power ∝ f^{−α}, α ≈ 1–2, for example quasars and stellar granulation). Testing periodicity against white noise produces false periods. Fit a bending power law or CARMA/damped-random-walk model (for example with `celerite2`) and test against that (Vaughan 2010, 2016). A single "period" in a red-noise series spanning fewer than about 5 cycles is not significant.
- [P] Contested periodicity and structure claims in plate data hinge on the null: Watters et al. show that a uniform-random background is false (edge and corner gradients, strips, clusters). The correct null distribution removes the claimed Earth-shadow deficit, and normalising by observing days removes the nuclear-test correlation. — [Watters et al. 2026, arXiv:2601.21946](https://arxiv.org/abs/2601.21946)

#### Injection–recovery and empirical nulls
- [P] Margot et al.: end-to-end injection recovery was 94.0% for their pipeline versus about 6% for a pipeline using the standard summed product. Efficiency and duty cycle are folded into prevalence limits and the Drake figure of merit. — [arXiv:2308.02712](https://arxiv.org/abs/2308.02712)
- [P] Kepler DR25: the Robovetter was run on real threshold-crossing events and on simulated datasets, giving a catalogue with measured completeness and reliability (8054 KOIs, 4034 candidates). [M] The simulated datasets are injected transits (completeness), and inverted and season-scrambled light curves (reliability). The inverted and scrambled sets provide false alarms with the real systematics but no real transits. — [Thompson et al. 2018, arXiv:1710.06758](https://arxiv.org/abs/1710.06758)
- [P] SETI@home calibrated with artificial "birdies" that model persistent ET signals. — [Anderson et al. 2025, arXiv:2506.14737](https://arxiv.org/abs/2506.14737)
- [P] Tellis & Marcy set laser thresholds by injecting synthetic lines and analysing them blind. — [arXiv:1704.02535](https://arxiv.org/abs/1704.02535)

#### Bayesian model comparison
- [P] Kipping & Gray model the Wow! signal as a stochastic repeater with a likelihood emulator built on the Big Ear logs. The maximum a posteriori solution has a likelihood of 32.3%; the 2σ credible interval for signal duration is 72 s < T < 77 min and for repeat rate 0.043 < λ < 59.8 per day. Adding 192 h of META, Hobart and ATA follow-up drops the peak likelihood to 1.78% (2.4σ tension), and about 62 more days of observation would be needed to reach 3σ. This is a worked example of quantifying "has it failed to repeat?". — [Kipping & Gray 2022, arXiv:2206.08374](https://arxiv.org/abs/2206.08374)
- [M] For candidate-against-null comparison, compute Bayes factors with nested sampling (`dynesty`, `UltraNest`). Read them on the Jeffreys or Kass–Raftery scales (ln B > 5 is "strong"), and report prior sensitivity. Bayes factors do not remove the look-elsewhere effect when the prior over the search space is broad; the prior volume plays the role of the trials factor.

#### Blind analysis and pre-registration
- [M] Blind analysis (Klein & Roodman 2005, Annu. Rev. Nucl. Part. Sci.) hides the signal region or adds a hidden offset until cuts are frozen. LIGO uses time slides and blind injections; the "Big Dog" hardware injection (2010) tested the collaboration's response end to end.
- [P] Example of why pre-registration alone is insufficient: a single-author preprint describes a "pre-registered empirical test" on VASCO transients (p = 0.00066) and combines tests to 3.6–4.7σ by Fisher's method while acknowledging dependence on the independence assumption. Pre-registration does not fix an invalid null or dependent tests. — [Cann 2026, arXiv:2604.06234](https://arxiv.org/abs/2604.06234)

#### Error control in an open-ended, growing search
- [P] Benjamini–Hochberg (1995) controls the false discovery rate (FDR), the expected proportion of false rejections; it is equal to the family-wise error rate when all nulls are true and smaller otherwise. — [JRSS-B 57, 289](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [P] **Alpha-investing** (Foster & Stine) is an adaptive sequential family of procedures that controls mFDR at every rejection: each test spends α-wealth, and each discovery earns wealth back. — [Foster & Stine 2008, JRSS-B 70, 429](https://doi.org/10.1111/j.1467-9868.2007.00643.x)
- [P] **Online FDR** (Javanmard & Montanari) controls FDR or mFDR for an ordered, possibly infinite stream of hypotheses, where each decision is made before later p-values are seen. [M] The LORD and LOND rules come from this work, and SAFFRON and ADDIS are later improvements; Python: `onlineFDR` is in R, and Python ports exist. — [arXiv:1502.06197](https://arxiv.org/abs/1502.06197); [arXiv:1603.09000](https://arxiv.org/abs/1603.09000)
- [P] **E-values** are random variables with expectation at most 1 under the null (betting scores, Bayes factors, likelihood ratios, stopped supermartingales); they combine by averaging. **e-BH** controls FDR under arbitrary dependence between e-values, unlike p-value BH. — [Vovk & Wang 2021, arXiv:1912.06116](https://arxiv.org/abs/1912.06116); [Wang & Ramdas 2022, arXiv:2009.02824](https://arxiv.org/abs/2009.02824); [Ramdas & Wang 2024 book, arXiv:2410.23614](https://arxiv.org/abs/2410.23614)
- [P] **Safe anytime-valid inference:** e-processes and confidence sequences remain valid at all stopping times, so data can be monitored continuously and analyses stopped or continued "for any reason". Built on test martingales, which are nonnegative martingales starting at 1. — [Ramdas et al. 2023, arXiv:2210.01948](https://arxiv.org/abs/2210.01948)
- [P] The first technosignature search I found applying FDR control is OVRO-LWA (2026). — [arXiv:2606.04304](https://arxiv.org/abs/2606.04304)

### Inferences
- Concrete Stargazer policy:
  1. Every detector emits a p-value calibrated on an **empirical** null (OFF scans or relabelled cadences for radio; inverted or scrambled or IAAFT surrogates for light curves; plate- and schedule-preserving permutations for images), plus its trials count.
  2. Within a run, apply BH (or e-BH when tests are dependent) at q = 0.01–0.05 to produce "flags for review".
  3. Across the programme's lifetime, run an online FDR rule (LORD++ or SAFFRON) or alpha-investing over the stream of runs, so that adding archives does not silently inflate false discoveries.
  4. For targets monitored repeatedly, use an e-process (a running likelihood ratio against the null), which may be checked after every new epoch without penalty.
- (computed) Scale: 10⁹ light curves × 10³ effective frequencies is 10¹² trials. A family-wise 3σ global significance (p = 1.35×10⁻³) then needs a local p of about 10⁻¹⁵ (about 8σ). Tool output should always state the global and local significance together.
- Freeze, version and publish the null-model choice, symbolisation, statistics list and thresholds before running on new data. Changing them after looking at a candidate should create a new pre-registered test on independent data.

### Gaps
- No technosignature paper was found using alpha-investing or e-values; the transfer is my inference.
- The Baluev formula, the Klein & Roodman reference and the Kepler inverted/scrambled details are [M].

## 9. Contaminants by domain

### Takeaway
In every domain, the dominant false positives are human-made or instrumental:
- **Radio:** terrestrial and satellite emitters (including unintended broadband emission from Starlink, now detectable in about 30% of some low-frequency images), intermodulation products locked to clock harmonics (BLC1), and local transient sources (perytons).
- **Optical:** satellite glints, cosmic rays, asteroids, emulsion and copy defects, and difference-imaging artefacts.
- **Kepler/TESS:** rolling-band, momentum-dump and scattered-light systematics, plus contamination by nearby variables (ephemeris matching).
- **Everywhere:** timing standards (minute-level JD errors), sampling aliases and quantisation create spurious periodicities.

### Cited Findings

#### Radio-frequency interference and satellite constellations
- [P] LOFAR detected 68 Starlink satellites at 110–188 MHz, with broadband features across the band and narrowband (< 12.2 kHz) lines at 125, 135, 143.05, 150 and 175 MHz. The 143.05 MHz line is reflected GRAVES radar; the rest is intrinsic and seen in 47 of 68 satellites. Flux densities are 0.1–10 Jy broadband and 10–500 Jy narrowband. — [Di Vruno et al. 2023, arXiv:2307.02316](https://arxiv.org/abs/2307.02316)
- [P] **2025:** EDA2 at the SKA-Low site analysed about 76 million all-sky images over about 29 days: 112,534 detections of 1,806 Starlink satellites. Up to about 30% of images in some datasets contain a satellite, emission falls in ITU-protected bands (13 satellites at 73.00–74.60 MHz; 703 at 150.05–153.00 MHz), and FM reflections appear at 99.7 MHz. — [Grigg et al. 2025, arXiv:2506.02831](https://arxiv.org/abs/2506.02831); [Grigg et al. 2023, arXiv:2309.15672](https://arxiv.org/abs/2309.15672)
- [P] **2026:** SNIFFLES I (Mopra and ATCA, 1–26 GHz) made 4,629 tracked observations totalling 375.9 h of Starlink, OneWeb, Amazon Leo and Guowang satellites, characterising intended, unwanted and unintended emission. — [Indermuehle & Lourenco 2026, arXiv:2608.12999](https://arxiv.org/abs/2608.12999)
- [P] Proxima search: the negative-drift excess among hits "is due to non-geosynchronous satellites". — [Smith et al. 2021, arXiv:2111.08007](https://arxiv.org/abs/2111.08007)
- [P] BLC1: an intermodulation product of local time-varying interferers, with look-alikes at clock-oscillator harmonics. Sheikh et al. also overlaid satellites above the horizon using space-track.org. — [Sheikh et al. 2021, arXiv:2111.06350](https://arxiv.org/abs/2111.06350)
- [P] Perytons were swept-frequency terrestrial transients mimicking dispersion. — [Petroff et al. 2015, arXiv:1504.02165](https://arxiv.org/abs/1504.02165)

#### Optical: glints, cosmic rays, artefacts
- [P] Geosynchronous glints: 0.1–0.3 s flashes at about 30–40 per day per deg² brighter than about 11 mag; some are point-like in long exposures and confound searches for fast transients. A GN-z11 "flash" is consistent with a glint. — [Nir et al. 2021, arXiv:2011.03497](https://arxiv.org/abs/2011.03497)
- [P] ZTF anomaly lists are dominated by image-subtraction artefacts and asteroid overlaps (68%). Plane and satellite tracks, bad columns and ghosts contaminate anomaly outputs. — [Malanchev et al. 2021, arXiv:2012.01419](https://arxiv.org/abs/2012.01419); [Semenikhin et al. 2024, arXiv:2409.10256](https://arxiv.org/abs/2409.10256)
- [P] Photographic plates: emulsion and copy-plate defects, and digitisation artefacts; checking against a second digitisation (SuperCOSMOS) removed 3,592 artefacts in one VASCO search. — [Hambly & Blair 2024, arXiv:2402.00497](https://arxiv.org/abs/2402.00497); [Solano et al. 2022, arXiv:2206.00907](https://arxiv.org/abs/2206.00907)
- [M] Cosmic rays in CCD images are removed by L.A.Cosmic (Laplacian edge detection; van Dokkum 2001) or `astroscrappy`, and by multi-exposure median combination. In spectra they mimic unresolved emission lines, hence the persistence requirement in laser searches.

#### Kepler and TESS systematics
- [P] Ephemeris matching identified **685 KOIs as false positives** (Q1–Q12) caused by contamination from eclipsing binaries, variables and other transiting planets, through direct PSF overlap, CCD crosstalk, column anomalies and reflections. — [Coughlin et al. 2014, arXiv:1401.1240](https://arxiv.org/abs/1401.1240)
- [M] Other Kepler and TESS systematics:
  - Kepler: rolling-band (Moiré) pattern noise, momentum dumps (about every 3 days in K2, with thruster firings every about 6 h), sudden pixel-sensitivity dropouts after cosmic-ray hits, quarter-to-quarter discontinuities and 90° rolls.
  - TESS: scattered Earth and Moon light near orbit perigee, momentum dumps (every few days), and a 13.7-day orbital gap.

  Cotrending basis vectors (PDC-MAP/PDC-msMAP) and the TESS quality flags must be applied, and every anomaly should be tested for coincidence with momentum-dump and quality-flag times across many stars on the same CCD.
- [P] Boyajian's Star dips were confirmed as not instrumental by checking Kepler data processing; the kind of check required is described in Boyajian et al. 2016. — [arXiv:1509.03622](https://arxiv.org/abs/1509.03622)

#### Aliasing, quantisation and timing
- [P] Lomb–Scargle aliasing: the observed peaks are the true frequency convolved with the window function; daily, lunar and annual sampling create aliases. — [VanderPlas 2018, arXiv:1703.09824](https://arxiv.org/abs/1703.09824)
- [P] Time standards adopted by different groups (HJD versus BJD, UTC versus TDB) have differed by as much as a minute. — [Eastman, Siverd & Gaudi 2010, arXiv:1005.4415](https://arxiv.org/abs/1005.4415)
- [P] Centimetre-per-second barycentric corrections require more than naive subtraction of velocities; photon-weighted midpoint errors are about 10 cm/s. This matters for laser-line wavelength stability and drift-rate prediction. — [Wright & Eastman 2014, arXiv:1409.4774](https://arxiv.org/abs/1409.4774); [Tronsgaard et al. 2019, arXiv:1908.00991](https://arxiv.org/abs/1908.00991)
- [P†] turboSETI blanks the DC bin of each coarse channel (a polyphase filterbank and quantisation artefact). — [turbo-seti source](https://pypi.org/project/turbo-seti/2.3.2/)
- [M] Coarse digitisation (2- or 8-bit samplers) and FFT scalloping generate spurious harmonics and spectral-leakage spurs at fixed fractions of the coarse-channel width. Integer-step photometry (low-count or saturated pixels) produces discrete flux levels that inflate symbolic-structure statistics. Detector readout periodicities (for example TESS 2 min and 20 s cadences, and Kepler's 29.4 min long-cadence sampling with a 1/LC Nyquist frequency of about 283 µHz) show up as artificial "periods".

### Inferences
- Stargazer needs a **contaminant database** checked before scoring:
  - Radio: ITU allocations, observatory RFI logs, and clock harmonics (10 MHz, 125 MHz and their multiples) and intermodulation sums and differences.
  - Satellite ephemerides: TLEs for the observation epoch, with the full constellation including unintended emission bands.
  - Optical: minor-planet positions and the times of known spacecraft systematics (momentum dumps, quality flags).
  - Time-domain: survey window functions and their aliases.
- Population-level co-occurrence (the same anomaly at the same time in many unrelated targets on the same detector, or at the same frequency in many pointings) is the most powerful generic systematics test. It is what found the BLC1 look-alikes and the Kepler ephemeris matches.
- Convert every timestamp to BJD_TDB before any period or synchronisation test; record the original standard.

### Gaps
- Satellite-constellation interference specifically at 1–2 GHz in Breakthrough Listen data (as opposed to low frequencies) was found only in the 2026 SNIFFLES abstract, without numbers.
- The Kepler and TESS systematic timescales listed above are [M].

## 10. Verification and reporting

### Takeaway
Verification has been formalised in four instruments:
- **Sheikh et al. 2021** give a ten-step checklist for narrowband candidates, built on BLC1.
- The **Rio Scale 2.0** communicates significance: J = A + B + C − 20, δ = 10^{(J−10)/2}, R = Q·δ. BLC1 scored R = 0.
- The **IAA Declaration of Principles** was revised (a draft in October 2025, an approval vote in December 2025–January 2026, and reported ratification in 2026). It requires independent confirmation by multiple facilities and organisations with different instruments before any claim, prompt retraction disclosure, open data in at least two geographically separate repositories, and peer-reviewed verification reports.
- The **St Andrews SETI Post-Detection Hub** (since 2022) coordinates post-detection research.

The 2024–2025 re-analyses of the Wow! signal (Arecibo Wow! I and II) refine its frequency to 1420.726 ± 0.005 MHz and its peak flux to more than 250 Jy. They favour an astrophysical H I-cloud flare over interference or a technosignature. This is a hypothesis, not a confirmation.

### Cited Findings

#### Sheikh et al. 2021 verification framework (verbatim steps, read in full text)
- [P] Applied "once known astrophysical origins have been ruled out", for narrowband signals of interest:
  1. Verify that all instrumentation was functioning correctly.
  2. Verify that the signal of interest was not present in the off-source observations at a lower S/N threshold.
  3. Check for catalogued RFI at the same frequency at the observatory.
  4. Compare the drift-rate evolution with known accelerational and electronic drifts from human technology.
  5. Compare the drift-rate evolution with the drift rates and periods expected in the target system and the Solar System.
  6. Search for other instances of the signal of interest in archival data from the same observatory.
  7. Search for similar signals at other frequencies within the same observation, and determine whether they look like RFI.
  8. Assess, qualitatively or quantitatively, whether the similar signals from steps 6 and 7 come from the same phenomenon.
  9. If they do, determine whether they are RFI using off-source observations, or steps 1–8 if necessary.
  10. Re-observe with the same instrument and with other instruments.

  The checklist is written for persistent narrowband signals from single-beam dishes. For multibeam instruments, other beams serve as reference pixels. Re-observation (step 10) can come earlier for periodic or transient signals (for example ones synchronised with a planet's transit). Further recommendations: simultaneous two-site observation, cold-load tests of receiver aliasing, noise injection on the look-alike population, and recording data during slews to study sidelobes. — [Sheikh et al. 2021, arXiv:2111.06350](https://arxiv.org/abs/2111.06350)

#### Rio Scale and Rio 2.0
- [P] The Rio Scale was modelled on the Torino scale for asteroid impacts, as a way to discuss and interpret any claimed discovery of ETI. — [Almár & Tarter 2011, Acta Astronautica 68, 358](https://doi.org/10.1016/j.actaastro.2009.07.007)
- [P] Rio 2.0 revises the scale to build cross-disciplinary consensus, provide a pedagogical tool for public understanding of the verification process, calibrate media expectations, and recommend consistent terminology. A critique ("Deconstructing the Rio Scale: problems of subjectivity and generalization") and a rebuttal were published in 2019. — [Forgan et al. 2019, IJA 18, 336](https://doi.org/10.1017/s1473550418000162); [rebuttal, IJA 18, 492](https://doi.org/10.1017/s1473550418000435)
- [P] Applied to BLC1 with the online calculator: A = 6 (amenable to study or repeat), B = 8 (sure it is not instrumental?), C = 1 (how likely a technosignature?). This gives J = A + B + C − 20 = 0 and δ = 10^{(J−10)/2} = 10⁻⁵; Q = 7, so R = Q·δ ≈ 0. Sheikh et al. note that the scale "does not provide a framework with which to analyse and determine if a signal is interference", lacks granularity between 0 and 1, and propose a logarithmic successor. — [Sheikh et al. 2021 supplementary, arXiv:2111.06350](https://arxiv.org/abs/2111.06350)

#### IAA Declaration of Principles (2010 → 2025 draft → 2026)
- [P] History: first formulated in 1989; reply protocols drafted in 1995; presented to UN COPUOS in 2000; streamlined in 2010. A task group was set up in 2022, a draft was presented at IAC 2024 in Milan, and a revised draft was published in October 2025. — [Garrett et al. 2025, arXiv:2510.14506](https://arxiv.org/abs/2510.14506)
- [P] Key revised principles in the October 2025 draft:
  - **(1) Handling candidate evidence:** the discoverer should authenticate the detection, ideally including "independent observations or other examinations by multiple facilities and by more than one organisation utilising different instrumentation and methods". Candidate information should be handled "with extreme care".
  - **(2) Communicating:** there is no obligation to disclose verification until a discovery is confirmed, but communication may be needed to dispel rumours. If a reported candidate proves not to be extraterrestrial, "this should be promptly disclosed". Speculative statements must be labelled as such. Researchers may decline media engagement, and institutions should protect them.
  - **(3) Communicating verification:** confirmed results go to the public, the scientific community and the UN Secretary-General. The verification report should be peer-reviewed and include the underlying data, analysis process and any information content, and be sent to the IAA, IAU, COSPAR, IISL, ITU, COPUOS and UNOOSA.
  - **(4) Archiving:** data, methods **and code** should be preserved in at least two repositories in different geographic locations, in open formats.
  - **(5) Data and frequency protection:** tamper-proof records, and ITU extraordinary procedures to protect the frequencies.
  - **(6)** Post-detection protocol; **(7)** no reply without international consultation; **(8)** ethical and legal considerations.
  - Scope: intelligent life only, not life in general and not UAP.

  — [Garrett et al. 2025, arXiv:2510.14506](https://arxiv.org/abs/2510.14506)
- [S] **2026:** search summaries report that the IAA SETI Committee approved the new Declaration in a vote in December 2025–January 2026, that it has been formally ratified by the IAA, and that a 5 June 2026 press release described it as the first major update in more than 15 years. The update reportedly creates a permanent Post-Detection Sub-Committee (social science, law, ethics) and addresses protection of researchers against harassment and doxxing, with presentation planned at the IAC in Türkiye in 2026. The ratified text was not read (iaaspace.org is unreachable). — [University of Manchester news (search summary)](https://www.manchester.ac.uk/about/news/beyond-disclosure-day-the-real-world-protocols/); [SETI Institute news, 5 June 2026 (search summary)](https://seti.org/news/beyond-disclosure-day/); [EPSC 2026 abstract (search summary)](https://meetingorganizer.copernicus.org/EPSC2026/EPSC2026-132.html)

#### St Andrews SETI Post-Detection Hub
- [S] Launched in 2022 at the University of St Andrews and coordinated by John Elliott, it is an initiative of the UK SETI Research Network, jointly hosted with the Centre for Global Law and Governance, Durham Law School and Jodrell Bank. It coordinates impact assessments, protocols and treaty work, and notes that no UN procedure exists for a radio detection, unlike for asteroid impacts. — [St Andrews news (search summary)](https://news.st-andrews.ac.uk/?p=22694); [Hub website (search summary)](https://seti.wp.st-andrews.ac.uk/)
- [P] A 2025 white paper (Genevieve, Kovacevic, Elliott, Dominik et al.) recommends that NASA fund post-detection readiness along six dimensions, guided by the Hub's framework. — [arXiv:2507.11587](https://arxiv.org/abs/2507.11587)

#### Wow! signal re-analyses (2022–2025)
- [P] Kipping & Gray 2022: the Wow! signal cannot be excluded as a stochastic repeater (Section 8). — [arXiv:2206.08374](https://arxiv.org/abs/2206.08374)
- [P] Arecibo Wow! I (2024): Arecibo drift scans at 1420 MHz (February–May 2020) found 10 kHz narrowband H I signals resembling the Wow! signal but about 100 times weaker, at multiple locations, identifiable as small cold H I clouds. The authors hypothesise that a transient radiation source (a magnetar flare or soft gamma repeater) suddenly brightened such a cloud through maser-like or superradiant emission. — [Méndez et al. 2024, arXiv:2408.08513](https://arxiv.org/abs/2408.08513)
- [P] Arecibo Wow! II (2025): from previously unpublished Ohio SETI archival data, they place the source in one of two adjacent fields at RA 19h25m02s ± 3s or 19h27m55s ± 3s, Dec −26°57′ ± 20′ (J2000). Peak flux density exceeds 250 Jy, and the frequency is 1420.726 ± 0.005 MHz, implying a higher radial velocity than assumed. Interference is disfavoured and an astrophysical origin supported. — [Méndez et al. 2025, arXiv:2508.10657](https://arxiv.org/abs/2508.10657)
- [P] Earlier claims of a cometary hydrogen source (comets 266P/Christensen and P/2008 Y2) were made and are not supported by the later analyses. — [Paris & Davies 2017, arXiv:1706.04642](https://arxiv.org/abs/1706.04642)

### Inferences
- **What Stargazer should do before labelling anything a "candidate".** Use a tiered vocabulary with clear gates: "hit" (crossed a threshold) → "flag" (passed automated contaminant filters and the within-run FDR) → "signal of interest" (passed the domain checklist below) → "candidate" (re-detected in independent data). Never use "technosignature" for unconfirmed results. This matches IAA Principle 2 ("speculative or unconfirmed statements… clearly identified").
- **Domain checklist** (generalising Sheikh's ten steps to archive data):
  1. Instrument and data-quality flags are clean at the event time.
  2. The event is absent from reference data (OFF scans, other beams, comparison stars on the same detector, a second digitisation) at a **lower** threshold.
  3. Cross-match against contaminant catalogues (RFI, satellites, asteroids, known variables, background galaxies, cosmic rays).
  4. Kinematic or temporal consistency with the target (drift against expected accelerations; a period consistent with no spacecraft or detector cycle).
  5. Search the whole archive for look-alikes at other times, frequencies and targets; if look-alikes exist, classify the population.
  6. Calibrated global significance (empirical null, trials factor, FDR stream position).
  7. Injection–recovery shows the pipeline would have detected and correctly characterised a similar injected signal.
  8. Independent re-detection in a different archive or instrument, ideally by another team.
  9. Rio 2.0 scoring, with stated A, B, C and Q.
- **Reporting:** auto-generate a provenance bundle (raw data identifiers, code version and hash, the parameters and nulls pre-registered, the full funnel counts, injection results) archived in two places (for example Zenodo plus a mirror). Include a default "retraction" pathway: if a flagged item is later explained, update its public record promptly (IAA Principle 2).
- Avoid public framing that encourages pareidolia. In particular, in the contested VASCO-type datasets, display the competing explanations alongside the flag.

### Gaps
- The ratified 2026 IAA text could not be read (iaaspace.org is not reachable); the post-October-2025 changes are known only from search summaries ([S]).
- The St Andrews Hub's own protocol documents and the 2023 Elliott et al. Hub paper were not read.
- The original "Deconstructing the Rio Scale" critique was not read; only the rebuttal's existence was confirmed.
- I did not check whether the Rio 2.0 online calculator is still live or whether a logarithmic successor scale has been published (no arXiv hit was found).
