# Recovering structure from signals and bitstreams of unknown origin and format (state of the art, 21 November 2026)

Scope: the structural layer only. The question is how to get from a raw measured signal, or a pile of bits already suspected to carry information, to its symbols, framing, layout, coding and grammar-like structure. Detection statistics, archives and semantic decipherment are covered by other researchers.

Evidence tags:
- **[P]**: read today in the primary source (arXiv/Zenodo text or abstract, PyPI JSON, repository file on raw.githubusercontent.com, Crossref record).
- **[S]**: taken from a web-search summary only.
- **[R]**: recalled from background knowledge and not re-checked today. Treat [R] items as leads to verify before they go into a specification.

Maturity labels:
- **ESTABLISHED**: textbook or widely deployed.
- **RESEARCH**: published and promising, but unproven on truly unknown data.
- **SPECULATIVE**: an idea, or a self-published claim.

Unreachable sources:
- gnuradio.org, asignin.space, seti.org, boingboing.net and destevez.net could not be fetched (proxy or DNS refusal).
- The GitHub API was not available, so commit dates are unknown. Maintenance status below comes from PyPI release dates, changelogs and READMEs.
- Semantic Scholar and OpenAlex returned HTTP 429 (rate-limited) for most queries.

---

## 1. Blind signal analysis (communications engineering and SIGINT): carrier, bandwidth, modulation, symbol rate, timing, constellation, equalisation

### Takeaway
The physical-layer chain is mature engineering. The steps are: PSD-based carrier and bandwidth estimation, then cyclostationary or power-law spectral lines for symbol rate and carrier, then timing recovery, then blind equalisation (CMA), then constellation clustering. Good open tools exist (GNU Radio, liquid-dsp, inspectrum, SigDigger, URH, SigMF). Automatic modulation classification (AMC), whether cumulant-based or deep-learning, only picks from a closed list of human modulation types. For Stargazer it should be a hypothesis generator, not an oracle. The parameter-estimation steps (rate, carrier, constellation) transfer to unknown signals better than the classifiers do.

### Cited Findings

**Method reference (algorithms to implement)**

- **Carrier and bandwidth estimation** (ESTABLISHED) [R]
  - Welch PSD, then a noise floor from a median or low percentile, then the occupied bandwidth at a −3 dB, −10 dB or 99%-power threshold.
  - The carrier offset of an M-PSK signal shows up as a spectral line at M·f_c in the spectrum of x^M (M-th power method; BPSK uses squaring).
  - Failure modes: low SNR, offset QAM, and spread-spectrum signals that hide below the noise floor.
- **Symbol-rate estimation** (ESTABLISHED) [R]
  - Linearly modulated signals are cyclostationary with cycle frequency α = 1/T_s.
  - Methods include |x|² (envelope) spectral lines, the delay-and-multiply line at R_s, and the spectral correlation function (SCF) or cyclic autocorrelation. FFT Accumulation Method (FAM) and Strip Spectral Correlation Analyser (SSCA) are the standard fast SCF estimators.
  - The line strength scales with excess bandwidth: roll-off β→0 makes the line vanish.
  - Elgenedy & Elezabi propose a family of autocorrelation zero-crossing symbol-rate estimators that work over a continuous range of rates with no candidate set. They are insensitive to roll-off, are extended to multipath channels, and were implemented in a complete DVB-S2 receiver. [P] [arXiv:2001.01692 (2020)](https://arxiv.org/abs/2001.01692)
- **Timing recovery** (ESTABLISHED) [R]
  - Gardner timing-error detector: 2 samples per symbol, works without the carrier.
  - Mueller–Müller detector: 1 sample per symbol, needs decisions.
  - Polyphase-filterbank clock sync (GNU Radio `symbol_sync`, liquid-dsp `symsync`). Oerder–Meyr square-law feed-forward estimator.
- **Blind equalisation and carrier tracking** (ESTABLISHED)
  - Godard's 1980 "self-recovering" equaliser is the origin of the constant-modulus algorithm (CMA). [P] [Godard 1980, IEEE Trans. Commun., doi:10.1109/TCOM.1980.1094608](https://doi.org/10.1109/TCOM.1980.1094608)
  - CMA tap update [R]: w ← w − μ·y(|y|²−R₂)·x*, with R₂ = E|s|⁴/E|s|². Follow-ons are the multi-modulus algorithm and decision-directed LMS, then a Costas or 4th-power PLL for the residual phase.
  - A semidefinite-relaxation CMA exists as research. [P] [arXiv:1808.07232](https://arxiv.org/abs/1808.07232)
- **Constellation recovery** (ESTABLISHED/RESEARCH) [R]
  - After timing and CMA: k-means or GMM on the I/Q symbol cloud, with the number of clusters chosen by BIC/MDL or the silhouette score, and the rotation ambiguity resolved by the 4th-power phase.
  - Gray or other labelling is then unknown, which leaves an M!-fold labelling ambiguity to resolve in the bit domain (see section 2).
- **Cumulant-based AMC** (ESTABLISHED)
  - Swami & Sadler's hierarchical classifier uses fourth-order cumulants (C40, C42) of the symbol-rate samples. [P] [Swami & Sadler 2000, IEEE Trans. Commun., doi:10.1109/26.837045](https://doi.org/10.1109/26.837045) (≈1,128 citations per Semantic Scholar)
  - Normalised C42 theoretical values [R]: BPSK −2, QPSK −1, 16-QAM −0.68, 64-QAM −0.62.
  - The method needs symbol-synchronous samples and some hundreds to thousands of symbols. Failure modes: frequency offset, multipath, and non-Gaussian noise.
- **Cyclostationary-feature AMC and detection** (ESTABLISHED)
  - Sutton et al. (2008) define cyclostationary "signatures" deliberately embedded in a signal and detected through cyclostationary analysis, for detection and network identification. [P] [doi:10.1109/jsac.2008.080103](https://doi.org/10.1109/jsac.2008.080103)
  - Câmara et al. (2019) compare cyclostationary AMC architectures under impulsive (alpha-stable) noise. [P] [doi:10.1109/access.2019.2943300](https://doi.org/10.1109/access.2019.2943300)
  - Inference: an *intentionally* embedded cyclic feature is a plausible design choice for a beacon, which makes the SCF a natural Stargazer probe.
- **Deep-learning AMC** (ESTABLISHED in benchmarks; RESEARCH for out-of-distribution signals)
  - O'Shea, Corgan & Clancy (2016) showed CNNs on raw I/Q outperform expert features, especially at low SNR. [P] [arXiv:1602.04105](https://arxiv.org/abs/1602.04105)
  - O'Shea, Roy & Clancy (2018, RadioML 2018) compare deep nets with a baseline of higher-order moments plus boosted trees. They cover carrier offset, symbol rate and multipath, include over-the-air tests, and list remaining problems. [P] [arXiv:1712.04578](https://arxiv.org/abs/1712.04578)
- **RadioML datasets**
  - RadioML 2018.01A: 24 digital and analog modulation classes, 1024-sample examples; 2,555,904 examples per papers. DeepSig's file name suggests about 2 million, so the counts disagree between sources. [S] [DeepSig datasets page](https://deepsig.squarespace.com/datasets-old); [arXiv:2508.20193](https://arxiv.org/pdf/2508.20193)
  - Licence: CC BY-NC-SA 4.0 (non-commercial), per DeepSig's "old" datasets page. Confirm the current terms. [S] [DeepSig](https://deepsig.squarespace.com/datasets-old)
  - RadioML 2016.10A is still used as a benchmark in 2026 (LLM-AUG augmentation work). [P] [arXiv:2604.17770](https://arxiv.org/abs/2604.17770)
  - Joint detection and AMC on a multi-signal coexisting synthetic set (CRML23). [P] [arXiv:2405.00736](https://arxiv.org/abs/2405.00736)
- **2025–2026 developments**
  - A July 2026 survey of AI-based modulation recognition for communications and radar. [P] [Huang et al., arXiv:2607.23014](https://arxiv.org/abs/2607.23014)
  - "Radio-FM" (August 2026): a family of RF foundation models with dual-channel I/Q processing and masked-reconstruction pretraining on 15 datasets, evaluated on 15 downstream tasks. [P] [arXiv:2608.05793](https://arxiv.org/abs/2608.05793)
  - Inference: self-supervised RF embeddings are the 2026 trend. They could be reused as generic "tokenisers" for RF, but they were pretrained on human modulations.

**Tools**

| Tool | Repo / licence | Status 2026 | Would do for Stargazer | Would not do |
|---|---|---|---|---|
| GNU Radio | github.com/gnuradio/gnuradio; GPL-3.0 (COPYING read [P]). GR4 core is MIT; blocks ported from GR3 keep GPLv3 [S] | 3.10.x is the stable line. GR4 reached RC1 on 22 Mar 2026 and RC3 by 13 Oct 2026; no 4.0.0 final found. Joined the Linux Foundation OCUDU foundation on 22 Sep 2026 [S] ([gnuradio.org news](https://www.gnuradio.org/news/2026-03-22-gr4-release-candidate-1/), [GR4 stewardship](https://www.gnuradio.org/news/2026-05-21-gr4-community-stewardship/)) | Full DSP chain: filters, PLLs, symbol sync, CMA/LMS equalisers, demodulators, file I/O | No automatic "what is this signal" analysis |
| liquid-dsp | github.com/jgaeddert/liquid-dsp; MIT (LICENSE read [P]) | Active: changelog 1.8.0 (2026-06-07) to 1.8.3 (2026-09-26); PyPI `liquid-dsp` 1.8.3, 27 Sep 2026 [P] ([PyPI](https://pypi.org/project/liquid-dsp/)) | Embeddable C library: symsync, equalisers (CMA/LMS/RLS), modems, FEC, framing; Python bindings | No GUI and no blind classification |
| inspectrum | github.com/miek/inspectrum; GPL-3.0 [P] | Packaged in most distributions, Homebrew and conda-forge [P] ([README](https://raw.githubusercontent.com/miek/inspectrum/main/README.md)) | Interactive spectrogram for 100 GB+ files; cursors to measure period and symbol rate and extract symbols; reads SigMF and cf32/cs16 and similar [P] | Manual only, no automation |
| SigDigger (+ suscan, sigutils) | github.com/BatchDrake/SigDigger; GPL-3.0 [P] | README mentions Qt6 builds and an AppImage "latest" development build; work now lands on `master` [P] ([README](https://raw.githubusercontent.com/BatchDrake/SigDigger/master/README.md)) | "Extract information of unknown radio signals": adjustable FSK/PSK/ASK demodulation, bursty-signal analysis, analog video [P] | Linux/macOS only; GUI-centred |
| Universal Radio Hacker (URH) | github.com/jopohl/urh; GPL-3.0 [P] | PyPI 2.10.0 released 17 Dec 2025 [P] ([PyPI](https://pypi.org/project/urh/)) | Demodulation with automatic detection of modulation parameters; custom decodings (for example CC1101 whitening); rule-based automatic protocol-field inference (WOOT'19); fuzzing [P] ([README](https://raw.githubusercontent.com/jopohl/urh/master/README.md)) | Assumes short-range ASK/FSK/PSK bursts, and its field inference assumes many messages |
| SigMF | github.com/sigmf/SigMF; sigmf-python is LGPL-3.0+ [P] | sigmf 1.14.0, 19 Oct 2026 [P] ([PyPI](https://pypi.org/project/sigmf/)) | The interchange format to adopt: JSON `global`/`captures`/`annotations` with `core:datatype`, `core:sample_rate`, `core:frequency` [P] | Metadata only |
| TorchSig | github.com/TorchDSP/torchsig; MIT [P] | 2.2.0, 31 Aug 2026; the PyPI history starts in April 2026 [P] ([PyPI](https://pypi.org/project/torchsig/)) | Synthetic RF dataset generation with impairments, spectrograms, PyTorch pipeline [P]. Useful for synthetic ground-truth benchmarks | Closed-set human modulations |
| gr-satellites | github.com/daniestevez/gr-satellites; GPL-3.0 (badge) [P] | Active, on conda-forge [P] | Reference implementations of many real deframers and FEC chains (CCSDS and others); a good "known format" negative control | Known formats only |
| SatDump | github.com/SatDump/SatDump; GPL-3.0 [P] | Active [P] | Similar: a library of known satellite decoders | Known formats only |

- Real-world example of the full chain on an authored "alien" signal: the A Sign in Space (ASIS) raw recordings are public as SigMF. The Green Bank Telescope (GBT) file is `ci8` at 1 MS/s, centre 8,410,272,812.5 Hz, start 2023-05-24T19:11:00Z, about 4.5 GB. Each of the ATA-X, ATA-Y, Medicina and GBT files is 4.2–4.8 GB, licensed CC-BY-4.0. [P] [Zenodo 7982936](https://zenodo.org/records/7982936)
- A spacecraft example of the layering to expect from an engineered link: ESCAPADE (2025) uses PCM/PSK/PM, a rate-1/2 turbo code, CCSDS TM frames and Space Packets, decoded by Daniel Estévez from an ATA recording. [S] (search summary of his Bluesky feed)

### Inferences
- A Stargazer physical-layer module should be implemented as estimators (PSD → bandwidth; SCF/cyclic autocorrelation → α lines → candidate symbol rates; x^M lines → carrier and modulation order; CMA → constellation; GMM+BIC → alphabet size). Use GNU Radio and liquid-dsp as engines. Treat AMC classifiers as an optional prior.
- Input requirements, rule of thumb [R]: cyclic features need roughly 10³–10⁵ symbols at moderate SNR; cumulants roughly 10³ symbols; CNN AMC is trained on 128–1024-sample windows and is unreliable below about 0 dB SNR.
- Licensing: GPL-3.0 tools (GNU Radio 3, URH, inspectrum, SigDigger) would force a GPL toolkit if linked. liquid-dsp (MIT), the GR4 core (MIT, [S]) and TorchSig (MIT) give permissive options.

### Gaps
- Could not read gnuradio.org directly, so the GR4 status rests on [S]. No GitHub commit dates (no API access) for inspectrum, SigDigger, Netzob and similar.
- Did not retrieve a quantitative comparison of open blind symbol-rate estimators on astronomy-like narrowband data.
- The ASIS physical-layer modulation (how the payload sat on the TGO carrier) could not be confirmed: destevez.net and asignin.space were unreachable and search snippets were silent.

---

## 2. Bitstream analysis: frame length and sync words, blind recognition of error-correcting codes, scramblers/LFSRs, interleavers, field boundaries

### Takeaway
There is a coherent algebraic toolbox. Each step tests a candidate parameter by a statistic that collapses when the parameter is right:
- frame period: autocorrelation or recurrence;
- block-code length: rank deficiency of reshaped bit matrices, or low-weight parity checks;
- LFSR: linear complexity via Berlekamp–Massey;
- interleaver: parity-check structure.

All of these assume *human-style linear codes over GF(2)*. Most need many frames (tens to thousands of codewords) and a modest bit-error rate (BER). In 2026 the research frontier is polar codes and subspace-based identification. Open-source code is scarce; most of it exists only inside papers.

### Cited Findings
- **Frame and sync detection** (ESTABLISHED) [R]
  - Compute the bitwise autocorrelation A(τ) = Σ (1−2b_i)(1−2b_{i+τ}), or the XOR-match count. A fixed sync word repeated every N bits produces peaks at τ = kN.
  - Then fold the stream at period N and compute per-column bias |p₁−0.5|. Constant columns form the sync word or fixed header; near-0.5 columns carry the payload.
  - A refinement for a drifting period is a recurrence plot or matrix profile on bits (section 5).
  - Data need: a few frames for a deterministic sync word. Failures: variable-length frames, bit slips, and scramblers that cover the sync word.
- **Block-code length and rank method** (ESTABLISHED in the literature) [R]
  - Reshape the bit stream into an m×n' matrix. For n' equal to a multiple of the true codeword length n, aligned, the GF(2) rank falls below n'. A k/n' deficit pattern recurs at multiples of n.
  - Gaussian elimination with noise-tolerant thresholds is in the Sicot, Houcke and Barbier, Valembois, Cluzeau and Marazin line of work.
  - Needs at least n' rows (roughly 2n'–5n'). Breaks down as BER rises past about 10⁻³–10⁻².
- **Low-weight parity-check search** (RESEARCH/ESTABLISHED)
  - Tixier (Inria, 2015) reconstructs an unknown block interleaver and recovers a convolutional code from M noisy interleaved codewords, without assuming the interleaver's structure.
  - Her method finds low-weight parity-check equations and says this is "more efficient than the methods calculating the rank of matrices, in particular when the data are noisy". Noise affects only the parity-check search stage. [P] [arXiv:1501.03715](https://arxiv.org/abs/1501.03715)
- **Deep-learning code recognition over candidate sets** (RESEARCH)
  - Dehdashtian, Hashemi & Salehkaleybar (2020) identify code parameters for LDPC, convolutional, turbo and polar codes over a candidate set, robustly to multipath, without SNR knowledge.
  - Their framing: "second-category" methods recover parameters only from a candidate set. Example candidate set: 802.11n's four LDPC rates {1/2, 2/3, 3/4, 5/6} at n = 648. [P] [arXiv:2009.07774](https://arxiv.org/abs/2009.07774)
- **2026 state of the art in blind code recognition** (RESEARCH)
  - Polar-code blind recognition by SCL decoding exploits the fact that frozen bits favour 0 while information bits are equiprobable in source-side log-likelihood ratios (LLRs), over multiple received vectors. [P] [Tu et al., arXiv:2605.13331 (May 2026)](https://arxiv.org/abs/2605.13331)
  - A hypothesis-testing analysis of the same problem. [P] [arXiv:2606.17705 (Jun 2026)](https://arxiv.org/abs/2606.17705)
  - A critical-set simplified version. [P] [arXiv:2607.17480 (Jul 2026)](https://arxiv.org/abs/2607.17480)
  - A "minimum denoised subspace discrepancy decoder" for code identification on the binary symmetric channel, with theoretical guarantees. The authors note that most existing methods rely on special code structure, are costly otherwise, and largely lack analytical guarantees. [P] [Singh, Krishnan & Yardi, arXiv:2601.15903 (Jan 2026)](https://arxiv.org/abs/2601.15903)
- **Scramblers and LFSRs**
  - Massey's 1969 paper "Shift-register synthesis and BCH decoding" gives the Berlekamp–Massey algorithm. [P] [doi:10.1109/TIT.1969.1054260](https://doi.org/10.1109/TIT.1969.1054260)
  - It finds the shortest LFSR generating a sequence in O(N²) and needs 2L bits for linear complexity L. [R]
  - For an additive (synchronous) scrambler, apply it to known-plaintext regions (for example idle or zero fill). A random sequence has linear complexity ≈ N/2, so a sharp plateau at L ≪ N/2 is the detection signature. [R]
  - Multiplicative (self-synchronising) scramblers: search for a sparse polynomial h(x) such that the descrambled output has minimum entropy or maximum bias. [R]
  - Failure modes: a single bit error breaks Berlekamp–Massey, so noise-tolerant variants or fast-correlation attacks are needed. Nonlinear scramblers defeat it.
- **Field boundaries in single messages**: see NEMESYS in section 3, which works on individual messages using bit-congruence changes.
- **Cautionary example of these tools used on a designed message** (SPECULATIVE)
  - A 2026 Zenodo preprint applies "Berlekamp–Massey LFSR analysis", finite-field enumeration, elliptic-curve point counting and "exhaustive parameter space exploration" to the ASIS payload.
  - It claims a self-referential "algebraic quine" over GF(625) with a chance probability of "approximately 5 × 10⁻¹⁹". [P] [Lacoche 2026, Zenodo 19218629](https://zenodo.org/records/19218629)
  - It is self-published and not peer-reviewed, and the project's public decode (section 4) does not support it. Treat it as a textbook garden-of-forking-paths case (section 6).

### Inferences
- Pipeline order for a suspected digital stream:
  1. hard-decision bits;
  2. autocorrelation and folding for period and sync;
  3. per-column bias map;
  4. Berlekamp–Massey linear complexity on sliding windows;
  5. rank/parity-check scans for (n,k) and convolutional rate;
  6. interleaver search;
  7. decode, and confirm by a syndrome or CRC pass rate.

  Every step yields a candidate *plus* a likelihood, which feeds the validation framework in section 6.
- An off-the-shelf open library for blind FEC recognition was not found. Stargazer would need to implement the rank and parity-check scans itself; liquid-dsp supplies the FEC decoders to *confirm* a hypothesis.
- For a non-human sender, GF(2)-linear coding is an assumption. Rank and Berlekamp–Massey tests still serve as cheap "is there linear redundancy at all" probes.

### Gaps
- Could not retrieve full texts of the classic survey papers on blind channel-code recognition (Semantic Scholar and OpenAlex rate-limited). Specific sample-size and BER thresholds above are [R].
- Found no maintained open-source blind interleaver, scrambler or code-recognition package. Absence on PyPI and arXiv code links only; GitHub search was not possible.

---

## 3. Reverse engineering unknown binary formats and network protocols; grammar inference

### Takeaway
Protocol reverse engineering (PRE) splits into two families:
- **trace-based**: many messages, alignment and clustering. Netzob, Discoverer, NetPlier, NEMESYS/NEMETYL.
- **execution-based**: needs the program binary. Polyglot, BinPRE, SEMPRE.

For a single message of unknown purpose only a few pieces transfer:
- NEMESYS-style intrinsic segmentation;
- entropy and byte-pattern visualisation (binwalk, Veles, scurve);
- grammar compression (Sequitur, Re-Pair) as a hierarchy probe.

Kaitai Struct is the right target language for *writing down* a recovered format, not for discovering one.

### Cited Findings
- **Polyglot** (CCS 2007), Caballero, Yin, Liang, Song (≈434 citations). It uses dynamic binary analysis ("shadowing"): it watches how the implementation processes a message. It notes that trace clustering is "limited by the lack of semantic information on network traces". **Requires an executable.** [P] [doi:10.1145/1315245.1315286](https://doi.org/10.1145/1315245.1315286)
- **Discoverer** (Cui, Kannan, Wang; USENIX Security 2007): automatic PRE from network traces by tokenisation, recursive clustering and type-based merging. Needs many messages. [S] [Microsoft Research page](https://www.microsoft.com/en-us/research/publication/discoverer-automatic-protocol-reverse-engineering-from-network-traces/); details [R]
- **NetPlier** (NDSS 2021, not USENIX), Ye, Zhang, Wang, Zhang, Xu (Purdue). doi:10.14722/ndss.2021.24531.
  - Multiple sequence alignment, then a probabilistic model whose posterior marginals pick the *keyword* (message-type) field.
  - On 10 protocols: 100% homogeneity and 97.9% completeness, versus about 92% and 52.3% for the best prior work. Case studies: Google Nest and a malware command-and-control protocol.
  - A 2023 critique shows it mis-merges and mis-splits fields (OSDP example). [S] [NDSS page](https://www.ndss-symposium.org/ndss-paper/netplier-probabilistic-network-protocol-reverse-engineering-from-message-traces); critique [arXiv:2305.11781](https://arxiv.org/abs/2305.11781)
- **NEMESYS** (Kleber, Kopp, Kargl; USENIX WOOT 2018) "derives field boundaries from the distribution of value changes throughout individual messages … exploits the intrinsic features of structure which are contained within each single message". The repository also contains:
  - NEMETYL message-type identification (INFOCOM 2020, [arXiv:2002.03391](https://arxiv.org/abs/2002.03391));
  - field-type classification by segment clustering;
  - PCA-refined segmentation (IEEE CNS 2022).

  Intended to be used as a Python library. MIT licence (2018). [P] [README](https://raw.githubusercontent.com/vs-uulm/nemesys/master/README.md)
  - Algorithm [R]: bit congruence between consecutive bytes, then a Gaussian-smoothed delta, with boundaries at inflection points.
  - This is **the PRE method that transfers best to a single unknown message**, though it assumes byte alignment.
- **Netzob**: GPLv3, Python 3 library. It infers message format ("vocabulary") and state machine ("grammar") using "bio-informatic and grammatical inferring algorithms", and can then generate and parse traffic. [P] [README](https://raw.githubusercontent.com/netzob/netzob/master/README.rst)
  - Last PyPI release 2.0.0, 8 Jan 2023, so effectively dormant on PyPI. [P] [PyPI](https://pypi.org/project/netzob/)
- **2024–2026 PRE research**
  - BinPRE (2024): binary-analysis field inference. [P] [arXiv:2409.01994](https://arxiv.org/abs/2409.01994)
  - Automatic protocol state-machine inference using fuzzy-DBSCAN format clustering and Needleman–Wunsch/K-Medoids session clustering (2024). [P] [arXiv:2412.02540](https://arxiv.org/abs/2412.02540)
  - **SEMPRE (22 Oct 2026)**: "recursive LLM reasoning grounded in fine-grained dynamic program analysis". Taint and loop traces guide an LLM tree-of-thought that slices messages into hierarchical fields, including nested sub-payloads. **Requires program execution.** [P] [arXiv:2610.31844](https://arxiv.org/abs/2610.31844)
- **Curated list** of PRE tools: [github.com/techge/PRE-list](https://www.github.com/techge/PRE-list) [S]
- **binwalk v3** is a Rust rewrite. Cargo.toml shows version 3.1.1, edition 2024, MIT licence. It identifies and extracts embedded files and offers "entropy analysis … to identify unknown compression or encryption". [P] ([README](https://raw.githubusercontent.com/ReFirmLabs/binwalk/master/README.md), [Cargo.toml](https://raw.githubusercontent.com/ReFirmLabs/binwalk/master/Cargo.toml), LICENSE MIT 2024)
  - The PyPI `binwalk` 2.1.0 (2015) is stale. [P] [PyPI](https://pypi.org/project/binwalk/)
- **Veles** (CodiLime): binary visualisation, "statistical visualizations let you find the important bits in a sea of binary data". It still points to IRC on freenode, so it looks unmaintained. [P] [README](https://raw.githubusercontent.com/codilime/veles/master/README.md)
- **scurve** (Cortesi): space-filling-curve byte visualisation; README read. [P] [repo](https://raw.githubusercontent.com/cortesi/scurve/master/README.md)
- **Kaitai Struct**: declarative format language (`.ksy`) compiled to parsers in many languages. Runtime `kaitaistruct` 0.11 (8 Sep 2025, MIT); the compiler is GPL-3.0. [P] ([PyPI](https://pypi.org/project/kaitaistruct/), compiler LICENSE)
  - Role in Stargazer: a machine-checkable *output* format for a recovered layout, which can be re-parsed to test it.
- **Grammar inference**
  - SEQUITUR (Nevill-Manning & Witten, 1997) infers a hierarchical grammar by replacing repeated phrases with rules, enforced by digram-uniqueness and rule-utility constraints. Linear time and space, incremental, about 50,000 symbols per second in their implementation. [P] [arXiv:cs/9709102](https://arxiv.org/abs/cs/9709102)
  - Re-Pair (Larsson & Moffat 2000, "Off-line dictionary-based compression") repeatedly replaces the most frequent pair. [P] [doi:10.1109/5.892708](https://doi.org/10.1109/5.892708)
  - Zenil & Uthamacumaran (2026) use "canonically relabelled Re-Pair grammar induction" as a structure statistic against matched nulls. On the Phaistos Disc it "saves 29 symbols against a null". [P] [arXiv:2610.08861](https://arxiv.org/abs/2610.08861)
  - Time-series grammar induction for anomaly detection (Sequitur on SAX words). [P] [Gao, Lin, Brif, arXiv:2001.11102](https://arxiv.org/abs/2001.11102)
  - `saxpy` (GPL-2.0) 2.0.1, 23 Jul 2026. [P] [PyPI](https://pypi.org/project/saxpy/)

### Inferences
- **Transfers to one unknown message:**
  - NEMESYS segmentation, if the message is byte-structured;
  - sliding-window entropy and bigram maps (binwalk, scurve style);
  - Sequitur/Re-Pair grammar size versus shuffled nulls;
  - fold-and-bias maps (section 2).
- **Assumes many samples:** Netzob, Discoverer, NetPlier, NEMETYL, URH field inference.
- **Assumes an executable:** Polyglot, BinPRE, SEMPRE and other taint-based tools. These are irrelevant for astronomical data except as inspiration for an "interventional" analogue: Zenil 2026 uses "executable controls" on known generators (section 6).
- A repeating beacon would supply "many messages" naturally. Alignment-based PRE then becomes usable once frames are cut (section 2).

### Gaps
- Could not read the Discoverer, Netzob (AsiaCCS 2014) or NEMESYS full papers for quantitative accuracy figures (rate limits; usenix.org not on the allowlist).
- No 2025–2026 LLM-based PRE tool found that works on messages alone, without a binary.

---

## 4. Discovering dimensionality and layout (raster images, arrays): Arecibo, Cosmic Call, A Sign in Space, Zenil et al., Heller's SETI Decrypt Challenge

### Takeaway
Every designed message so far has hidden its layout in a semiprime or composite length, or in repeating page frames:
- Arecibo: 1,679 = 23×73.
- Heller's challenge: 1,902,341 = 7×359×757.
- Cosmic Call: 127×127 pages.
- ASIS: a 256×256 bitmap with 10-byte header and footer.

Width recovery is solved for such cases by factorisation plus autocorrelation or complexity scans. Zenil et al. show algorithmic-complexity (Block Decomposition Method, BDM) scans recover Arecibo's width. Their own 2026 follow-up shows the *raw* BDM minimum is at the wrong width (11) and needs shape normalisation plus shuffle calibration. The hard part of ASIS was not layout but a hidden second layer: a reversible cellular automaton run for 6,625 generations. Citizen scientists needed about 10 days to extract the bitmap and about 1 year to decode it.

### Cited Findings

**Arecibo message**
- 1,679 bits "meant to be arranged into 23 columns of 73 rows". Arranging it as 23 rows × 73 columns instead scrambles the picture. [P] [Zenil et al., arXiv:2303.16045](https://arxiv.org/abs/2303.16045)
- Sent 16 Nov 1974 [P] ([arXiv:2411.09790](https://arxiv.org/abs/2411.09790)), at about 2,380 MHz using frequency-shift keying, about 10 bit/s, about 3 minutes [R].

**Zenil, Adams, Abrahão & Ozelim**, "An Optimal, Universal and Agnostic Decoding Method for Message Reconstruction, Bio and Technosignature Detection" (arXiv:2303.16045, v4 10 Jan 2025) [P] ([arXiv](https://arxiv.org/abs/2303.16045)). The author list is Zenil, Adams, Abrahão and Ozelim; the brief's "Hernández-Orozco" is a co-author of the earlier BDM paper.
- Method: reshape the 1-D stream into every candidate (m,n) partition and compute three measures:
  - Shannon (block) entropy;
  - BDM (CTM complexity of 4×4 blocks summed with log multiplicities);
  - zlib length (b64 of zlib output).

  Candidate dimensions are the **downward complexity spikes**.
- Noise:
  - At 3% of bits randomly replaced (≈1.5% actually flipped) the method still finds 23×73.
  - BDM is more robust than compression and entropy. Compression loses sensitivity at about 10% bits flipped.
  - At 60% pixels replaced (≈30% flipped) BDM is still sensitive while the others are about 50% sensitive.
  - zlib "may not perform well" on small strings because of overhead; images were up-scaled 6×6 per bit to help zlib.
- Generalisation test on 102 Caltech-101 images (about 300×200 pixels): 29,511 downward BDM peaks plus eligible neighbours.
  - The correct dimension was the lowest-BDM candidate in about 80% of cases; the rest were blamed on border effects.
  - A depth-6 decision tree separated all peaks, but the authors say it "could not be used to extrapolate to any other image".
  - The correct size is sometimes an immediate neighbour of the peak. Spikes also appear at multiples or divisors of the true width (for example 250 for 500).

**Zenil & Uthamacumaran**, "Reconstruction of Geometric and Structural Aspects of Undeciphered Messages by Information-theoretic Methods" (5 Oct 2026) — a 2026 development [P] ([arXiv:2610.08861](https://arxiv.org/abs/2610.08861))
- New PyBDM calculation on the 1,679-bit Arecibo message over widths 10–100.
  - **The raw BDM minimum is at width 11.** Shape-normalised BDM selects 23.
  - Calibrating each width with 99 bit-count-preserving shuffles gives the largest standard score, z = 10.174, at width 23. Comparing with 199 shuffled full-width scans gives a search-corrected p = 0.005.
  - The authors call this an estimator-sensitivity analysis, "not a preregistered replication": normalisation was chosen after an exploratory scan.
- Noise table (exact bit inversion, 99 replicates):

  | Bits inverted | Width 23 ranked first | Width 23 in top 3 | Median rank |
  |---|---|---|---|
  | 1% | 47% | 100% | — |
  | 3% | 54% | 93% | — |
  | 5% | 42% | 85% | — |
  | 10% | 18% | 44% | 5 |
  | 15% | 9% | 19% | 15 |

- Extended to the Phaistos Disc, khipu, Rongorongo, Indus and Voynich with matched nulls. For example, Phaistos retains 0.714 bits of excess adjacent-sign mutual information after controls; the matched null mean was already 2.522 bits, so the uncorrected figure was 3.236. The abstract stresses "distinguishing recoverable organisation from semantic decipherment".

**Block Decomposition Method (BDM)** (Zenil, Hernández-Orozco, Kiani, Soler-Toscano, Rueda-Toicen; arXiv:1609.00110) [P] ([arXiv](https://arxiv.org/abs/1609.00110))
- Extends the Coding Theorem Method (CTM) by decomposing an object into blocks whose complexity is taken from small Turing-machine output frequencies.
- It "performs like Shannon entropy" when block-level CTM information runs out, so it falls back to entropy on large objects.
- PyBDM: 1D and 2D binary arrays and 1D non-binary arrays; the alphabet must be mapped to integers; boundary leftovers are handled by the chosen partition algorithm; includes perturbation analysis. MIT, but only release 0.1.0 (22 Sep 2019). [P] ([README](https://raw.githubusercontent.com/sztal/pybdm/master/README.rst), [PyPI](https://pypi.org/project/pybdm/))

**Heller's SETI Decrypt Challenge** (2016) [P] ([Heller, arXiv:1706.00653](https://arxiv.org/abs/1706.00653); submitted to the International Journal of Astrobiology, 19 Oct 2017)
- Posted to Twitter and Facebook on 26 April 2016; submission deadline 3 June 2016.
- 1,902,341 bits = 7 × 359 × 757, all primes: seven pages of 359 (width) × 757 (height) pixels.
  - Page 1 is a configuration frame defining the page size.
  - Page 2 is a count from 0 to 756 in *little-endian* binary.
  - Page 3 lists the first 757 primes (largest 5,749).
  - Later pages hold images, with headers giving spatial and temporal yardsticks.
- 305 total submissions (solutions, hint requests and so on) and **66 correct submissions from 71 people**. The first solutions arrived the next day.
- About half were derived fully independently; the other half benefited from online discussion and a Reddit spoiler. Most submissions after the spoiler were correct.
- The frame size was recognised by most respondents. Interpreting the header yardsticks as physical units was where people diverged.

**Cosmic Call (Dutil–Dumas message, 1999; DDM2 2003)** [S] ([Wikipedia](https://en.wikipedia.org/wiki/Cosmic_Call); [plover.com analysis series](https://blog.plover.com/aliens/dd/p00.html))
- 23 pages of 127×127 bitmaps (127 is prime). The same bits can be read as 127 × 2,921.
- Content is mostly 5×7 = 35-bit glyphs, filtered to differ in enough bit positions to tolerate errors and spread across lines against burst errors.
- Each page has a black-then-white border frame. Page-number glyphs tolerate changes to any 4 of 9 pixels.
- The scientific part was sent three times at 100 bit/s, the public part once at 2,000 bit/s.
- Known content error: neutron mass given as 1.67392 instead of 1.67492, corrected in DDM2.
- No formal FEC code (such as Reed–Solomon) was confirmed. Redundancy is by design of glyphs and layout.

**A Sign in Space (ASIS)**
- Transmitted 24 May 2023 at 19:00 UTC by ESA's ExoMars Trace Gas Orbiter. Received by the Allen Telescope Array, the Green Bank Telescope and Medicina; some accounts add the VLA. [S] ([SETI Institute](https://www.seti.org/news/first-contact-global-team-simulates-message-from-extraterrestrial-intelligence-to-earth/); [sci.news](https://www.sci.news/astronomy/extraterrestrial-message-11938.html))
- Raw data published as SigMF. [P] [Zenodo 7982936](https://zenodo.org/records/7982936)
- **Extraction:** "a community of 5000 citizen scientists extracted the signal within 10 days". [S] ([ESA](https://www.esa.int/ESA_Multimedia/Images/2024/10/Alien_signal_decoded); [ESA blog 23 Oct 2024](https://blogs.esa.int/to-mars-and-back/2024/10/23/father-daughter-team-decodes-alien-signal/))
- **Decoding:** solved by Ken and Keli Chaffin (father and daughter). The project site gives 7 June 2024, 21:50 CEST, with names changed for privacy, which makes it about a year after transmission. The SETI Institute headline calls it "A Year-Long Journey". [S] ([asignin.space 12 Jun 2024](https://asignin.space/12-june-2024-finally-the-message-we-designed-for-a-sign-in-space-has-been-successfully-decoded/); [SETI Institute](https://www.seti.org/news/decoding-the-message-from-a-sign-in-space-a-year-long-journey/); [CNN 11 Nov 2024](https://www.cnn.com/2024/11/11/science/seti-alien-signal-decoded-sign-in-space))
- **Format:**
  - The payload is a 256×256 one-bit bitmap ("starmap"). The file `data17.bin` has a 10-byte header, a 65,536-bit body (MSB-first) with **exactly 625 set pixels**, and a 10-byte footer. [P] [Lacoche, Zenodo 21032854](https://zenodo.org/records/21032854) (methodology PDF)
  - The hidden layer is a **reversible Margolus 2×2 block cellular automaton** (single-rotation rule, counter-clockwise, phase = generation mod 2, toroidal) run for **6,625 generations**. [P] (same PDF, adopting the community decode)
  - Ken Chaffin "suspected the data came from a cellular automaton" and ran simulations "for hours and days on end". [S] (ESA)
  - The output shows five amino acids, drawn as blocks of 1, 6, 7 or 8 pixels standing for the atomic numbers of H, C, N and O. [S] (search summaries of asignin.space and ESA)
- 2024 follow-on interpretation study: 285 active interpreters. [P] [de Paulis & Jensvold, Zenodo 18156989 (Jan 2026)](https://zenodo.org/records/18156989)
- **2026 developments**
  - "A Sign in Space part 2" (arXiv 6 Oct 2026) is about a global telescope-network detection exercise using an ESA spacecraft, NVIDIA Holoscan and IAA post-detection protocols. It is not a new decoding puzzle description. [P] [de Paulis, Melis, Estévez, arXiv:2610.10587](https://arxiv.org/abs/2610.10587)
  - An independent 2026 preprint gives a "ledger" of what a blind pipeline can determine from ASIS. It re-verifies exact reversibility under the Margolus rule and denoises by the automaton's own dynamics.
    - It shows {1,6,7,8} → {H,C,N,O} is the unique admissible blob-size alphabet *given chemistry priors*.
    - Its own MDL test *retracts* a claimed third dimension: a 2D lattice fits with fewer parameters.
    - It concludes that "self-consistency is not external validation". [P] [Zenodo 21032854](https://zenodo.org/records/21032854) (self-published preprint)

**Autocorrelation row-length detection** (ESTABLISHED) [R]
- For a raster of width W sent row by row, the bit autocorrelation has peaks at τ = W, 2W, …, because vertically adjacent pixels are correlated.
- Variants:
  - 2-D FFT of the folded image versus width;
  - minimising the column-wise XOR between consecutive rows;
  - matrix-profile periodicity.
- Works when images have vertical coherence. Fails for text-like or pre-compressed content, which has no inter-row correlation.
- Prime-factorising the length gives a tiny candidate set (Arecibo: {23,73}) when the sender intends it.

### Inferences
- A Stargazer "layout" module should run several detectors and require agreement:
  1. factorise N, and N minus small header/footer sizes (ASIS needed the 10+10-byte trim);
  2. autocorrelation and folding scans over W;
  3. normalised, shuffle-calibrated BDM or compression scans with a search-corrected p-value, following Zenil 2026;
  4. 2-D autocorrelation and FFT on the best candidates;
  5. a recursion test for 3-D (pages: Heller's 7, Cosmic Call's 23).
- Lessons from the episodes:
  - Layout is fast to find (Heller: one day; ASIS extraction: about 10 days).
  - Hidden *transforms* are slow (ASIS cellular automaton: about 1 year).
  - The interpretation of units and headers is where solvers diverge (Heller).
  - So the toolkit needs a "transform hypothesis library": XOR with an LFSR, cellular-automaton rules, bit-order reversals (Heller's little-endian), Gray code, run-length coding. Each should be scored by complexity drop.

### Gaps
- The exact ASIS physical layer (modulation and framing on the TGO downlink) and the extraction tooling the community used were not confirmed.
- The Cosmic Call primary technical report (Dutil & Dumas) was not read.
- Exact Arecibo radio parameters are [R].

---

## 5. Algorithmic-information and compression probes; time-series and audio tokenisation

### Takeaway
Compression length is the most general structure score and the natural currency for comparing decoding hypotheses (MDL). It is biased on short inputs because of compressor overhead and finite windows: zlib needs kilobytes, LLMs a context of roughly 2 kB tokens, and model size must be counted. BDM/CTM targets very short strings but degrades to entropy at scale. For continuous signals, the matrix profile (STUMPY), SAX plus Sequitur, HMMs, ε-machines (CSSR) and change-point detection (ruptures) are mature, permissively licensed tools for turning a waveform into tokens. HuBERT-style discrete units are powerful but speech-specific.

### Cited Findings

**Compression and complexity measures**
- **NCD** (Cilibrasi & Vitányi): NCD(x,y) = [C(xy) − min(C(x),C(y))] / max(C(x),C(y)), followed by hierarchical clustering. It is a similarity metric that "approximates universality" for a "normal compressor"; the normalised information distance it approximates is uncomputable. [P] ([arXiv:cs/0312044](https://arxiv.org/abs/cs/0312044); IEEE Trans. IT 2005, [doi:10.1109/TIT.2005.844059](https://doi.org/10.1109/TIT.2005.844059))
  - Known limits [R]: the compressor window (gzip 32 KB, bzip2 900 KB blocks) breaks NCD when |x|+|y| exceeds it; small inputs are dominated by header overhead.
  - NCD on large malware files was studied. [P] [arXiv:1509.00689](https://arxiv.org/abs/1509.00689)
- **Neural NCD** (Hurwitz, Nicholas, Raff, 2024) turns LLMs into lossless compressors for NCD and finds "classification accuracy is not predictable by compression rate alone". So a better compressor is not automatically a better structure probe. [P] [arXiv:2410.15280](https://arxiv.org/abs/2410.15280)
- **"Language Modeling Is Compression"** (Delétang et al., 2023) [P] ([arXiv:2309.10668](https://arxiv.org/abs/2309.10668))
  - Chinchilla 70B compresses ImageNet patches to 43.4% and LibriSpeech to 16.4% of raw size, beating PNG (58.5%) and FLAC (30.3%).
  - Transformers "can only compress a few kilobytes"; contexts of 2048 tokens, versus gzip's 32 KB and LZMA2's virtually unlimited window. Datasets were chunked into 2048-byte sequences.
  - With the **adjusted** rate (model parameters counted in the compressed output), "scaling beyond a certain point will deteriorate the compression performance".
  - Any compressor such as gzip can be turned into a conditional generative model.
- **MDL / two-part codes** (ESTABLISHED) [R]
  - Score hypothesis H by L(H) + L(data | H) (Rissanen). In practice L(data | H) comes from arithmetic coding under a fitted model (for example adaptive context or n-gram per field) and L(H) from a parameter-count code.
  - Zenil 2026 and Lacoche 2026 both use MDL to reject an extra dimension [P] (sections 4 and 6).
- **Context-mixing compressors**: paq8px is GPL-licensed and actively versioned (archive extension `.paq8px217` in its README). It has strong specialised models (image and audio) that can bias the comparison. [P] [README](https://raw.githubusercontent.com/hxim/paq8px/master/README.md)
- **Lempel–Ziv** [R]: LZ76 complexity (Kaspar–Schuster normalisation c·log₂n/n) is a fast, parameter-free probe. It is biased upward for n below about 10³ symbols; compare against shuffled surrogates of the same length.

**Time-series and audio tokenisation**
- **SAX** (Lin, Keogh, Wei, Lonardi 2007): z-normalise, Piecewise Aggregate Approximation (PAA) with w segments, then map to an a-letter alphabet with Gaussian breakpoints. Its distance lower-bounds Euclidean distance. [P] ([DMKD, doi:10.1007/s10618-007-0064-z](https://doi.org/10.1007/s10618-007-0064-z); properties noted in [arXiv:2010.00730](https://arxiv.org/abs/2010.00730))
  - SAX misses within-segment trends. [P] [arXiv:1905.00421](https://arxiv.org/abs/1905.00421)
  - Libraries: `pyts` 0.14.1 (BSD-3, 21 Oct 2026), `saxpy` (GPL-2.0), `tslearn` 0.10.0 (BSD-2, 22 Oct 2026), `sktime` 1.3.0 (BSD-3, 8 Nov 2026). [P] (PyPI JSON)
- **Matrix profile and STUMPY**: STUMPY (JOSS 2019, [doi:10.21105/joss.01504](https://doi.org/10.21105/joss.01504)) covers motif and discord discovery, semantic segmentation (FLUSS), snippets, Dask-distributed and GPU versions. BSD-3, release 1.14.1 on 8 Feb 2026. [P] ([README](https://raw.githubusercontent.com/stumpy-dev/stumpy/main/README.rst); [PyPI](https://pypi.org/project/stumpy/))
  - The window m must be chosen, or scanned with pan-matrix-profile. Cost is O(n²) [R].
  - Robustness to noise is uncertain; the conventional pre-filtering cannot be used unsupervised. [P] [Hehir & Smeaton, arXiv:2306.10151](https://arxiv.org/abs/2306.10151)
- **ε-machines and CSSR** (Shalizi & Shalizi 2004): CSSR builds the minimal optimal predictive HMM (causal states) from a discrete sequence; at least comparable to cross-validated HMMs. [P] [arXiv:cs/0406011](https://arxiv.org/abs/cs/0406011)
  - **2026 astronomy application:** ε-machines on repeating-FRB waiting times. The data were symbolised into k = 2–5 bins, with history length L = 5 and permutation surrogates. Validation used synthetic sequences of N = 100–2,000 symbols (200 realisations each). Sequences: 1,613, 877 and 193 waiting times. Both FAST sources show about 1 bit of temporal memory (p ≤ 0.01); FRB 20220912A is memoryless. [P] [Kimpson & O'Leary, arXiv:2607.14421](https://arxiv.org/abs/2607.14421)
  - Earlier: ε-machines as agnostic planetary-complexity measures. [P] [arXiv:2202.03699](https://arxiv.org/abs/2202.03699)
- **HMMs and change points**: `hmmlearn` 0.3.3 (BSD, Oct 2024; limited-maintenance mode [R]); `ruptures` 1.1.10 (BSD-2, 27 Nov 2025) with PELT, binary segmentation, window and dynamic-programming searches [P] (PyPI JSON; algorithm list [R]).
- **Zero-resource acoustic units**
  - HuBERT: offline k-means targets with masked prediction. The fairseq release has HuBERT Base (~95M parameters, LibriSpeech 960 h) with "L9 km500" units (500 k-means clusters on layer-9 features). Label rate is 50 Hz for HuBERT features and 100 Hz for MFCC. [P] ([arXiv:2106.07447](https://arxiv.org/abs/2106.07447); [fairseq README](https://raw.githubusercontent.com/facebookresearch/fairseq/main/examples/hubert/README.md))
  - ZeroSpeech 2021: CPC, then k-means, then an LSTM/BERT language model, giving above-chance results on phonetic, lexical, syntactic and semantic probes. [P] [arXiv:2011.11588](https://arxiv.org/abs/2011.11588)
  - Unit discovery by VQ-VAE / VQ-CPC [P] ([arXiv:2005.09409](https://arxiv.org/abs/2005.09409)) and Bayesian subspace HMMs [P] ([arXiv:1904.03876](https://arxiv.org/abs/1904.03876)).
  - Spoken-term discovery by aligning discrete-unit sequences (DUSTED, 2024). [P] [arXiv:2408.14390](https://arxiv.org/abs/2408.14390)
  - 2026: HuBERT-style masked unit prediction is still the template, for example for the Interspeech 2026 "Unsupervised Speech in the Wild" challenge (BiMamba2, 47.88M parameters, 250 h, 67 languages). [P] [arXiv:2609.28758](https://arxiv.org/abs/2609.28758)
  - The completeness of k-means units is debated. [P] [arXiv:2409.06109](https://arxiv.org/abs/2409.06109)

### Inferences
- Use compression as a *relative* score: the same compressor, same length, compared against shuffled or surrogate baselines. Never use it as an absolute verdict on short data.
  - Below about 1 kbit, prefer BDM/CTM or LZ76 with surrogate calibration.
  - From 1 kbit to 10 Mbit, use an ensemble of zlib, xz and paq8px (or cmix).
  - For LLM compressors, account for model size and the context limit (about 2 kB).
- For tokenising an unknown waveform, a defensible open stack is:
  - change points (ruptures) and matrix-profile motifs (STUMPY) to propose segment boundaries and units;
  - SAX or k-means/GMM on segment features to discretise;
  - CSSR or HMM to model the resulting symbol process;
  - Sequitur/Re-Pair for hierarchy.

  Choose the alphabet size and window by MDL, as in the FRB paper's k-sweep.
- HuBERT/wav2vec units carry speech priors. For a non-speech astrophysical signal, self-supervised training *on the target data itself* (or RF foundation models, section 1) is the closer analogue. That is RESEARCH, not established.

### Gaps
- The cmix repository could not be found at the path tried. The licence and status of cmix (and of a pip-installable CSSR package referenced in Kimpson 2026) are unconfirmed.
- No benchmark comparing these tokenisers on synthetic "alien" signals was found.

---

## 6. Validation: telling real structure from analyst-imposed structure; false-decoding episodes

### Takeaway
Every recovered structure should be scored against matched surrogates, with the search itself corrected for (look-elsewhere). It should be confirmed by at least one of three things:
- independent methods agreeing;
- held-out predictability or compression gain;
- intrinsic redundancy checks passing (CRC/syndrome, exact reversibility, round-trip re-encoding).

The cautionary cases show how this fails:
- Bible codes: extensive search yields "significant" patterns.
- An ASIS "quine" claim at p ≈ 5×10⁻¹⁹ after exhaustive exploration.
- Zenil's own raw-BDM minimum landing at the wrong width.
- Heller's solvers converging through spoilers rather than independently.

### Cited Findings
- **Shuffle and surrogate calibration with search correction**
  - Zenil & Uthamacumaran 2026 calibrate each width with 99 bit-count-preserving shuffles, then compare the maximum statistic with 199 shuffled full scans: search-corrected p = 0.005 at width 23.
  - They flag that normalisation was chosen after exploration ("not a preregistered replication").
  - Phaistos nulls preserve face and role (4,999 draws; 5,000 for statistical and dictionary measures, 1,000 for CTM-BDM). The excess over the null is reported, not the raw value: 0.714 versus a raw 3.236 bits.
  - Row-alignment tests "yield no comparable evidence after correction". [P] [arXiv:2610.08861](https://arxiv.org/abs/2610.08861)
- **Permutation surrogates for memory:** 200 surrogates and false-discovery-rate adjustment across sources in the FRB ε-machine work; FRB 20220912A came out consistent with no memory. [P] [arXiv:2607.14421](https://arxiv.org/abs/2607.14421)
- **Redundancy and invertibility confirmation**
  - The ASIS decode is confirmed by exact reversibility under the Margolus rule, by denoising through the automaton's own dynamics, and by a round-trip re-projection test.
  - The author nonetheless stresses "self-consistency is not external validation".
  - An MDL test retracted a 3-D reading: a 2D lattice gives the same residual with fewer parameters, with b₃ ≈ 3·b₁ + b₂.
  - Two independent molecule-recognition networks (DECIMER, MolNexTR) agreed on the identities. They were calibrated at 100% on standard depictions but 0% on the raw pixel notation, so they corroborate only the graph-to-identity step. [P] [Zenodo 21032854](https://zenodo.org/records/21032854)
- **False or over-claimed decoding episodes**
  - **ASIS "algebraic quine"** (2026, self-published): exhaustive search over finite fields, LFSRs and elliptic curves. It claims chance probability "≈5 × 10⁻¹⁹; under uniformity assumptions ≈10⁻⁷⁶", with a "null" of only 120 alternative inputs. The fact that the structure "breaks with any single parameter change" is offered as support, but is equally the signature of a tuned fit. [P] [Zenodo 19218629](https://zenodo.org/records/19218629)
    - Lesson: a p-value computed after an unconstrained search over transform families is meaningless without simulating the *whole search* on null data.
    - The same author's later preprint ([Zenodo 21032854](https://zenodo.org/records/21032854)) is far more cautious.
  - **Bible codes**: "Solving the Bible Code Puzzle" (Bar-Hillel, Bar-Natan, Kalai, McKay; Statistical Science 1999) refuted equidistant-letter-sequence claims [P] ([doi:10.1214/ss/1009212243](https://doi.org/10.1214/ss/1009212243)). Their argument [R] is that undisclosed researcher choices (spellings, names, dates) let a similar effect be produced in *War and Peace*.
    - Lesson: freeze the analysis pipeline and its choices before looking, and run it on control texts.
  - **Heller's challenge**: about half the "correct" solutions were not independent (spoilers and online discussion). Convergence of solvers is not independent confirmation unless isolation is enforced. Even the author notes a real message "might not even foresee a single correct interpretation". [P] [arXiv:1706.00653](https://arxiv.org/abs/1706.00653)
  - **Estimator artefacts**: raw BDM picks width 11 for Arecibo. Shape normalisation and calibration are necessary. [P] [arXiv:2610.08861](https://arxiv.org/abs/2610.08861)
  - **Compression versus classification:** better compression does not mean better downstream structure recovery. [P] [arXiv:2410.15280](https://arxiv.org/abs/2410.15280)
  - **Pulsar "LGM-1" (1967)** and the **perytons** (2015, traced to microwave ovens) are the classic cases of apparent artificiality from a natural or local source [R]. They belong mainly to the detection layer covered by another researcher.
- **Synthetic ground truth**
  - TorchSig generates impaired RF datasets [P].
  - Zenil 2023 used 102 Caltech-101 images with known dimensions [P].
  - ASIS is a "deliberately authored alien message with a knowable … intended meaning", now public (SigMF raw data plus the bitmap) [P].
  - Arecibo, the Heller challenge (bitstream published with the paper's description) and Cosmic Call are other known-answer messages.

### Inferences
Validation protocol for Stargazer (proposed):
1. Pre-register the hypothesis family and search ranges.
2. For every statistic, compute a null by re-running *the whole search* on surrogates that preserve length, bit count, and (for time series) the power spectrum (IAAFT) or block structure.
3. Report a search-corrected p-value and the effect size over the null.
4. Require two methods from different families to agree, for example autocorrelation width and calibrated BDM width.
5. Require held-out predictability: fit on half the frames, measure code length on the other half.
6. Require a hard check where possible: CRC, syndrome or parity pass rate, exact invertibility, a Kaitai re-parse round trip.
7. Run the identical pipeline on known-natural controls (pulsars, RFI, scintillating sources) and on known-answer messages (Arecibo, Heller, Cosmic Call, ASIS) as positive controls.
8. Keep human interpreters blind to one another's results, which is the lesson of Heller's challenge.

### Gaps
- The McKay et al. full text was not read (only the Crossref record); the details of the refutation are [R].
- No systematic catalogue of "false decodes" in SETI was found in the literature searched.

---

## 7. Practical demonstrations by amateurs and researchers of decoding real unknown-format artefacts

### Takeaway
Successful amateur decodes used commodity tools (Audacity, Python, cellular-automaton simulators, SDR viewers) plus human pattern recognition. Each had a strong prior: the record cover's instructions, prime-number layouts, the guess of a cellular automaton. The time scale runs from one day (Heller) to about 10 days (ASIS extraction) to about 1 year (ASIS hidden transform).

### Cited Findings
- **Voyager Golden Record images** (1977): Ron Barry decoded the analog-encoded images by following the cover diagram's instructions in Audacity, using a 32-bit, 384 kHz WAV supplied by David Pescovitz; Boing Boing, 5 Sep 2017. [S] [Boing Boing](https://boingboing.net/2017/09/05/how-to-decode-the-images-on-th.html)
  - Follow-ons include:
    - Marc Bäuerle's browser-based decoder of all 116 images, written because details in Barry's write-up were unclear [S] ([GitHub](https://github.com/MarcBaeuerle/Golden-Record));
    - a blog re-decode from the raw recording [S] ([hackcorrelation](https://hackcorrelation.blogspot.com/2017/10/decoding-voyager-golden-record-images.html));
    - the `voyagerimb` Python tool [S] ([GitHub](https://GitHub.com/aizquier/voyagerimb));
    - a University of Luxembourg course project [S] ([uni.lu](https://basv53.uni.lu/projects/07_voyager_decoder.html)).
  - NASA notes the images are copyright-protected. [S] [NASA](https://science.nasa.gov/mission/voyager/golden-record-contents/images/)
  - Format [R]: scan-line analog video in the audio channels, about 512 lines per image, with line sync pulses; the cover etching gives the waveform and the line count.
- **ASIS:** about 5,000 citizen scientists extracted the bitmap from the raw radio data within about 10 days. Ken and Keli Chaffin found the cellular-automaton layer after about a year by running simulations "for hours and days on end". [S] ([ESA](https://www.esa.int/ESA_Multimedia/Images/2024/10/Alien_signal_decoded))
- **Heller's challenge:** first correct solutions the day after posting; 66 correct from 71 people over about 5 weeks. [P] [arXiv:1706.00653](https://arxiv.org/abs/1706.00653)
- **Spacecraft telemetry by amateurs:** Daniel Estévez's gr-satellites (GPL-3.0) and blog decodes of real spacecraft (for example ESCAPADE, rate-1/2 turbo code, CCSDS) show the full SDR → FEC → frame → packet chain done with open tools. [P] for the repository ([README](https://raw.githubusercontent.com/daniestevez/gr-satellites/main/README.md)); [S] for the ESCAPADE detail.

### Inferences
- Stargazer could formalise what these amateurs did by hand into hypothesis libraries plus scoring: measure periods, factorise lengths, try transforms, look at the image.
- The ASIS cellular-automaton layer shows the transform library must include generative and dynamical decoders (cellular-automaton rules, iterated maps), scored by complexity drop and by exact reversibility.

### Gaps
- Could not read the Boing Boing, uni.lu or GitHub write-ups directly (hosts not allowed), so the Voyager technical parameters are [R] or [S].
- No first-hand write-up by the Chaffins was reachable.
- No primary account of the ASIS extraction tooling (which SDR software or scripts the community used) was found.
