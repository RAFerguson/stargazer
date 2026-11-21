# Research phase

This folder holds the research that `docs/PLAN.md` and `docs/INNOVATION_AGENDA.md` are built on.

- `reports/Celestial signal search data and methods.md` is the consolidated research report.
- `research_notes/Celestial signal search data and methods/` holds eight per-topic notes. Each is organised as key
  questions with a takeaway, cited findings with links, inferences and gaps:

  | Note | Covers |
  |---|---|
  | `optical_highenergy_timedomain_archives.md` | Optical, infrared, X-ray, gamma-ray, neutrino and gravitational-wave archives and how to query them, with live tests of the main services |
  | `radio_seti_archives.md` | Public radio SETI data, hit tables and software, with a timing benchmark |
  | `pulsars_frbs_radio_transients.md` | Pulsars, fast radio bursts and long-period radio transients: catalogues, methods, natural behaviour, SETI proposals |
  | `quasar_agn_variability_periodicity.md` | How quasars vary, what happened to every periodicity claim, the statistics, the data |
  | `technosignature_detection_methods.md` | Search strategies, information-theoretic tests, anomaly detection, statistics, contaminants, verification protocols |
  | `decoding_structure_of_unknown_signals.md` | Recovering symbols, framing, layout and coding from an unknown signal or bit stream |
  | `decipherment_and_translation_without_a_key.md` | Ancient scripts, computational decipherment, animal communication, interstellar message tests, language models |
  | `prior_art_and_engineering.md` | Existing projects, reusable software, compute, data volumes, lessons |

## How the research was done, and its limits

The research ran in two rounds from a cloud sandbox.

- **14 November 2026.** The sandbox could reach package registries and public cloud buckets but no paper or
  archive site, and the web-search budget ran out part-way. Three notes (pulsars, quasars, prior art) were written
  from search summaries and package metadata, and say so finding by finding.
- **21 November 2026.** Access was opened to arXiv, three paper indexes and the core astronomy archives. The other
  five notes were written or completed then, from primary sources and from live requests to the archive services.
  Publisher sites and several archives stayed out of reach; each note lists what it could not read.

Every note explains its own tags at the top. In all of them a finding is marked as read in a primary source, taken
from a search summary, or recalled and not re-checked.

## Citation check

The notes cite 516 distinct arXiv identifiers. `citation_check/arxiv_citation_check.py` looks each one up through
the arXiv API, records the paper's title, first author and date, and counts the vocabulary its title and abstract
share with the line that cites it, so that weak matches can be read by a person. The result is in
`citation_check/citations.tsv`.

Status on 21 November 2026:

- **286 identifiers checked**: every identifier in the notes as they stood that morning (the three notes from the
  first round and the first half of the archive note). All 286 exist. The 31 citations with the weakest vocabulary
  match were read by hand, and none was misattributed.
- **230 identifiers not yet checked.** These were added later the same day, in the notes written with direct arXiv
  access, whose authors tag each as read in the source. The check stalled when arXiv began rate-limiting the
  session. The script is resumable: rerunning it fills in the rest.

One item the consolidated report lists as unverified has since been confirmed directly: arXiv:2509.10601 is by
El-Badry, Hogg and Rix, and its abstract states the result quoted (all 116 testable Gaia periodic candidates failed
on later data).

The facts that still rest on search summaries or recall, and that the design depends on, are listed in Appendix F of
the plan and at the end of the consolidated report.
