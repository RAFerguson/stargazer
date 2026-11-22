"""The written report of a blind periodicity run, generated from its `summary.json` and nothing else.

Numbers in the report are copied from the summary, never retyped, so the report cannot drift from the results.
The wording follows the outcome: a run with nothing promoted says so in its first sentence and then says what
that does and does not rule out.
"""
from __future__ import annotations

import json
import pathlib


def _pct(value: float) -> str:
    return f"{100 * value:.0f}%"


def _probability(p: float) -> str:
    if p is None:
        return "n/a"
    if p >= 0.01:
        return f"{p:.2f}"
    return f"{p:.1e}" if p < 1e-4 else f"{p:.4f}"


def periodicity_report(run_dir: str | pathlib.Path, *, title: str, question: str, figures: dict | None = None,
                       notes: list[str] | None = None) -> str:
    """Markdown for one run. `figures` maps "calibration" and "completeness" to (light path, dark path) relative
    to the report; `notes` are extra paragraphs for the caveats section."""
    run_dir = pathlib.Path(run_dir)
    s = json.loads((run_dir / "summary.json").read_text())
    n, hits, promoted = s["searched"], s["hits"], s["promoted"]
    estimate = s["decoy_estimate"]
    sample, prep = s.get("sample", {}), s["preparation"]
    uniform = s["raw_p_uniformity"]
    null = s["empirical_null"]["real"]
    lines = [f"# {title}", "", f"**Question.** {question}", ""]

    # ---- the answer --------------------------------------------------------------------------------------
    if promoted["real"] == 0:
        answer = (f"**Answer.** No. Of {n['real']:,} quasars searched, none was promoted. {hits['real']} reached the "
                  f"hit level and all of them failed vetting or fell short of the promotion threshold.")
    else:
        answer = (f"**Answer.** {promoted['real']} of {n['real']:,} quasars searched were promoted to vetted hits. "
                  "A vetted hit is a periodic component that red noise as modelled, the calendar and the sampling "
                  "do not explain. It is not a confirmed signal.")
    answer += (f" Of the {n['decoy']:,} decoy twins, which contain nothing by construction, {promoted['decoy']} "
               f"were promoted, so about {estimate['expected_false_among_real']:.0f} false promotions are expected "
               "among the real objects.")
    planted = [c for c in s["completeness"] if c["trials"]]
    recovered = sum(c["recovered"] for c in planted)
    trials = sum(c["trials"] for c in planted)
    answer += (f" Of {trials} signals planted in real light curves, {recovered} ({_pct(recovered / trials)}) were "
               "promoted at the right period; the table below says which kinds.")
    lines += [answer, ""]

    # ---- the sample --------------------------------------------------------------------------------------
    lines += ["## The sample", ""]
    if sample:
        lines.append(f"- Gaia DR3 quasar candidates in the sky cell: {sample['gaia_candidates']:,}; after parallax "
                     f"and proper-motion cuts: {sample['after_astrometric_cuts']:,}; with a ZTF light curve: "
                     f"{sample['with_ztf_light_curve']:,}.")
    dropped = ", ".join(f"{count} for {reason}" for reason, count in prep["dropped"].items() if count)
    lines.append(f"- After the fixed quality rules: {prep['objects_kept']:,} kept"
                 + (f" ({dropped} dropped)." if dropped else "."))
    scales = ", ".join(f"{band} {value:.2f}" for band, value in prep["pooled_error_scale"].items() if value)
    if scales:
        lines.append(f"- Scatter of same-night pairs relative to the quoted error bars, pooled by band: {scales}. "
                     "Values below 1 mean the quoted errors are, if anything, too large; errors were never shrunk.")
    lines.append(f"- Trial periods examined across all real objects: {s['hypotheses_examined']:,}.")
    lines.append("")

    # ---- calibration -------------------------------------------------------------------------------------
    lines += ["## Did the search behave?", ""]
    if figures and "calibration" in figures:
        light, dark = figures["calibration"]
        lines += ["<picture>", f'  <source media="(prefers-color-scheme: dark)" srcset="{dark}">',
                  '  <img alt="Share of real objects and of decoy twins at or below each probability, on logarithmic '
                  f'axes, with the line a calibrated null follows." src="{light}">', "</picture>", ""]
    lines += ["| | Objects | At or below 0.05 | Expected | At or below 0.01 | Expected |", "|---|---|---|---|---|---|"]
    for kind, label in (("real", "Real objects"), ("decoy", "Decoy twins")):
        u = uniform[kind]
        lines.append(f"| {label} | {u['n']:,} | {u['below_0.05']} | {u['expected_below_0.05']:.0f} | "
                     f"{u['below_0.01']} | {u['expected_below_0.01']:.0f} |")
    lines += ["", f"The population of real objects sets its own null: fitted centre {null['centre']:+.2f} and width "
              f"{null['width']:.2f} on the z-score scale, where a perfectly modelled population gives 0 and 1. "
              + ("The fitted width is not above 1, so no recalibration was applied."
                 if null["width_used"] == 1.0 and null["centre_used"] == 0.0 else
                 f"Probabilities were recalibrated with centre {null['centre_used']:.2f} and width "
                 f"{null['width_used']:.2f}, which can only make them less significant."), ""]
    rates = s.get("decoy_gate_pass_rate") or {}
    if rates:
        lines.append("Share of decoy hits that each gate let through: "
                     + ", ".join(f"{gate} {_pct(rate)}" for gate, rate in rates.items()) + ".")
        lines.append("")
    pile = s.get("period_distribution_real_vs_decoy")
    if pile:
        lines.append("Do real objects prefer particular periods more than their decoys do? Kolmogorov–Smirnov "
                     f"probability {pile['ks_p']:.2f} (a small value would point to an instrumental period).")
        lines.append("")

    # ---- the most significant real objects -----------------------------------------------------------------
    lines += ["## The most significant real objects", "",
              "| Object | Best period (days) | Probability | After recalibration | Gates failed | Promoted |",
              "|---|---|---|---|---|---|"]
    for row in s["top_real"]:
        lines.append(f"| {row['name']} | {row['period_days']:.1f} | {_probability(row['p_raw'])} | "
                     f"{_probability(row['p_final'])} | {', '.join(row['failed_gates']) or 'none'} | "
                     f"{'yes' if row['promoted'] else 'no'} |")
    lines += ["", "Each of the five most significant, and every promoted object, has a dossier in `dossiers/`.", ""]

    # ---- completeness ------------------------------------------------------------------------------------
    lines += ["## What the search would have found", ""]
    if figures and "completeness" in figures:
        light, dark = figures["completeness"]
        lines += ["<picture>", f'  <source media="(prefers-color-scheme: dark)" srcset="{dark}">',
                  '  <img alt="Three grids, one per signal shape, giving the share of planted signals recovered for each '
                  f'amplitude and period." src="{light}">', "</picture>", ""]
    lines += ["| Shape | Amplitude (mag) | Period (days) | Planted | Right period | Promoted | 95% interval | "
              "Share of quasars that could carry it unseen |", "|---|---|---|---|---|---|---|---|"]
    for c in planted:
        # Conservative: the lower end of the recovery interval, not the point estimate, sets the limit.
        limit = min(1.0, 1 - 0.05 ** (1 / (n["real"] * c["low95"]))) if c["low95"] > 0 else None
        lines.append(f"| {c['label']} | {c['amplitude_mag']:.2f} | {c['period_days']:.0f} | {c['trials']} | "
                     f"{c['period_found']} | {c['recovered']} ({_pct(c['completeness'])}) | "
                     f"{_pct(c['low95'])}–{_pct(c['high95'])} | "
                     + (f"below {100 * limit:.1f}%" if limit is not None and promoted["real"] == 0 and limit < 1
                        else "not constrained")
                     + " |")
    lines += ["", "The last column is an upper limit, at 95% confidence, on the fraction of quasars in this sample that "
              "carry such a signal, given that none was promoted. It uses the lower end of the recovery interval, "
              "so it errs on the side of a weaker limit. Where the recovery interval reaches zero, the search says "
              "nothing about that kind of signal.", ""]

    # ---- caveats -----------------------------------------------------------------------------------------
    lines += ["## What this does and does not mean", "",
              "- The search covers strict periods from 40 days to a third of each light curve's length, away from "
              "the year and its fractions. Shorter periods, and periods near a year, are not tested.",
              "- The noise model is one damped random walk per object. The comparison of real objects with their "
              "decoy twins above is the check on that model.",
              "- A null result here is a statement about these objects, this period range and these amplitudes. "
              "It is not a statement that nothing periodic exists."]
    lines += [f"- {note}" for note in (notes or [])]
    lines += ["", "## Record", "",
              f"- Run `{s['run']}`, configuration hash `{s['config_hash']}`.",
              f"- Hash of the hidden label table, committed before the analysis: `{s['sealed_commitment']}`.",
              f"- Compute time: {s['seconds'] / 60:.0f} minutes.", ""]
    return "\n".join(lines)
