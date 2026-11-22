"""The quasar periodicity search, end to end and blind.

Order of events, which is the point of the design:

1. **Register.** The full configuration is hashed and written to the trials ledger before any light curve is
   analysed.
2. **Prepare.** Fixed quality rules: rescale error bars from same-night pairs, drop wild points, require enough
   nights and seasons.
3. **Twin.** Every object gets a decoy twin (a simulation of its own fitted red noise at its own sampling), and
   a fraction get a plant (the real light curve plus a known periodic signal).
4. **Seal.** Real objects, decoys and plants are shuffled under random labels. A hash of the label table is
   recorded.
5. **Analyse.** The same code searches every sealed subject and runs the vetting gates on the hits. It cannot
   tell which subjects are real.
6. **Unseal and score.** Probabilities are recalibrated against the population (an empirical null), the
   pre-registered promotion rule is applied, and three numbers come out: how many real objects were promoted,
   how many decoys were (the false-alarm estimate), and what fraction of plants were (the sensitivity, and
   from it upper limits).
7. **Record.** Whatever happened goes into the ledger, including "nothing found". Promoted real objects are put
   to the project-wide online test and get a dossier and a dated prediction.
"""
from __future__ import annotations

import itertools
import pathlib
import time as clock
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from ..core import Series, content_hash, rng_for
from ..detect import periodicity
from ..detect.periodicity import SearchConfig, band_scale, merge_bands
from ..ledger.predictions import PredictionLedger
from ..ledger.trials import TrialsLedger
from ..nulls import drw
from ..report import dossier
from ..trust import blind, empirical_null, inject
from ..trust.blind import Subject
from ..trust.decoys import uniformity
from ..trust.inject import RecoveryCell
from ..vet import gates as vet_gates
from ..vet import quality


@dataclass(frozen=True)
class RunConfig:
    """Everything that defines a run. Its hash is registered before the data are analysed."""

    name: str = "quasar-periodicity"
    seed: int = 1
    search: SearchConfig = field(default_factory=lambda: SearchConfig(max_decoys=1000))   # see SearchConfig
    min_points_band: int = 20        # a band with fewer points is dropped
    min_nights: int = 60             # an object needs this many nightly points over all bands
    min_span_days: float = 1000.0    # and this much time between its first and last night
    clip_sigma: float = 6.0
    pair_gap_days: float = 0.5       # measurements closer than this calibrate the error bars
    gate_below_p: float = 0.05       # gates are run on every subject below this raw probability
    hit_p: float = 0.01
    promote_p: float = 1e-3          # promotion needs this (after recalibration) and no failed gate
    deep_decoys: int = 20000         # promoted real objects are re-calibrated with this many simulated nulls
    plant_fraction: float = 0.4
    plant_amplitudes_mag: tuple[float, ...] = (0.05, 0.10, 0.20)
    plant_periods_days: tuple[float, ...] = (150.0, 300.0, 600.0)
    plant_shapes: tuple[str, ...] = ("sinusoid", "sawtooth", "flare")
    period_tolerance: float = 0.05   # a plant counts as recovered only if its period is found this closely
    workers: int = 4

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------------------------------------------------------
# The sample: Gaia quasar candidates in one sky cell, with their ZTF light curves
# ------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class SampleConfig:
    """Which objects are searched. Fixed and registered before any light curve is analysed."""

    healpix_order: int = 2
    healpix_pixel: int = 100
    g_max: float = 19.0                  # Gaia G magnitude limit
    min_quasar_prob: float = 0.9         # Gaia DR3 combined classifier probability of being a quasar
    max_parallax_over_error: float = 3.0   # a significant parallax means a star
    max_proper_motion: float = 3.0       # mas per year; quasars do not move
    match_radius_arcsec: float = 1.5
    max_objects: int | None = None       # None: all that qualify

    def to_dict(self) -> dict:
        return asdict(self)


def fetch_sample(sample: SampleConfig, *, cache=None, progress: bool = False
                 ) -> tuple[dict[str, list[Series]], dict[str, dict], dict]:
    """Gaia DR3 quasar candidates in the sample cell and their ZTF light curves.

    Returns (objects, info, report): light curves by name ("gaia:<source_id>"), what is known about each
    object, and a count of what was selected and matched. Objects without a ZTF light curve are counted and
    left out. Light curves are fetched one ZTF partition at a time to bound memory.
    """
    from ..sources import gaia, healpix_pixel, ztf

    candidates = gaia.quasar_candidates_in_pixel(sample.healpix_order, sample.healpix_pixel, g_max=sample.g_max,
                                                 min_quasar_prob=sample.min_quasar_prob, cache=cache)
    kept = []
    for c in candidates:
        parallax, motion = c["parallax_over_error"], c["pm"]
        if parallax is not None and parallax > sample.max_parallax_over_error:
            continue
        if motion is not None and motion > sample.max_proper_motion:
            continue
        kept.append(c)
    kept.sort(key=lambda c: c["source_id"])
    if sample.max_objects is not None:
        kept = kept[:sample.max_objects]
    cell = healpix_pixel(4, np.array([c["ra"] for c in kept]), np.array([c["dec"] for c in kept]))
    objects: dict[str, list[Series]] = {}
    info: dict[str, dict] = {}
    started = clock.time()
    for count, pixel in enumerate(np.unique(cell), 1):
        group = [c for c, p in zip(kept, cell) if p == pixel]
        curves = ztf.fetch_many([(f"gaia:{c['source_id']}", c["ra"], c["dec"]) for c in group],
                                radius_arcsec=sample.match_radius_arcsec, cache=cache, fragments="all")
        for c in group:
            name = f"gaia:{c['source_id']}"
            bands = curves.get(name) or []
            if not bands:
                continue
            objects[name] = bands
            info[name] = {"Gaia DR3 source_id": c["source_id"], "RA (deg)": round(c["ra"], 6),
                          "Dec (deg)": round(c["dec"], 6), "Gaia G (mag)": round(c["g_mag"], 2),
                          "Gaia quasar probability": round(c["quasar_prob"], 3),
                          "ZTF light curves (band: points)": ", ".join(f"{b.band}: {len(b)}" for b in bands),
                          "ZTF object ids": ", ".join(str(b.meta.get("objectid")) for b in bands)}
        if progress:
            print(f"  fetched cell {count} of {np.unique(cell).size}: {len(objects)} objects with light curves "
                  f"({clock.time() - started:.0f} s)", flush=True)
    report = {"gaia_candidates": len(candidates), "after_astrometric_cuts": len(kept),
              "with_ztf_light_curve": len(objects)}
    return objects, info, report


def register(config: "RunConfig", ledger_dir: str | pathlib.Path, *, sample: SampleConfig | None = None,
             description: str = "", code_version: str = "") -> dict:
    """Write the run's full definition to the trials ledger. Do this, and publish the ledger, before `run`."""
    trials = TrialsLedger(pathlib.Path(ledger_dir) / "trials.jsonl")
    definition = {"run": config.to_dict(), "sample": sample.to_dict() if sample else None,
                  "code_version": code_version}
    return trials.preregister(config.name, definition, description=description)


# ------------------------------------------------------------------------------------------------------------
# Prepare
# ------------------------------------------------------------------------------------------------------------

def prepare(objects: dict[str, list[Series]], config: RunConfig) -> tuple[list[Subject], dict]:
    """Apply the fixed quality rules. Returns the subjects that qualify and a report of what was done."""
    pooled: dict[str, list[np.ndarray]] = {}
    for bands in objects.values():
        for band in bands:
            pooled.setdefault(band.band, []).append(quality.pair_residuals(band, config.pair_gap_days))
    pooled_scale = {name: quality.error_scale(np.concatenate(parts), min_pairs=200) for name, parts in pooled.items()}
    subjects, dropped = [], {"too few nights": 0, "too short a time span": 0, "no usable band": 0}
    scales_used: dict[str, list[float]] = {}
    for name, bands in objects.items():
        kept = []
        for band in bands:
            own = quality.error_scale(quality.pair_residuals(band, config.pair_gap_days))
            scale = own if own is not None else pooled_scale.get(band.band)
            band = quality.rescale_errors(band, scale if scale is not None else 1.0)
            band = quality.clip_outliers(band, config.clip_sigma)
            scales_used.setdefault(band.band, []).append(band.meta["error_scale"])
            if len(band) >= config.min_points_band:
                kept.append(band)
        if not kept:
            dropped["no usable band"] += 1
            continue
        nightly = merge_bands(kept).binned(config.search.bin_days or 1.0)
        if len(nightly) < config.min_nights:
            dropped["too few nights"] += 1
        elif nightly.span < config.min_span_days:
            dropped["too short a time span"] += 1
        else:
            subjects.append(Subject(name, kept))
    report = {"objects_in": len(objects), "objects_kept": len(subjects), "dropped": dropped,
              "pooled_error_scale": pooled_scale,
              "median_error_scale_applied": {k: float(np.median(v)) for k, v in scales_used.items()}}
    return subjects, report


# ------------------------------------------------------------------------------------------------------------
# Twins
# ------------------------------------------------------------------------------------------------------------

def decoy_twin(subject: Subject, rng: np.random.Generator, config: RunConfig) -> Subject:
    """A simulation of this object's fitted red noise: one realisation, seen by every band at that band's own
    times with that band's own error bars, in that band's own units."""
    bands = subject.bands
    merged = merge_bands(bands)
    nightly = merged.binned(config.search.bin_days) if config.search.bin_days > 0 else merged
    sigma_grid, tau_grid = drw.default_grids(nightly.time, nightly.value, nightly.error,
                                             n_sigma=config.search.n_sigma, n_tau=config.search.n_tau)
    sigma, tau = drw.fit(nightly, sigma_grid=sigma_grid, tau_grid=tau_grid).draw(1, rng)
    scales = [band_scale(band) for band in bands]
    time = np.concatenate([band.time for band in bands])
    error = np.concatenate([band.error / scale for band, (_, scale) in zip(bands, scales)])
    order = np.argsort(time, kind="stable")
    simulated = np.empty(time.size)
    simulated[order] = drw.simulate(time[order], error[order], sigma, tau, rng)[0]
    out, start = [], 0
    for band, (centre, scale) in zip(bands, scales):
        part = simulated[start:start + len(band)]
        start += len(band)
        out.append(Series(band.time, centre + scale * part, band.error, band.band, band.source_id, band.kind,
                          band.time_standard))
    return Subject(subject.name, out, "decoy", {"sigma": float(sigma[0]), "tau_days": float(tau[0])})


def _flare(series: Series, period_days: float, amplitude: float, phase_cycles: float) -> Series:
    return inject.periodic_flares(series, period_days, amplitude, 0.15, phase_cycles)


_SHAPES = {"sinusoid": inject.sinusoid, "sawtooth": inject.sawtooth, "flare": _flare}


def plant_twin(subject: Subject, parameters: dict) -> Subject:
    """The real object with a known periodic signal of the same size (in magnitudes) added to every band."""
    shape = _SHAPES[parameters["shape"]]
    bands = [shape(band, parameters["period_days"], parameters["amplitude_mag"], parameters["phase_cycles"])
             for band in subject.bands]
    return Subject(subject.name, bands, "plant", dict(parameters))


def make_twins(subjects: list[Subject], config: RunConfig) -> list[Subject]:
    """One decoy for every subject and a plant for a fixed fraction, cycling evenly through the declared grid
    of planted signals."""
    out = []
    grid = list(itertools.product(config.plant_amplitudes_mag, config.plant_periods_days, config.plant_shapes))
    chooser = rng_for(config.seed, "plants")
    planted = chooser.permutation(len(subjects))[:int(round(config.plant_fraction * len(subjects)))]
    cells = chooser.permutation(len(grid))
    for subject in subjects:
        out.append(decoy_twin(subject, rng_for(config.seed, "decoy", subject.name), config))
    for count, index in enumerate(planted):
        amplitude, period, shape = grid[cells[count % len(grid)]]
        parameters = {"amplitude_mag": amplitude, "period_days": period, "shape": shape,
                      "phase_cycles": float(chooser.random())}
        out.append(plant_twin(subjects[index], parameters))
    return out


# ------------------------------------------------------------------------------------------------------------
# Analyse (blind) and score
# ------------------------------------------------------------------------------------------------------------

def analyse(subject: Subject, config: RunConfig) -> dict:
    """Search one subject and, if it is a hit, run the gates. Knows nothing of what the subject is."""
    merged = merge_bands(subject.bands)
    result = periodicity.search(merged, rng_for(config.seed, "analyse", subject.name), config.search)
    out = {"label": subject.name, "search": result.to_dict(), "gates": []}
    if result.status == "searched" and result.p <= config.gate_below_p:
        out["gates"] = [g.to_dict() for g in vet_gates.run_gates(merged, subject.bands, result, config.search)]
    return out


def _analyse_packed(args: tuple[Subject, RunConfig]) -> dict:
    return analyse(*args)


def probability_used(search: dict) -> float:
    """The simulated probability; when no simulated null came close, the fitted-tail extrapolation, capped at
    the simulation's own resolution so that an extrapolation can only refine, never invent, significance."""
    decoy = search.get("decoy")
    if not decoy:
        return float("nan")
    tail = search.get("p_tail", float("nan"))
    if decoy["exceedances"] == 0 and np.isfinite(tail):
        return float(min(decoy["p"], max(tail, 1e-300)))
    return float(decoy["p"])


def _write_results(rows: list[dict], out_dir: pathlib.Path) -> None:
    """Every subject's full record, compressed, and a small table of the columns that matter."""
    import csv
    import gzip
    import json

    with gzip.open(out_dir / "results.json.gz", "wt") as handle:
        json.dump(dossier._jsonable(rows), handle)
    gate_names = ("calendar", "window", "seasons", "bands", "holdout", "coherence", "alias")
    with open(out_dir / "results.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "kind", "name", "status", "nights", "period_days", "statistic", "p_raw",
                         "p_final", "coherence_p", "harmonic_fraction", *gate_names, "promoted", "planted"])
        for r in sorted(rows, key=lambda r: (r["kind"], r["p_final"] if np.isfinite(r["p_final"]) else 2.0)):
            s = r["search"]
            verdicts = {g["gate"]: g["passed"] for g in r["gates"]}
            planted = (f"{r['parameters']['shape']} {r['parameters']['amplitude_mag']} mag "
                       f"{r['parameters']['period_days']} d") if r["kind"] == "plant" else ""
            writer.writerow([r["label"], r["kind"], r["name"], s["status"], s["n_points"],
                             f"{s['period_days']:.2f}", f"{s['statistic']:.2f}", f"{r['p_raw']:.3g}",
                             f"{r['p_final']:.3g}", f"{s['coherence_p']:.3g}", f"{s['harmonic_fraction']:.3f}",
                             *["" if verdicts.get(g) is None else ("pass" if verdicts[g] else "fail")
                               for g in gate_names], "yes" if r["promoted"] else "no", planted])


def score(analysed: list[dict], key: dict[str, dict], config: RunConfig) -> dict:
    """Unsealed scoring: recalibrate, apply the promotion rule, count by kind, build the completeness table."""
    rows = []
    for item in analysed:
        truth = key[item["label"]]
        rows.append({**item, **truth, "searched": item["search"]["status"] == "searched",
                     "p_raw": probability_used(item["search"])})
    searched = [r for r in rows if r["searched"]]
    nulls = {}
    for kind in ("real", "decoy"):
        p = np.array([r["p_raw"] for r in searched if r["kind"] == kind])
        nulls[kind] = empirical_null.fit(p)
    for r in searched:
        null = nulls["decoy"] if r["kind"] == "decoy" else nulls["real"]
        r["p_final"] = float(null.recalibrate(np.array([r["p_raw"]]))[0])
        passed, failed = vet_gates.verdict([vet_gates.GateResult(**g) for g in r["gates"]])
        r["failed_gates"] = failed
        r["hit"] = r["p_final"] <= config.hit_p
        r["promoted"] = bool(r["p_final"] <= config.promote_p and r["gates"] and passed)
    for r in rows:
        r.setdefault("p_final", float("nan"))
        r.setdefault("hit", False)
        r.setdefault("promoted", False)
        r.setdefault("failed_gates", [])

    def count(kind: str, flag: str) -> int:
        return sum(1 for r in searched if r["kind"] == kind and r[flag])

    n = {kind: sum(1 for r in searched if r["kind"] == kind) for kind in ("real", "decoy", "plant")}
    summary: dict[str, Any] = {
        "searched": n,
        "not_searched": sum(1 for r in rows if not r["searched"]),
        "hits": {kind: count(kind, "hit") for kind in n},
        "promoted": {kind: count(kind, "promoted") for kind in n},
        "empirical_null": {kind: null.to_dict() for kind, null in nulls.items()},
        "raw_p_uniformity": {kind: uniformity(np.array([r["p_raw"] for r in searched if r["kind"] == kind]))
                             for kind in ("real", "decoy") if n[kind]},
        "decoy_estimate": blind.decoy_estimate(count("real", "promoted"), count("decoy", "promoted"), n["real"],
                                               n["decoy"]),
    }
    gated = [r for r in searched if r["kind"] == "decoy" and r["gates"]]
    summary["decoy_gate_pass_rate"] = {
        gate: float(np.mean([g["passed"] is not False for r in gated for g in r["gates"] if g["gate"] == gate]))
        for gate in ("calendar", "window", "seasons", "bands", "holdout", "coherence", "alias")} if gated else {}
    # Completeness: a plant is recovered if it is promoted with the right period.
    cells = []
    for amplitude, period, shape in itertools.product(config.plant_amplitudes_mag, config.plant_periods_days,
                                                      config.plant_shapes):
        mine = [r for r in searched if r["kind"] == "plant" and r["parameters"]["amplitude_mag"] == amplitude
                and r["parameters"]["period_days"] == period and r["parameters"]["shape"] == shape]
        right = [abs(r["search"]["period_days"] / period - 1) <= config.period_tolerance for r in mine]
        recovered = sum(1 for r, ok in zip(mine, right) if ok and r["promoted"])
        cell = RecoveryCell(shape, {"amplitude_mag": amplitude, "period_days": period}, len(mine), recovered).to_dict()
        cell["period_found"] = sum(right)
        cell["significant_before_gates"] = sum(1 for r, ok in zip(mine, right) if ok and r["p_final"] <= config.promote_p)
        if n["real"] and recovered:
            # With no real object promoted, the 95% upper limit on the fraction of objects carrying this signal.
            cell["fraction_upper95_if_none_found"] = float(min(1.0, 1 - 0.05 ** (1 / (n["real"] * cell["completeness"]))))
        cells.append(cell)
    summary["completeness"] = cells
    # Do real objects prefer particular periods more than their decoys do? A pile-up is instrumental.
    from scipy.stats import ks_2samp

    real_periods = np.array([r["search"]["period_days"] for r in searched if r["kind"] == "real"])
    decoy_periods = np.array([r["search"]["period_days"] for r in searched if r["kind"] == "decoy"])
    if real_periods.size > 20 and decoy_periods.size > 20:
        test = ks_2samp(np.log(real_periods), np.log(decoy_periods))
        summary["period_distribution_real_vs_decoy"] = {"ks_statistic": float(test.statistic),
                                                        "ks_p": float(test.pvalue)}
    return {"rows": rows, "summary": summary}


# ------------------------------------------------------------------------------------------------------------
# The whole run
# ------------------------------------------------------------------------------------------------------------

def run(objects: dict[str, list[Series]], config: RunConfig, out_dir: str | pathlib.Path, *,
        info: dict[str, dict] | None = None, ledger_dir: str | pathlib.Path | None = None,
        description: str = "", progress: bool = False, search_id: str | None = None) -> dict:
    """Run the blind search on `objects` (name -> bands) and write everything to `out_dir`.

    Returns the summary. With `ledger_dir`, the run is recorded to the trials ledger there, promoted real
    objects are put to the project-wide online test, and each gets a dated prediction. If the run was
    registered beforehand with `register`, pass the entry's hash as `search_id`: the configuration is then
    checked against what was registered, and the run refuses to start if they differ.
    """
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    started = clock.time()
    config_hash = content_hash(config.to_dict())
    trials = TrialsLedger(pathlib.Path(ledger_dir) / "trials.jsonl") if ledger_dir else None
    registration = None
    if trials and search_id:
        registration = next((e for e in trials.entries() if e["hash"] == search_id and e["type"] == "preregister"),
                            None)
        if registration is None:
            raise KeyError("no registered search with that id in the ledger")
        if registration["config"].get("run") != dossier._jsonable(config.to_dict()):
            raise ValueError("the configuration differs from the one registered; register the new one first")
    elif trials:
        registration = trials.preregister(config.name, {"run": config.to_dict(), "sample": None,
                                                        "code_version": ""}, description=description,
                                          planned_hypotheses=len(objects))
    dossier.dump_json({"config": config.to_dict(), "config_hash": config_hash}, out_dir / "config.json")

    subjects, preparation = prepare(objects, config)
    everything = subjects + make_twins(subjects, config)
    sealed = blind.seal(everything, rng_for(config.seed, "seal"))
    (out_dir / "sealed_commitment.txt").write_text(sealed.commitment + "\n")
    if trials:
        trials.note("sealed set committed", search_id=registration["hash"], commitment=sealed.commitment,
                    n_subjects=len(sealed.subjects))

    jobs = [(subject, config) for subject in sealed.subjects]
    if config.workers > 1 and len(jobs) > 8:
        with ProcessPoolExecutor(max_workers=config.workers) as pool:
            analysed = []
            for done, item in enumerate(pool.map(_analyse_packed, jobs, chunksize=4), 1):
                analysed.append(item)
                if progress and done % 200 == 0:
                    print(f"  analysed {done} of {len(jobs)} sealed subjects ({clock.time() - started:.0f} s)",
                          flush=True)
    else:
        analysed = [_analyse_packed(job) for job in jobs]

    # Keep the sealed results before anything else happens, so that a failure after this point costs nothing.
    import gzip
    import json

    with gzip.open(out_dir / "sealed_results.json.gz", "wt") as handle:
        json.dump(dossier._jsonable(analysed), handle)

    # The analysis is final. Only now is the label table opened.
    scored = score(analysed, sealed.unseal(), config)
    summary = scored["summary"]
    if config.deep_decoys:
        # A fixed rule that can only demote: every promoted real object is calibrated again, with a fresh and
        # much larger set of simulated nulls, and must still reach the promotion threshold.
        import dataclasses

        null = empirical_null.EmpiricalNull(**summary["empirical_null"]["real"])
        by_name = {subject.name: subject for subject in subjects}
        deep_config = dataclasses.replace(config.search, max_decoys=config.deep_decoys, coherence_below_p=-1.0)
        for r in scored["rows"]:
            if r["kind"] == "real" and r["promoted"]:
                deep = periodicity.search(merge_bands(by_name[r["name"]].bands),
                                          rng_for(config.seed, "deep", r["name"]), deep_config).to_dict()
                p_deep = float(null.recalibrate(np.array([probability_used(deep)]))[0])
                r["deep"] = {"p_raw": probability_used(deep), "p_final": p_deep, "draws": deep["decoy"]["draws"],
                             "exceedances": deep["decoy"]["exceedances"]}
                if p_deep > config.promote_p:
                    r["promoted"] = False
                    r["failed_gates"] = r["failed_gates"] + ["deep calibration"]
        promoted_real = sum(1 for r in scored["rows"] if r["kind"] == "real" and r["promoted"])
        summary["promoted_before_deep_calibration"] = dict(summary["promoted"])
        summary["promoted"]["real"] = promoted_real
        summary["decoy_estimate"] = blind.decoy_estimate(promoted_real, summary["promoted"]["decoy"],
                                                         summary["searched"]["real"], summary["searched"]["decoy"])
    summary.update({"run": config.name, "config_hash": config_hash, "preparation": preparation,
                    "sealed_commitment": sealed.commitment, "seconds": clock.time() - started})
    n_frequencies = sum(r["search"]["n_frequencies"] for r in scored["rows"] if r["kind"] == "real")
    summary["hypotheses_examined"] = int(n_frequencies)
    dossier.dump_json(sealed.proof(), out_dir / "sealed_proof.json")
    _write_results(scored["rows"], out_dir)

    context = {"run": config.name, "configuration hash": config_hash,
               "real objects searched": summary["searched"]["real"],
               "decoy twins promoted": f"{summary['promoted']['decoy']} of {summary['searched']['decoy']}"}
    real = sorted((r for r in scored["rows"] if r["kind"] == "real" and r["searched"]), key=lambda r: r["p_final"])
    written = []
    for r in real:
        if not (r["promoted"] or len(written) < 5):
            continue
        record = dossier.periodicity_dossier(r["name"], r["search"], r["gates"], promoted=r["promoted"],
                                             p_final=r["p_final"], info=(info or {}).get(r["name"]),
                                             context=context)
        if r["promoted"] and trials:
            formal = trials.test(f"{config.name}: {r['name']} is periodic at {r['search']['period_days']:.1f} d",
                                 r["p_final"], search_id=registration["hash"])
            record["project_wide_test"] = {k: formal[k] for k in ("p", "level", "discovery")}
            ledger = PredictionLedger(pathlib.Path(ledger_dir) / "predictions.jsonl")
            statement = (f"{r['name']}: a sinusoid of period {r['search']['period_days']:.2f} d and phase "
                         f"{r['search']['phase_cycles']:.3f} cycles at time zero, fitted to data up to the end of "
                         "this run's light curve, will give a positive log-likelihood gain over red noise alone "
                         "on data taken after it.")
            record["prediction"] = ledger.commit_open(statement, due="next ZTF data release",
                                                      metadata={"run": config.name, "object": r["name"]})
        written.append(str(dossier.write(record, out_dir / "dossiers")))
    summary["dossiers"] = written
    summary["promoted_real"] = [{"name": r["name"], "period_days": r["search"]["period_days"],
                                 "p_final": r["p_final"]} for r in real if r["promoted"]]
    summary["top_real"] = [{"name": r["name"], "period_days": r["search"]["period_days"], "p_raw": r["p_raw"],
                            "p_final": r["p_final"], "failed_gates": r["failed_gates"],
                            "promoted": r["promoted"]} for r in real[:10]]
    dossier.dump_json(summary, out_dir / "summary.json")
    if trials:
        trials.record(registration["hash"], n_objects=summary["searched"]["real"], n_hypotheses=n_frequencies,
                      summary={k: summary[k] for k in ("searched", "hits", "promoted", "decoy_estimate",
                                                       "empirical_null", "sealed_commitment")})
    return summary
