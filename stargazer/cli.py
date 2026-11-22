"""Command line: `stargazer <command>`.

    stargazer sources check                 which archives answer from here
    stargazer fetch ztf RA DEC -o lc.json   one object's ZTF light curves
    stargazer search lc.json                periodicity search and vetting gates on saved light curves
    stargazer quasar register|run           the blind quasar search (register first, publish, then run)
    stargazer decode FILE                   content ladder for a bit stream
    stargazer bench                         the benchmark scorecard
    stargazer ledger                        what has been tested so far, and whether the record is intact
"""
from __future__ import annotations

import json
import pathlib
import subprocess

import click
import numpy as np

from .core import BitStream, Series, rng_for


def _code_version() -> str:
    try:
        root = pathlib.Path(__file__).resolve().parent.parent
        return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return ""


@click.group()
def main() -> None:
    """Search public sky data for signals that natural models do not explain, and say honestly what was found."""


# ---- sources --------------------------------------------------------------------------------------------------

@main.group()
def sources() -> None:
    """Data sources."""


@sources.command("check")
@click.option("--timeout", default=20.0, show_default=True)
def sources_check(timeout: float) -> None:
    """Probe every registered archive and say which are reachable."""
    from .sources import registry

    for item in registry.check(timeout=timeout):
        click.echo(f"{item.status:12s} {item.id:22s} {item.host:34s} {item.detail}")


@main.group()
def fetch() -> None:
    """Fetch data for one object."""


@fetch.command("ztf")
@click.argument("ra", type=float)
@click.argument("dec", type=float)
@click.option("--radius", default=2.0, show_default=True, help="match radius in arcseconds")
@click.option("-o", "--out", type=click.Path(), required=True)
def fetch_ztf(ra: float, dec: float, radius: float, out: str) -> None:
    """ZTF light curves at RA DEC (degrees), one series per band, saved as JSON."""
    from .sources import ztf

    bands = ztf.fetch_lightcurves(ra, dec, radius)
    pathlib.Path(out).write_text(json.dumps([band.to_dict() for band in bands]))
    click.echo(f"{len(bands)} bands: " + ", ".join(f"{b.band} {len(b)} points" for b in bands) + f" -> {out}")


# ---- search ---------------------------------------------------------------------------------------------------

@main.command()
@click.argument("light_curves", type=click.Path(exists=True))
@click.option("--seed", default=1, show_default=True)
@click.option("--max-decoys", default=2000, show_default=True)
def search(light_curves: str, seed: int, max_decoys: int) -> None:
    """Periodicity search and vetting gates on light curves saved by `fetch`."""
    from .detect import periodicity
    from .vet import gates

    bands = [Series.from_dict(d) for d in json.loads(pathlib.Path(light_curves).read_text())]
    config = periodicity.SearchConfig(max_decoys=max_decoys)
    merged = periodicity.merge_bands(bands)
    result = periodicity.search(merged, rng_for(seed, "cli-search", light_curves), config)
    click.echo(f"status: {result.status}; {result.n_points} nightly points over {result.span_days:.0f} days")
    if result.status != "searched":
        return
    click.echo(f"best period {result.period_days:.1f} d ({result.cycles_covered:.1f} cycles), statistic "
               f"{result.statistic:.1f}, probability against red noise {result.p:.2g} "
               f"({result.decoy.exceedances} of {result.decoy.draws} decoys)")
    click.echo(f"white-noise test would have said {result.classical['white_noise_fap']:.1e} (not evidence)")
    for flag in result.flags:
        click.echo(f"note: {flag}")
    for gate in gates.run_gates(merged, bands, result, config):
        verdict = "n/a " if gate.passed is None else ("pass" if gate.passed else "FAIL")
        click.echo(f"  [{verdict}] {gate.gate}: {gate.sentence}")


# ---- the blind quasar search ----------------------------------------------------------------------------------

@main.group()
def quasar() -> None:
    """The blind quasar periodicity search."""


def _configs(name: str, seed: int, order: int, pixel: int, max_objects: int | None, workers: int):
    from .pipelines import quasar as pipeline

    return (pipeline.RunConfig(name=name, seed=seed, workers=workers),
            pipeline.SampleConfig(healpix_order=order, healpix_pixel=pixel, max_objects=max_objects))


_COMMON = [click.option("--name", required=True, help="a name for the run"),
           click.option("--seed", default=1, show_default=True),
           click.option("--order", default=2, show_default=True, help="HEALPix order of the sky cell"),
           click.option("--pixel", default=100, show_default=True, help="nested HEALPix pixel of the sky cell"),
           click.option("--max-objects", default=None, type=int),
           click.option("--workers", default=4, show_default=True),
           click.option("--ledger", "ledger_dir", default="ledger", show_default=True, type=click.Path())]


def _common(function):
    for option in reversed(_COMMON):
        function = option(function)
    return function


@quasar.command("register")
@_common
@click.option("--description", default="")
def quasar_register(name, seed, order, pixel, max_objects, workers, ledger_dir, description) -> None:
    """Write the run's definition to the trials ledger. Commit and publish the ledger before `run`."""
    from .pipelines import quasar as pipeline

    config, sample = _configs(name, seed, order, pixel, max_objects, workers)
    entry = pipeline.register(config, ledger_dir, sample=sample, description=description,
                              code_version=_code_version())
    click.echo(f"registered {name}: search id {entry['hash']}")
    click.echo(f"configuration hash {entry['config_hash']}; code version {entry['config']['code_version'] or 'unknown'}")


@quasar.command("run")
@_common
@click.option("--search-id", required=True, help="the id printed by `register`")
@click.option("--out", "out_dir", required=True, type=click.Path())
@click.option("--cache", "cache_dir", default=None, type=click.Path())
def quasar_run(name, seed, order, pixel, max_objects, workers, ledger_dir, search_id, out_dir, cache_dir) -> None:
    """Fetch the registered sample and run the blind search on it."""
    from .pipelines import quasar as pipeline
    from .sources import Cache

    config, sample = _configs(name, seed, order, pixel, max_objects, workers)
    objects, info, report = pipeline.fetch_sample(sample, cache=Cache(cache_dir) if cache_dir else None,
                                                  progress=True)
    click.echo(f"sample: {report}")
    summary = pipeline.run(objects, config, out_dir, info=info, ledger_dir=ledger_dir, progress=True,
                           search_id=search_id)
    summary["sample"] = report
    from .report import dossier

    dossier.dump_json(summary, pathlib.Path(out_dir) / "summary.json")
    click.echo(json.dumps({k: summary[k] for k in ("searched", "hits", "promoted", "decoy_estimate",
                                                   "empirical_null")}, indent=1, default=float))


# ---- decode ---------------------------------------------------------------------------------------------------

@main.command()
@click.argument("path", type=click.Path(exists=True))
@click.option("--bytes", "as_bytes", is_flag=True, help="read the file as raw bytes (most significant bit first)")
@click.option("--seed", default=0, show_default=True)
@click.option("-o", "--out", type=click.Path(), default=None, help="write the reading as Markdown here")
def decode(path: str, as_bytes: bool, seed: int, out: str | None) -> None:
    """Climb the content ladder for a bit stream: a text file of 0s and 1s ('?' for unknown), or raw bytes."""
    from .decode import workbench

    if as_bytes:
        stream = BitStream(np.unpackbits(np.frombuffer(pathlib.Path(path).read_bytes(), dtype=np.uint8)))
    else:
        stream = BitStream.from_text("".join(pathlib.Path(path).read_text().split()))
    reading = workbench.analyse(stream, name=pathlib.Path(path).name, seed=seed)
    text = workbench.to_markdown(reading)
    if out:
        pathlib.Path(out).write_text(text)
        click.echo(f"content tier {reading.tier} of 5: {reading.tier_phrase} -> {out}")
    else:
        click.echo(text)


# ---- bench ----------------------------------------------------------------------------------------------------

@main.command()
@click.option("--seed", default=0, show_default=True)
@click.option("--size", default="small", type=click.Choice(["small", "full"]), show_default=True)
@click.option("--workers", default=4, show_default=True)
@click.option("-o", "--out", type=click.Path(), default=None, help="directory for scorecard.md and scorecard.json")
def bench(seed: int, size: str, workers: int, out: str | None) -> None:
    """Run the decoding workbench on every benchmark signal and compare with the known answers."""
    from .bench import score
    from .report import dossier

    card = score.run(seed=seed, size=size, workers=workers)
    text = score.to_markdown(card)
    click.echo(text)
    if out:
        directory = pathlib.Path(out)
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "scorecard.md").write_text(text)
        dossier.dump_json(card, directory / "scorecard.json")


# ---- ledger ---------------------------------------------------------------------------------------------------

@main.command()
@click.option("--dir", "ledger_dir", default="ledger", show_default=True, type=click.Path())
def ledger(ledger_dir: str) -> None:
    """What has been tested so far, how many discoveries, and whether the record is intact."""
    from .ledger.predictions import PredictionLedger
    from .ledger.trials import TrialsLedger

    directory = pathlib.Path(ledger_dir)
    click.echo(json.dumps(TrialsLedger(directory / "trials.jsonl").summary(), indent=1))
    predictions = PredictionLedger(directory / "predictions.jsonl")
    click.echo(json.dumps({"predictions": predictions.hit_rate(), "chain_intact": predictions.verify()}, indent=1))


if __name__ == "__main__":
    main()
