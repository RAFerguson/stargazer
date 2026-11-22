"""Tests for stargazer.sources.

Unmarked tests run offline against small recorded excerpts in tests/fixtures/sources and files the tests write.
Tests marked `network` read the live archives (deselected by default; run with `pytest -m network`).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from stargazer.core import Series
from stargazer.sources import SOURCES, Cache, HatsCatalog, check, fermi, gaia, kepler, ztf
from stargazer.sources.cache import cache_key
from stargazer.sources.hats import (angular_separation_arcsec, cell_range, cone_cells, healpix_pixel,
                                    parse_partition_info, parse_s3_listing, partitions_for_cells, prune_row_groups)

FIXTURES = Path(__file__).parent / "fixtures" / "sources"
RA, DEC = 196.38754, -10.555389     # PG 1302-102


# ------------------------------------------------------------------------------------------------ offline


def test_healpix_pixels_match_catalogue_numbering():
    assert int(healpix_pixel(4, RA, DEC)[0]) == 1608
    assert int(healpix_pixel(9, RA, DEC)[0]) == 1647540
    lo, hi = cell_range(4, 1608, 9)
    assert lo <= 1647540 <= hi and hi - lo + 1 == 4 ** 5
    assert cell_range(9, 1647540, 4) == (1608, 1608)


def test_cone_cells_small_and_edge():
    order, cells = cone_cells(RA, DEC, 2.0)
    assert order == 9 and 1647540 in cells
    # A cone straddling a pixel boundary must include both pixels.
    hp_size = 3600 * math.degrees(math.sqrt(4 * math.pi / (12 * 4 ** 9)))
    ra2 = RA
    for step in np.linspace(0, hp_size / 3600 * 2, 400):   # walk east until the order-9 pixel changes
        if int(healpix_pixel(9, RA + step, DEC)[0]) != 1647540:
            ra2 = RA + step
            break
    edge_ra = ra2 - 0.5 / 3600
    order, cells = cone_cells(edge_ra, DEC, 2.0)
    assert {int(healpix_pixel(9, edge_ra, DEC)[0]), int(healpix_pixel(9, ra2 + 0.5 / 3600, DEC)[0])} <= cells
    # Large cones drop to a coarser order.
    order, cells = cone_cells(RA, DEC, 3600.0)
    assert order < 9 and int(healpix_pixel(order, RA, DEC)[0]) in cells


def test_angular_separation():
    assert angular_separation_arcsec(10.0, 0.0, np.array([10.0]), np.array([1.0]))[0] == pytest.approx(3600.0)
    assert angular_separation_arcsec(0.0, 89.0, np.array([180.0]), np.array([89.0]))[0] == pytest.approx(7200.0)


def test_partition_lookup_from_recorded_excerpt():
    parts = parse_partition_info((FIXTURES / "ztf_dr24_partition_info_excerpt.csv").read_text())
    assert (4, 1608, 770558) in parts
    index: dict[int, np.ndarray] = {}
    for o, p, _ in parts:
        index.setdefault(o, []).append(p)  # type: ignore[arg-type]
    index = {o: np.array(sorted(p)) for o, p in index.items()}
    order, cells = cone_cells(RA, DEC, 2.0)
    assert partitions_for_cells(index, order, cells) == [(4, 1608)]
    # A coarse cell picks up every finer partition inside it.
    assert partitions_for_cells(index, 3, [402]) == [(4, 1608), (4, 1609), (4, 1610), (4, 1611)]
    # Without a num_rows column (the Gaia layout) rows are -1.
    assert parse_partition_info("Norder,Npix\n2,0\n2,1\n") == [(2, 0, -1), (2, 1, -1)]


def test_catalog_metadata_parsing(monkeypatch):
    texts = {"hats.properties": "#HATS\nhats_col_ra=ra\nhats_npix_suffix=/\nhats_creation_date=2025-06-13T19\\:56UTC\n",
             "partition_info.csv": (FIXTURES / "ztf_dr24_partition_info_excerpt.csv").read_text()}
    cat = HatsCatalog("https://example.invalid/cat", "ra", "dec")
    monkeypatch.setattr(cat, "_get_text", lambda rel: texts[rel])
    props = cat.properties()
    assert props["pixel_is_directory"] is True and props["hats_creation_date"] == "2025-06-13T19:56UTC"
    assert cat.partitions_for_cone(RA, DEC, 2.0) == [(4, 1608)]


def test_s3_listing_parser():
    xml = (b'<?xml version="1.0" encoding="UTF-8"?><ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">'
           b"<IsTruncated>true</IsTruncated><NextContinuationToken>abc</NextContinuationToken>"
           b"<Contents><Key>a/part0.snappy.parquet</Key></Contents><Contents><Key>a/x.txt</Key></Contents>"
           b"</ListBucketResult>")
    assert parse_s3_listing(xml) == (["a/part0.snappy.parquet", "a/x.txt"], "abc")


def _write_sorted_parquet(path: Path) -> None:
    """Rows sorted by _healpix_9 across 8 row groups, so that pruning has something to prune."""
    pix9 = np.sort(np.random.default_rng(1).integers(1646592, 1647616, 8000)).astype(np.int32)
    table = pa.table({"_healpix_9": pix9, "_healpix_29": pix9.astype(np.int64) << 40,
                      "ra": np.linspace(0, 1, pix9.size), "dec": np.zeros(pix9.size), "value": np.arange(pix9.size)})
    pq.write_table(table, path, row_group_size=1000)


def test_row_group_pruning(tmp_path):
    path = tmp_path / "sorted.parquet"
    _write_sorted_parquet(path)
    metadata = pq.read_metadata(path)
    assert metadata.num_row_groups == 8
    keep = prune_row_groups(metadata, [(9, 1647540)])
    assert 1 <= len(keep) <= 2
    table = pq.ParquetFile(path).read_row_groups(keep)
    all_rows = pq.read_table(path)
    assert pa.compute.sum(pa.compute.equal(table["_healpix_9"], 1647540)).as_py() == \
        pa.compute.sum(pa.compute.equal(all_rows["_healpix_9"], 1647540)).as_py() > 0
    # A whole order-4 cell keeps everything; a cell outside the file keeps nothing.
    assert prune_row_groups(metadata, [(4, 1608)]) == list(range(8))
    assert prune_row_groups(metadata, [(4, 1700)]) == []
    # Extra bounds prune on other columns' statistics.
    assert prune_row_groups(metadata, [(4, 1608)], {"value": (None, 1500)}) == [0, 1]


def test_row_group_pruning_without_stats(tmp_path):
    path = tmp_path / "nostats.parquet"
    table = pa.table({"_healpix_9": np.arange(4000, dtype=np.int32), "ra": np.zeros(4000), "dec": np.zeros(4000)})
    pq.write_table(table, path, row_group_size=1000, write_statistics=False)
    assert prune_row_groups(pq.read_metadata(path), [(9, 5)]) == [0, 1, 2, 3]


def test_cache_round_trip(tmp_path):
    cache = Cache(tmp_path / "c")
    key = cache_key("https://example.invalid/file", (0, 99))
    assert not cache.has(key) and cache.get(key) is None
    sha = cache.put(key, b"hello")
    assert sha == "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
    assert cache.has(key) and cache.get(key) == b"hello"
    info = cache.info(key)
    assert info["key"] == key and info["size"] == 5 and info["sha256"] == sha and info["stored_utc"].endswith("Z")
    assert Cache(tmp_path / "c").get(key) == b"hello"


def test_cache_default_root(monkeypatch, tmp_path):
    monkeypatch.setenv("STARGAZER_CACHE", str(tmp_path / "env"))
    assert Cache().root == tmp_path / "env"
    monkeypatch.delenv("STARGAZER_CACHE")
    assert Cache().root == Path.home() / ".stargazer" / "cache"


def test_fetch_bytes_uses_cache(tmp_path):
    from stargazer.sources import fetch_bytes
    cache = Cache(tmp_path)
    url = "https://example.invalid/never-fetched"
    cache.put(cache_key(url), b"cached")
    assert fetch_bytes(url, cache=cache)[0] == b"cached"


def test_registry():
    hosts = set()
    for line in (Path(__file__).parents[1] / "spikes" / "first-light" / "sources.txt").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            hosts.add(line.split("|")[0].strip())
    assert hosts <= {s.host for s in SOURCES.values()}
    assert "10.26131/IRSA598" in SOURCES["ztf_dr24_lc"].citation
    assert "Gaia" in SOURCES["gaia_dr3_hats"].citation and "DPAC" in SOURCES["gaia_dr3_hats"].citation
    assert "MAST" in SOURCES["kepler_lc"].citation
    assert "Light Curve Repository" in SOURCES["fermi_lcr"].citation
    assert all(s.route in ("cloud", "api", "code") for s in SOURCES.values())
    with pytest.raises(KeyError):
        check(["no_such_source"])


def test_irsa_csv_parser():
    text = (FIXTURES / "irsa_ztf_lc_196.38754_-10.555389_excerpt.csv").read_text()
    curves = ztf.parse_irsa_csv(text, RA, DEC)
    assert {(c.objectid, c.filterid) for c in curves} == {(372106200003324, 1), (372206200012063, 2)}
    for c in curves:
        assert c.separation_arcsec < 1.0 and 58000 < c.hmjd.min() < 61000
    query = {"ra": RA, "dec": DEC, "route": "irsa"}
    kept = ztf.build_series(curves, query, keep_flagged=False)
    every = ztf.build_series(curves, query, keep_flagged=True)
    assert [s.band for s in kept] == ["g", "r"]
    rows = {}
    for line in text.splitlines()[1:]:
        rows[int(line.split(",")[0])] = rows.get(int(line.split(",")[0]), 0) + 1
    for k, e in zip(kept, every):
        n = rows[e.meta["objectid"]]
        assert len(e) == e.meta["nepochs"] == n
        assert len(k) + k.meta["n_flagged_dropped"] == n and k.meta["n_flagged_dropped"] > 0
        assert sum(f != 0 for f in e.meta["catflags"]) == k.meta["n_flagged_dropped"]
        assert k.time_standard == ztf.TIME_STANDARD and k.kind == "mag" and k.source_id.startswith("ztf:")
        assert k.provenance.source == "ztf_dr24_lc" and len(k.provenance.content_sha256) == 64
    first = text.splitlines()[1].split(",")
    assert kept[0].time[0] == pytest.approx(float(first[2]) - 2400000.5)


def test_fermi_json_parser():
    payload = json.loads((FIXTURES / "fermi_lcr_4FGL_J1555.7+1111_weekly_excerpt.json").read_text())
    parsed = fermi.parse_lcr_json(payload)
    assert parsed["n_bins"] == 28 and parsed["n_upper_limits"] == 4 and parsed["n_not_converged"] == 0
    assert parsed["met"].size == 24 and np.all(parsed["error"] > 0)
    assert parsed["flux"][0] == pytest.approx(6.39e-08)
    assert parsed["error"][0] == pytest.approx((7.87e-08 - 4.91e-08) / 2)
    # A detection whose fit did not converge is dropped and counted.
    met0 = payload["flux"][0][0]
    payload["fit_convergence"] = [[m, 1 if m == met0 else c] for m, c in payload["fit_convergence"]]
    parsed = fermi.parse_lcr_json(payload)
    assert parsed["n_not_converged"] == 1 and parsed["met"].size == 23
    assert fermi.met_to_mjd(239889601) == pytest.approx(54686.5 + 7.428704e-4 + 1 / 86400, abs=1e-9)


# ------------------------------------------------------------------------------------------------ live


@pytest.mark.network
def test_live_ztf_cloud():
    curves = {s.band: s for s in ztf.fetch_lightcurves(RA, DEC, 2.0, keep_flagged=True)}
    assert len(curves["g"]) == 254 and len(curves["r"]) == 306
    assert curves["g"].source_id == "ztf:372106200003324" and curves["r"].source_id == "ztf:372206200012063"
    clean = {s.band: s for s in ztf.fetch_lightcurves(RA, DEC, 2.0)}
    for band in ("g", "r"):
        assert len(clean[band]) + clean[band].meta["n_flagged_dropped"] == len(curves[band])
        assert clean[band].meta["separation_arcsec"] < 1.0
        assert "doi:10.26131/IRSA598" in clean[band].provenance.citation


@pytest.mark.network
def test_live_ztf_irsa_matches_cloud():
    irsa = {s.band: s for s in ztf.fetch_lightcurves(RA, DEC, 2.0, route="irsa", keep_flagged=True)}
    assert len(irsa["g"]) == 254 and len(irsa["r"]) == 306
    cloud = {s.band: s for s in ztf.fetch_lightcurves(RA, DEC, 2.0, keep_flagged=True)}
    for band in ("g", "r"):
        assert irsa[band].source_id == cloud[band].source_id
        # The CSV's hjd carries fewer digits than the Parquet hmjd: differences up to 5e-6 d (0.4 s) were seen.
        assert np.allclose(irsa[band].time, cloud[band].time, atol=2e-5)
        assert np.allclose(irsa[band].value, cloud[band].value, atol=1e-5)


@pytest.mark.network
def test_live_ztf_fetch_many():
    positions = [("pg1302", RA, DEC), ("offset", RA + 0.01, DEC + 0.01), ("nothing", 0.0, -89.0)]
    out = ztf.fetch_many(positions, 2.0, keep_flagged=True)
    assert {s.band: len(s) for s in out["pg1302"]}["g"] == 254
    assert out["nothing"] == []
    assert all(isinstance(s, Series) for v in out.values() for s in v)


@pytest.mark.network
def test_live_gaia_finds_the_quasar():
    found = gaia.quasar_candidates(RA, DEC, 10 / 3600)
    assert [q["source_id"] for q in found] == [3622979843899961472]
    assert found[0]["quasar_prob"] >= 0.9 and found[0]["g_mag"] < 19 and found[0]["separation_arcsec"] < 1


@pytest.mark.network
def test_live_gaia_pixel():
    found = gaia.quasar_candidates_in_pixel(9, 1647540)
    assert 3622979843899961472 in {q["source_id"] for q in found}
    assert all(q["quasar_prob"] >= 0.9 and q["g_mag"] <= 19.0 for q in found)


@pytest.mark.network
def test_live_kepler():
    series = kepler.fetch_lightcurve(8462852)
    assert len(series) > 50000 and len(series.meta["quarters"]) == 18
    assert series.kind == "flux" and abs(np.median(series.value) - 1) < 0.01
    assert 54900 < series.time[0] < 55000


@pytest.mark.network
def test_live_fermi():
    series = fermi.fetch_lightcurve("4FGL J1555.7+1111", cadence="weekly")
    assert series.span > 15 * 365.25
    assert series.meta["n_bins"] == len(series) + series.meta["n_upper_limits_dropped"] + \
        series.meta["n_not_converged_dropped"]
    with pytest.raises(LookupError):
        fermi.fetch_lightcurve("PG 1553+113")
