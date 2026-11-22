"""ZTF DR24 light curves, from the public HATS copy on S3 (route "cloud") or the IRSA light-curve service.

Both routes return the same structure: one `Series` per band (g, r, i), magnitudes in the ZTF native system,
time in heliocentric modified Julian date (HJD - 2400000.5, the catalogue's `hmjd`).
"""
from __future__ import annotations

import csv
import io
import logging
import urllib.parse
from dataclasses import dataclass
from typing import Any

import numpy as np
import pyarrow as pa

from ..core import Provenance, Series, content_hash
from .cache import Cache, fetch_bytes, utc_now
from .hats import SEPARATION_COLUMN, HatsCatalog, angular_separation_arcsec
from .registry import ZTF_CITATION

log = logging.getLogger(__name__)

SOURCE_ID = "ztf_dr24_lc"
HATS_URL = "https://ipac-irsa-ztf.s3.amazonaws.com/ztf/enhanced/dr24/lc/hats/ztf_dr24_lc-hats"
IRSA_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
TIME_STANDARD = "HJD - 2400000.5 (ZTF hmjd)"
BANDS = {1: "g", 2: "r", 3: "i"}
FILTERCODES = {"zg": 1, "zr": 2, "zi": 3}
COLUMNS = ["objectid", "filterid", "nepochs", "lightcurve.hmjd", "lightcurve.mag", "lightcurve.magerr",
           "lightcurve.catflags"]
IRSA_TIMEOUT = 300.0   # the service took 150 s for a single 2-arcsecond cone on 22 Nov 2026


@dataclass
class _Curve:
    """One object's light curve in one band, as read from either route."""

    objectid: int
    filterid: int
    hmjd: np.ndarray
    mag: np.ndarray
    magerr: np.ndarray
    catflags: np.ndarray
    separation_arcsec: float
    nepochs: int


def catalog(cache: Cache | None = None) -> HatsCatalog:
    """The ZTF DR24 light-curve HATS catalogue on S3 (position columns objra, objdec)."""
    return HatsCatalog(HATS_URL, "objra", "objdec", cache=cache)


def fetch_lightcurves(ra: float, dec: float, radius_arcsec: float = 2.0, *, route: str = "cloud",
                      keep_flagged: bool = False, cache: Cache | None = None,
                      fragments: str = "longest") -> list[Series]:
    """ZTF DR24 light curves of the object(s) within radius_arcsec of (ra, dec) (ICRS degrees).

    Returns one Series per band found (band "g", "r" or "i"), ordered g, r, i. If several objects in a band lie
    within the radius, the closest is used and the others' ids go in meta["other_objectids"]. value and error are
    ZTF PSF magnitudes (kind "mag"); time is hmjd (TIME_STANDARD). By default only epochs with catflags == 0 are
    kept; keep_flagged=True keeps all and stores the flags in meta["catflags"]. Epochs with a non-finite value or
    a non-positive error are always dropped (counted in meta["n_invalid_dropped"]). A band left with no epochs
    is omitted. route is "cloud" (HATS Parquet on S3) or "irsa" (the IRSA light-curve service).
    """
    query = {"ra": ra, "dec": dec, "radius_arcsec": radius_arcsec, "route": route, "keep_flagged": keep_flagged}
    if route == "cloud":
        table = catalog(cache).read_cone(ra, dec, radius_arcsec, COLUMNS)
        curves = curves_from_table(table)
        query["catalog"] = HATS_URL
    elif route == "irsa":
        url = irsa_url(ra, dec, radius_arcsec)
        data, _ = fetch_bytes(url, cache=cache, timeout=IRSA_TIMEOUT)
        curves = parse_irsa_csv(data.decode("utf-8"), ra, dec)
        query["url"] = url
    else:
        raise ValueError(f"route must be 'cloud' or 'irsa', not {route!r}")
    return build_series(curves, query, keep_flagged, fragments)


def fetch_many(positions: list[tuple[str, float, float]], radius_arcsec: float = 2.0, *,
               keep_flagged: bool = False, cache: Cache | None = None,
               fragments: str = "longest") -> dict[str, list[Series]]:
    """Light curves for many (name, ra, dec) targets from the cloud copy, reading each row group at most once.

    Returns {name: list of Series} with the same per-target structure as fetch_lightcurves (an empty list when
    nothing lies within the radius). Read costs (row groups, bytes, seconds) are logged at INFO by the
    stargazer.sources.hats logger.
    """
    reader = catalog(cache)
    tables = reader.read_positions(positions, radius_arcsec, COLUMNS)
    out = {}
    for name, ra, dec in positions:
        query = {"name": name, "ra": ra, "dec": dec, "radius_arcsec": radius_arcsec, "route": "cloud",
                 "keep_flagged": keep_flagged, "catalog": HATS_URL}
        out[name] = build_series(curves_from_table(tables[name]), query, keep_flagged, fragments)
    return out


def curves_from_table(table: pa.Table) -> list[_Curve]:
    """Per-object curves from a cone read of the HATS catalogue (columns as in COLUMNS plus the separation)."""
    curves = []
    for row in table.to_pylist():
        lc = row["lightcurve"]
        curves.append(_Curve(int(row["objectid"]), int(row["filterid"]), np.asarray(lc["hmjd"], dtype=float),
                             np.asarray(lc["mag"], dtype=float), np.asarray(lc["magerr"], dtype=float),
                             np.asarray(lc["catflags"], dtype=np.int64), float(row[SEPARATION_COLUMN]),
                             int(row["nepochs"])))
    return curves


def irsa_url(ra: float, dec: float, radius_arcsec: float) -> str:
    """The IRSA ZTF light-curve service request (CSV) for a cone; radius is sent in degrees."""
    params = {"POS": f"CIRCLE {ra!r} {dec!r} {radius_arcsec / 3600.0:.6g}", "BANDNAME": "g,r,i", "FORMAT": "csv"}
    return f"{IRSA_URL}?{urllib.parse.urlencode(params, quote_via=urllib.parse.quote)}"


def parse_irsa_csv(text: str, ra: float, dec: float) -> list[_Curve]:
    """Per-object curves from the IRSA light-curve service's CSV output.

    Rows are grouped by oid. time is hjd - 2400000.5, matching the cloud copy's hmjd. The separation is that of
    the median per-epoch position from (ra, dec).
    """
    reader = csv.DictReader(io.StringIO(text))
    needed = {"oid", "hjd", "mag", "magerr", "catflags", "filtercode", "ra", "dec"}
    if reader.fieldnames is None or not needed <= set(reader.fieldnames):
        raise ValueError(f"IRSA light-curve CSV lacks columns {needed - set(reader.fieldnames or [])}")
    groups: dict[tuple[int, int], list[dict[str, str]]] = {}
    for row in reader:
        code = row["filtercode"].strip()
        if code not in FILTERCODES:
            raise ValueError(f"unknown ZTF filtercode {code!r}")
        groups.setdefault((int(row["oid"]), FILTERCODES[code]), []).append(row)
    curves = []
    for (oid, filterid), rows in groups.items():
        col = {k: np.array([float(r[k]) for r in rows]) for k in ("hjd", "mag", "magerr", "ra", "dec")}
        sep = float(angular_separation_arcsec(ra, dec, np.median(col["ra"]), np.median(col["dec"])))
        curves.append(_Curve(oid, filterid, col["hjd"] - 2400000.5, col["mag"], col["magerr"],
                             np.array([int(r["catflags"]) for r in rows], dtype=np.int64), sep, len(rows)))
    return curves


def _usable(curve: _Curve, keep_flagged: bool) -> np.ndarray:
    valid = np.isfinite(curve.hmjd) & np.isfinite(curve.mag) & np.isfinite(curve.magerr) & (curve.magerr > 0)
    return valid if keep_flagged else valid & (curve.catflags == 0)


def build_series(curves: list[_Curve], query: dict[str, Any], keep_flagged: bool,
                 fragments: str = "longest") -> list[Series]:
    """Series from per-object curves, flags handled as documented in fetch_lightcurves, with a filled Provenance.

    ZTF stores one light curve per (field, chip quadrant, filter), so a source that falls in two overlapping
    fields has two or more fragments in the same band, often of very different length (a primary field with
    hundreds of epochs and a secondary one with a few dozen). `fragments` says what to return for each band:
    "longest" (default) the fragment with the most usable epochs, "closest" the one nearest the requested
    position, or "all" every fragment as its own Series, longest first. Fragments have independent zero points
    and must not be concatenated without removing an offset; `detect.periodicity.merge_bands` does that.
    """
    if fragments not in ("longest", "closest", "all"):
        raise ValueError('fragments must be "longest", "closest" or "all"')
    out = []
    retrieved = utc_now()
    for filterid in sorted(BANDS):
        in_band = [c for c in curves if c.filterid == filterid]
        if not in_band:
            continue
        if fragments == "closest":
            in_band.sort(key=lambda c: c.separation_arcsec)
        else:
            in_band.sort(key=lambda c: (-int(np.sum(_usable(c, keep_flagged))), c.separation_arcsec))
        for curve in (in_band if fragments == "all" else in_band[:1]):
            valid = np.isfinite(curve.hmjd) & np.isfinite(curve.mag) & np.isfinite(curve.magerr) & (curve.magerr > 0)
            good = curve.catflags == 0
            keep = valid if keep_flagged else valid & good
            meta: dict[str, Any] = {
                "objectid": curve.objectid, "nepochs": curve.nepochs,
                "n_flagged_dropped": 0 if keep_flagged else int(np.sum(valid & ~good)),
                "n_invalid_dropped": int(np.sum(~valid)), "separation_arcsec": curve.separation_arcsec,
            }
            others = [c.objectid for c in in_band if c is not curve]
            if others:
                meta["other_objectids"] = others
            if keep_flagged:
                order = np.argsort(curve.hmjd[keep], kind="stable")
                meta["catflags"] = curve.catflags[keep][order].tolist()
            if not np.any(keep):
                log.info("ZTF object %d band %s has no usable epochs; omitted", curve.objectid, BANDS[filterid])
                continue
            t, v, e = curve.hmjd[keep], curve.mag[keep], curve.magerr[keep]
            series = Series(t, v, e, band=BANDS[filterid], source_id=f"ztf:{curve.objectid}", kind="mag",
                            time_standard=TIME_STANDARD, meta=meta)
            series.provenance = Provenance(SOURCE_ID, dict(query, band=BANDS[filterid]), retrieved,
                                           content_hash({"time": series.time, "value": series.value,
                                                         "error": series.error}), ZTF_CITATION)
            out.append(series)
    return out
