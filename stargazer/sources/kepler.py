"""Kepler long-cadence light curves from the public MAST copy on S3 (s3://stpubdata/kepler/public/lightcurves)."""
from __future__ import annotations

import concurrent.futures as futures
import io
import logging
import urllib.parse
from typing import Any

import numpy as np
from astropy.io import fits

from ..core import Provenance, Series, content_hash
from .cache import Cache, fetch_bytes, utc_now
from .hats import parse_s3_listing
from .registry import KEPLER_CITATION

log = logging.getLogger(__name__)

SOURCE_ID = "kepler_lc"
BUCKET_URL = "https://stpubdata.s3.amazonaws.com"
PREFIX = "kepler/public/lightcurves"
TIME_STANDARD = "BJD_TDB - 2400000.5"


def kic_prefix(kic: int) -> str:
    """The S3 key prefix holding a target's files: <first four digits>/<nine-digit zero-padded id>/."""
    if kic <= 0 or kic >= 10 ** 9:
        raise ValueError(f"KIC id out of range: {kic}")
    padded = f"{int(kic):09d}"
    return f"{PREFIX}/{padded[:4]}/{padded}/"


def list_long_cadence(kic: int, cache: Cache | None = None) -> list[str]:
    """HTTPS URLs of the long-cadence (_llc.fits) files of one KIC target, in time order."""
    prefix = kic_prefix(kic)
    keys: list[str] = []
    token = None
    while True:
        query = {"list-type": "2", "prefix": prefix}
        if token:
            query["continuation-token"] = token
        data, _ = fetch_bytes(f"{BUCKET_URL}/?{urllib.parse.urlencode(query)}", cache=cache)
        page, token = parse_s3_listing(data)
        keys += page
        if not token:
            break
    files = sorted(k for k in keys if k.endswith("_llc.fits"))
    if not files:
        raise LookupError(f"no Kepler long-cadence files for KIC {kic} under {BUCKET_URL}/{prefix}")
    return [f"{BUCKET_URL}/{k}" for k in files]


def read_quarter(data: bytes, flux_column: str = "PDCSAP_FLUX") -> dict[str, Any]:
    """One quarter from a Kepler light-curve FITS file's bytes.

    Returns a dict with quarter, time (BJD_TDB - 2400000.5), value and error (flux divided by the quarter's
    median flux, so dimensionless), median_flux (e-/s), n_raw, n_dropped. Points with SAP_QUALITY != 0 or a
    non-finite time, flux or error, or a non-positive error, are dropped.
    """
    with fits.open(io.BytesIO(data), memmap=False) as hdul:
        quarter = int(hdul[0].header["QUARTER"])
        table = hdul["LIGHTCURVE"].data
        header = hdul["LIGHTCURVE"].header
        offset = float(header["BJDREFI"]) + float(header["BJDREFF"]) - 2400000.5
        time = np.asarray(table["TIME"], dtype=float) + offset
        flux = np.asarray(table[flux_column], dtype=float)
        error = np.asarray(table[f"{flux_column}_ERR"], dtype=float)
        quality = np.asarray(table["SAP_QUALITY"])
    keep = (quality == 0) & np.isfinite(time) & np.isfinite(flux) & np.isfinite(error) & (error > 0)
    if not np.any(keep):
        median = float("nan")
    else:
        median = float(np.median(flux[keep]))
        if not median > 0:
            raise ValueError(f"quarter {quarter}: median {flux_column} is {median}, cannot normalise")
    return {"quarter": quarter, "time": time[keep], "value": flux[keep] / median, "error": error[keep] / median,
            "median_flux": median, "n_raw": int(time.size), "n_dropped": int(np.sum(~keep))}


def fetch_lightcurve(kic: int, *, flux_column: str = "PDCSAP_FLUX", cache: Cache | None = None,
                     workers: int = 4) -> Series:
    """All long-cadence quarters of a Kepler target, each divided by its own median flux and joined.

    Returns a Series of relative flux (kind "flux", dimensionless, about 1), time in BJD_TDB - 2400000.5
    (Kepler BKJD + 54832.5), error from the matching "<flux_column>_ERR" column scaled the same way. meta holds
    kic, flux_column, quarters, points_dropped (total and per quarter) and normalisation (median e-/s per
    quarter). flux_column may be "PDCSAP_FLUX" (default) or "SAP_FLUX".
    """
    urls = list_long_cadence(kic, cache)
    with futures.ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        blobs = list(pool.map(lambda url: fetch_bytes(url, cache=cache)[0], urls))
    quarters = [read_quarter(blob, flux_column) for blob in blobs]
    quarters = [q for q in quarters if q["time"].size]
    if not quarters:
        raise LookupError(f"KIC {kic}: no usable points in {len(urls)} long-cadence files")
    meta: dict[str, Any] = {
        "kic": int(kic), "flux_column": flux_column, "quarters": [q["quarter"] for q in quarters],
        "points_dropped": int(sum(q["n_dropped"] for q in quarters)),
        "points_dropped_by_quarter": {str(q["quarter"]): q["n_dropped"] for q in quarters},
        "normalisation": {"method": "divide each quarter by its median flux",
                          "median_flux_e_per_s": {str(q["quarter"]): q["median_flux"] for q in quarters}},
        "files": len(urls),
    }
    series = Series(np.concatenate([q["time"] for q in quarters]), np.concatenate([q["value"] for q in quarters]),
                    np.concatenate([q["error"] for q in quarters]), band="Kepler", source_id=f"kic:{int(kic)}",
                    kind="flux", time_standard=TIME_STANDARD, meta=meta)
    series.provenance = Provenance(
        SOURCE_ID, {"kic": int(kic), "flux_column": flux_column, "prefix": f"{BUCKET_URL}/{kic_prefix(kic)}"},
        utc_now(), content_hash({"time": series.time, "value": series.value, "error": series.error}),
        KEPLER_CITATION)
    log.info("Kepler KIC %d: %d quarters, %d points (%d dropped)", kic, len(quarters), len(series),
             meta["points_dropped"])
    return series
