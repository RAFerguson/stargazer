"""Gaia DR3 quasar candidates from the HATS copy of gaia_source on S3 (s3://stpubdata/gaia/gaia_dr3/public/hats).

That copy is partitioned by HEALPix at orders 2 to 6, one ".parquet" file per partition, sorted within each
file by source_id (and so by `_healpix_29`: Gaia source ids encode an order-12 pixel), which lets row-group
statistics prune by sky position. Column names checked against dataset/_common_metadata on 22 Nov 2026:
in_qso_candidates (bool), classprob_dsc_combmod_quasar (float), phot_g_mean_mag (float), parallax_over_error,
pm (total proper motion, mas/yr), phot_variable_flag (string), ra, dec (degrees, epoch 2016.0).
"""
from __future__ import annotations

import logging
from typing import Any

import pyarrow as pa
import pyarrow.compute as pc

from .cache import Cache
from .hats import HatsCatalog

log = logging.getLogger(__name__)

SOURCE_ID = "gaia_dr3_hats"
HATS_URL = "https://stpubdata.s3.amazonaws.com/gaia/gaia_dr3/public/hats/gaia"
COLUMNS = ["source_id", "ra", "dec", "phot_g_mean_mag", "in_qso_candidates", "classprob_dsc_combmod_quasar",
           "parallax_over_error", "pm", "phot_variable_flag"]


def catalog(cache: Cache | None = None) -> HatsCatalog:
    """The Gaia DR3 gaia_source HATS catalogue on S3 (position columns ra, dec)."""
    return HatsCatalog(HATS_URL, "ra", "dec", cache=cache)


def _selection(g_max: float, min_quasar_prob: float) -> pc.Expression:
    return (pc.field("in_qso_candidates")
            & (pc.field("classprob_dsc_combmod_quasar") >= min_quasar_prob)
            & (pc.field("phot_g_mean_mag") <= g_max))


def _bounds(g_max: float, min_quasar_prob: float) -> dict[str, tuple[Any, Any]]:
    return {"in_qso_candidates": (True, True), "classprob_dsc_combmod_quasar": (min_quasar_prob, None),
            "phot_g_mean_mag": (None, g_max)}


def _records(table: pa.Table) -> list[dict[str, Any]]:
    """Selected rows as dicts; null measurements become None."""
    out = []
    for row in table.to_pylist():
        out.append({"source_id": int(row["source_id"]), "ra": row["ra"], "dec": row["dec"],
                    "g_mag": row["phot_g_mean_mag"], "quasar_prob": row["classprob_dsc_combmod_quasar"],
                    "parallax_over_error": row["parallax_over_error"], "pm": row["pm"],
                    "phot_variable_flag": row["phot_variable_flag"]})
    return out


def quasar_candidates(ra: float, dec: float, radius_deg: float, *, g_max: float = 19.0,
                      min_quasar_prob: float = 0.9, cache: Cache | None = None) -> list[dict[str, Any]]:
    """Gaia DR3 quasar candidates within radius_deg of (ra, dec) (ICRS degrees), nearest first.

    Selects in_qso_candidates, classprob_dsc_combmod_quasar >= min_quasar_prob and phot_g_mean_mag <= g_max.
    Each dict has source_id, ra, dec (degrees), g_mag (Gaia G, mag), quasar_prob (DSC combined probability),
    parallax_over_error, pm (total proper motion, mas/yr), phot_variable_flag ("VARIABLE" or "NOT_AVAILABLE"),
    and separation_arcsec.
    """
    table = catalog(cache).read_cone(ra, dec, radius_deg * 3600.0, COLUMNS)
    table = table.filter(_selection(g_max, min_quasar_prob))
    table = table.sort_by("_separation_arcsec")
    records = _records(table)
    for record, sep in zip(records, table.column("_separation_arcsec").to_pylist()):
        record["separation_arcsec"] = sep
    return records


def quasar_candidates_in_pixel(order: int, pixel: int, *, g_max: float = 19.0, min_quasar_prob: float = 0.9,
                               cache: Cache | None = None) -> list[dict[str, Any]]:
    """Gaia DR3 quasar candidates in one nested HEALPix cell (order, pixel), with the same selection and dict
    fields as quasar_candidates (without separation_arcsec), ordered by source_id.

    Only the nine needed columns are read. Row groups are pruned by sky position and by the statistics of the
    selection columns; every remaining row group of the covering partition(s) is downloaded. Measured on
    22 Nov 2026: order-4 cell 1608 (inside order-2 partition 100) 10.8 MB, 1-2.4 s; the whole order-2 partition
    100 (1.19 million rows, 1,426 candidates) 53 MB, 3.2 s; the order-4 partition 112 (3 row groups) 24 MB, 2.1 s.
    """
    table = catalog(cache).read_pixel(order, pixel, COLUMNS, row_filter=_selection(g_max, min_quasar_prob),
                                      bounds=_bounds(g_max, min_quasar_prob))
    return _records(table.sort_by("source_id"))
