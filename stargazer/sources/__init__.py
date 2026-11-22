"""Data access: connectors to public archives, a HATS catalogue reader, a download cache and a source registry.

Every connector returns `stargazer.core.Series` objects carrying a filled `Provenance` (source id, query,
retrieval time, content hash of the returned arrays, citation).
"""

from . import fermi, frb, gaia, kepler, ztf
from .cache import Cache, FetchError, fetch_bytes
from .hats import HatsCatalog, ReadStats, cone_cells, healpix_pixel, prune_row_groups
from .registry import SOURCES, Reachability, SourceInfo, check

__all__ = [
    "Cache", "FetchError", "HatsCatalog", "ReadStats", "Reachability", "SOURCES", "SourceInfo", "check",
    "cone_cells", "fermi", "fetch_bytes", "frb", "gaia", "healpix_pixel", "kepler", "prune_row_groups", "ztf",
]
