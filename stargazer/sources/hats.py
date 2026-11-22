"""Read HATS catalogues (Parquet partitioned by HEALPix pixel) straight from HTTPS object storage.

A HATS catalogue is a directory holding `hats.properties`, `partition_info.csv` and
`dataset/Norder=<o>/Dir=<d>/Npix=<p><suffix>`, where each partition covers one nested HEALPix pixel of
varying order and `<suffix>` is ".parquet" (one file) or "/" (a directory of part files). Rows carry
`_healpix_29` (and in some catalogues `_healpix_9`) pixel numbers; the Parquet row-group statistics on those
columns let a reader download only the row groups that can hold a given patch of sky.

This module never uses the `hats` or `lsdb` packages: pixel arithmetic is done with astropy_healpix, and the
files are opened with fsspec over HTTPS (HTTP range requests through the proxy in the environment).

Reads go in two stages. First the position columns of every candidate row group are read (cheap: a few bytes
per row) and the angular separation is computed; only row groups that actually hold a match are then read for
the requested columns, and only matching rows are kept. This matters because some catalogues (ZTF DR24) are not
spatially sorted inside a partition, so row-group statistics prune little and the heavy light-curve columns are
the bulk of every row group.
"""
from __future__ import annotations

import concurrent.futures as futures
import csv
import io
import logging
import math
import threading
import time
import urllib.parse
import xml.etree.ElementTree as ElementTree
from dataclasses import dataclass, field
from typing import Any, Iterable

import astropy.units as u
import fsspec
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from astropy_healpix import HEALPix

from .cache import Cache, FetchError, fetch_bytes
from .registry import USER_AGENT

log = logging.getLogger(__name__)

MAX_ORDER = 29
FINE_ORDER = 9            # the order of the `_healpix_9` column, and the finest order used for cone cells
SEPARATION_COLUMN = "_separation_arcsec"

Cell = tuple[int, int]    # (order, nested pixel)


# ----------------------------------------------------------------------------------------------- pixel geometry


def healpix_pixel(order: int, ra: float | np.ndarray, dec: float | np.ndarray) -> np.ndarray:
    """Nested HEALPix pixel number(s) at the given order for RA, Dec in degrees (ICRS), as an int64 array."""
    hp = HEALPix(nside=2 ** order, order="nested")
    pix = hp.lonlat_to_healpix(np.atleast_1d(ra) * u.deg, np.atleast_1d(dec) * u.deg)
    return np.asarray(pix, dtype=np.int64)


def angular_separation_arcsec(ra1: float, dec1: float, ra2: np.ndarray, dec2: np.ndarray) -> np.ndarray:
    """Great-circle separation in arcseconds between one position and arrays of positions (all in degrees)."""
    ra1r, dec1r = math.radians(ra1), math.radians(dec1)
    ra2r, dec2r = np.radians(np.asarray(ra2, dtype=float)), np.radians(np.asarray(dec2, dtype=float))
    hav = np.sin((dec2r - dec1r) / 2) ** 2 + math.cos(dec1r) * np.cos(dec2r) * np.sin((ra2r - ra1r) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0)))) * 3600.0


def _pixel_size_arcsec(order: int) -> float:
    return math.degrees(math.sqrt(4 * math.pi / (12 * 4 ** order))) * 3600.0


def cone_cells(ra: float, dec: float, radius_arcsec: float) -> tuple[int, set[int]]:
    """HEALPix cells (order, set of nested pixels) that together cover a cone; all angles in degrees/arcsec.

    The order is the finest (at most 9) whose pixels are at least twice the radius across. The cells are those
    touched by the centre and by four points one radius away (north, south, east, west), as in the first-light
    prototype, which is enough when the cone is small next to a pixel; astropy_healpix's conservative
    `cone_search_lonlat` is added to the set so that cones near pixel corners, and large cones, are covered too.
    """
    if radius_arcsec <= 0:
        raise ValueError("radius must be positive")
    order = FINE_ORDER
    while order > 0 and _pixel_size_arcsec(order) < 2 * radius_arcsec:
        order -= 1
    step = radius_arcsec / 3600.0
    stretch = max(math.cos(math.radians(dec)), 1e-6)
    points = [(ra, dec), (ra + step / stretch, dec), (ra - step / stretch, dec),
              (ra, min(dec + step, 90.0)), (ra, max(dec - step, -90.0))]
    rr, dd = zip(*points)
    cells = {int(p) for p in healpix_pixel(order, np.mod(rr, 360.0), dd)}
    hp = HEALPix(nside=2 ** order, order="nested")
    cells.update(int(p) for p in hp.cone_search_lonlat(ra * u.deg, dec * u.deg, radius=radius_arcsec * u.arcsec))
    return order, cells


def cell_range(order: int, pixel: int, at_order: int) -> tuple[int, int]:
    """Inclusive range of pixel numbers at `at_order` (>= order) that lie inside the cell (order, pixel)."""
    if at_order < order:
        shift = 2 * (order - at_order)
        return pixel >> shift, pixel >> shift
    shift = 2 * (at_order - order)
    return pixel << shift, ((pixel + 1) << shift) - 1


# ----------------------------------------------------------------------------------------------- row groups


def _leaf_index(metadata: pq.FileMetaData, name: str) -> int | None:
    for k in range(metadata.num_columns):
        if metadata.schema.column(k).path == name:
            return k
    return None


def _stats_range(metadata: pq.FileMetaData, group: int, leaf: int) -> tuple[Any, Any] | None:
    stats = metadata.row_group(group).column(leaf).statistics
    if stats is None or not stats.has_min_max:
        return None
    return stats.min, stats.max


def prune_row_groups(metadata: pq.FileMetaData, cells: Iterable[Cell],
                     bounds: dict[str, tuple[Any, Any]] | None = None) -> list[int]:
    """Indices of the row groups that can hold rows inside any of the HEALPix cells.

    A row group is dropped when the min/max statistics of `_healpix_9` or `_healpix_29` (whichever the file has)
    show that none of its rows can lie in any cell. A file with neither column, or without statistics, keeps
    every row group. `bounds` adds further {column: (low, high)} tests (None for an open end): a row group whose
    [min, max] for that column cannot meet [low, high] is dropped too.
    """
    cells = list(cells)
    tests: list[tuple[int, list[tuple[Any, Any]]]] = []
    for name, at_order in (("_healpix_9", FINE_ORDER), ("_healpix_29", MAX_ORDER)):
        leaf = _leaf_index(metadata, name)
        if leaf is not None:
            tests.append((leaf, [cell_range(o, p, at_order) for o, p in cells]))
    for name, (low, high) in (bounds or {}).items():
        leaf = _leaf_index(metadata, name)
        if leaf is not None:
            tests.append((leaf, [(low, high)]))
    keep = []
    for group in range(metadata.num_row_groups):
        ok = True
        for leaf, ranges in tests:
            stats = _stats_range(metadata, group, leaf)
            if stats is None:
                continue
            lo, hi = stats
            if not any((a is None or hi >= a) and (b is None or lo <= b) for a, b in ranges):
                ok = False
                break
        if ok:
            keep.append(group)
    return keep


# ----------------------------------------------------------------------------------------------- files


class _CountingFile(io.RawIOBase):
    """A read-only file wrapper that counts the bytes read through it (that is, downloaded)."""

    def __init__(self, inner: Any, counter: "ReadStats", url: str) -> None:
        self._inner = inner
        self._counter = counter
        self._url = url

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, offset: int, whence: int = 0) -> int:
        return self._inner.seek(offset, whence)

    def tell(self) -> int:
        return self._inner.tell()

    def read(self, size: int = -1) -> bytes:
        data = self._inner.read(size)
        self._counter.add_bytes(self._url, len(data))
        return data

    def readinto(self, buffer: Any) -> int:
        data = self.read(len(buffer))
        buffer[:len(data)] = data
        return len(data)

    def close(self) -> None:
        self._inner.close()
        super().close()


@dataclass
class ReadStats:
    """What a read cost: files opened, row groups available/scanned/read in full, bytes downloaded, seconds."""

    files: int = 0
    row_groups_total: int = 0
    row_groups_scanned: int = 0
    row_groups_read: int = 0
    bytes_read: int = 0
    seconds: float = 0.0
    cache_hits: int = 0
    by_file: dict[str, int] = field(default_factory=dict)
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False, compare=False)

    def add_bytes(self, url: str, n: int) -> None:
        with self._lock:
            self.bytes_read += n
            self.by_file[url] = self.by_file.get(url, 0) + n

    def add(self, **counts: int) -> None:
        with self._lock:
            for name, n in counts.items():
                setattr(self, name, getattr(self, name) + n)

    def summary(self) -> str:
        return (f"{self.files} file(s), {self.row_groups_scanned} of {self.row_groups_total} row groups scanned, "
                f"{self.row_groups_read} read in full, {self.bytes_read / 1e6:.1f} MB downloaded, "
                f"{self.cache_hits} cache hit(s), {self.seconds:.1f} s")


def _table_to_bytes(table: pa.Table) -> bytes:
    sink = pa.BufferOutputStream()
    pq.write_table(table, sink, compression="zstd")
    return sink.getvalue().to_pybytes()


def _table_from_bytes(data: bytes) -> pa.Table:
    return pq.read_table(pa.BufferReader(data))


def _unique(names: Iterable[str]) -> list[str]:
    out: list[str] = []
    for name in names:
        if name not in out:
            out.append(name)
    return out


class HatsCatalog:
    """One HATS catalogue served over HTTPS.

    base_url is the catalogue directory (the one holding hats.properties); ra_column and dec_column name the
    position columns (degrees). cache, if given, keeps metadata files, Parquet footers, position columns and
    matched rows, so repeating a query needs no network. workers is how many row groups of one file are
    downloaded at once (each holds one decompressed row group in memory; about 0.5 GB for ZTF light curves).
    """

    def __init__(self, base_url: str, ra_column: str, dec_column: str, cache: Cache | None = None,
                 timeout: float = 120.0, workers: int = 3) -> None:
        self.base_url = base_url.rstrip("/")
        self.ra_column = ra_column
        self.dec_column = dec_column
        self.cache = cache
        self.timeout = timeout
        self.workers = max(1, workers)
        self.last_stats = ReadStats()
        self._properties: dict[str, Any] | None = None
        self._partitions: list[tuple[int, int, int]] | None = None
        self._by_order: dict[int, np.ndarray] | None = None
        self._files: dict[Cell, list[str]] = {}
        self._sizes: dict[str, int] = {}
        self._empty: dict[tuple[str, ...], pa.Table] = {}
        self._fs = fsspec.filesystem("https", client_kwargs={"trust_env": True,
                                                             "headers": {"User-Agent": USER_AGENT}})

    # ---- metadata

    def _get_text(self, relative: str) -> str:
        data, _ = fetch_bytes(f"{self.base_url}/{relative}", cache=self.cache, timeout=self.timeout)
        return data.decode("utf-8")

    def properties(self) -> dict[str, Any]:
        """The catalogue's properties as a dict of strings (hats.properties, else the older `properties` file).

        Adds two derived keys: "hats_npix_suffix" (defaulting to ".parquet" when the file omits it) and
        "pixel_is_directory" (True when the suffix is "/", meaning each partition is a directory of part files).
        """
        if self._properties is None:
            try:
                text = self._get_text("hats.properties")
            except FetchError:
                text = self._get_text("properties")
            props: dict[str, Any] = {}
            for line in text.splitlines():
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                props[key.strip()] = value.strip().replace("\\:", ":")
            props.setdefault("hats_npix_suffix", ".parquet")
            props["pixel_is_directory"] = props["hats_npix_suffix"] == "/"
            self._properties = props
        return self._properties

    def partitions(self) -> list[tuple[int, int, int]]:
        """(order, pixel, rows) for every partition, from partition_info.csv; rows is -1 when the file omits it."""
        if self._partitions is None:
            self._partitions = parse_partition_info(self._get_text("partition_info.csv"))
            self._by_order = None
        return self._partitions

    def _index(self) -> dict[int, np.ndarray]:
        if self._by_order is None:
            by_order: dict[int, list[int]] = {}
            for order, pixel, _ in self.partitions():
                by_order.setdefault(order, []).append(pixel)
            self._by_order = {o: np.array(sorted(p), dtype=np.int64) for o, p in by_order.items()}
        return self._by_order

    def partitions_for_cells(self, order: int, pixels: Iterable[int]) -> list[Cell]:
        """Partitions (order, pixel) that overlap any of the given cells at `order`, sorted."""
        return partitions_for_cells(self._index(), order, pixels)

    def partitions_for_cone(self, ra: float, dec: float, radius_arcsec: float) -> list[Cell]:
        """Partitions (order, pixel) that can hold rows within radius_arcsec of (ra, dec) in degrees."""
        return self.partitions_for_cells(*cone_cells(ra, dec, radius_arcsec))

    def partition_files(self, order: int, pixel: int) -> list[str]:
        """HTTPS URLs of the Parquet file(s) of one partition.

        For a ".parquet" suffix this is a single file. For a directory suffix ("/") the part files are listed
        with the S3 ListObjectsV2 API when the host is an S3 bucket.
        """
        key = (order, pixel)
        if key not in self._files:
            stem = f"{self.base_url}/dataset/Norder={order}/Dir={(pixel // 10000) * 10000}/Npix={pixel}"
            if self.properties()["pixel_is_directory"]:
                self._files[key] = self._list_directory(stem + "/")
            else:
                self._files[key] = [stem + self.properties()["hats_npix_suffix"]]
        return self._files[key]

    def _list_directory(self, url: str) -> list[str]:
        parts = urllib.parse.urlsplit(url)
        if not parts.netloc.endswith(".s3.amazonaws.com"):
            raise FetchError(f"cannot list partition directory {url}: not an S3 bucket host")
        prefix = parts.path.lstrip("/")
        files: list[str] = []
        token = None
        while True:
            query = {"list-type": "2", "prefix": prefix}
            if token:
                query["continuation-token"] = token
            listing = f"{parts.scheme}://{parts.netloc}/?{urllib.parse.urlencode(query)}"
            data, _ = fetch_bytes(listing, cache=self.cache, timeout=self.timeout)
            keys, token = parse_s3_listing(data)
            files += [f"{parts.scheme}://{parts.netloc}/{k}" for k in keys if k.endswith(".parquet")]
            if not token:
                break
        if not files:
            raise FetchError(f"no Parquet files found under {url}")
        return sorted(files)

    def empty_table(self, columns: list[str]) -> pa.Table:
        """A zero-row table with the schema of the given columns, from dataset/_common_metadata."""
        key = tuple(columns)
        if key not in self._empty:
            data, _ = fetch_bytes(f"{self.base_url}/dataset/_common_metadata", cache=self.cache, timeout=self.timeout)
            table = pq.ParquetFile(pa.BufferReader(data)).read(columns=list(columns))
            self._empty[key] = table.append_column(SEPARATION_COLUMN, pa.array([], type=pa.float64())) \
                if SEPARATION_COLUMN not in table.column_names else table
        return self._empty[key]

    # ---- Parquet access

    def _open(self, url: str, stats: ReadStats, metadata: pq.FileMetaData | None = None,
              size: int | None = None) -> tuple[pq.ParquetFile, Any]:
        """Open one remote Parquet file, reusing a parsed or cached footer when there is one."""
        footer_key, size_key = f"hats-footer:{url}", f"hats-size:{url}"
        first_open = metadata is None
        if metadata is None and self.cache is not None and (blob := self.cache.get(footer_key)) is not None:
            metadata = pq.read_metadata(pa.BufferReader(blob))
            stats.add(cache_hits=1)
        if size is None and self.cache is not None and (blob := self.cache.get(size_key)) is not None:
            size = int(blob)   # with the size known, fsspec opens the file without a HEAD request
        try:
            # cache_type "none": every read is one exact range request (pyarrow's pre-buffering coalesces them).
            handle = self._fs.open(url, mode="rb", block_size=1 << 20, cache_type="none", size=size)
        except Exception as error:  # fsspec raises a variety of aiohttp and OS errors
            raise FetchError(f"could not open {url}: {error}") from error
        if self.cache is not None and url not in self._sizes and not self.cache.has(size_key):
            self.cache.put(size_key, str(handle.size).encode())
        self._sizes[url] = handle.size
        counted = _CountingFile(handle, stats, url)
        try:
            parquet = pq.ParquetFile(counted, metadata=metadata, pre_buffer=True)
        except Exception as error:
            counted.close()
            raise FetchError(f"could not read the Parquet footer of {url}: {error}") from error
        if self.cache is not None and first_open and metadata is None:
            sink = pa.BufferOutputStream()
            parquet.metadata.write_metadata_file(sink)
            self.cache.put(footer_key, sink.getvalue().to_pybytes())
        if first_open:
            stats.add(files=1)
        return parquet, counted

    def _read_group(self, parquet: pq.ParquetFile, url: str, group: int, columns: list[str],
                    stats: ReadStats, cache_key: str | None = None) -> pa.Table:
        if cache_key is not None and self.cache is not None and (blob := self.cache.get(cache_key)) is not None:
            stats.add(cache_hits=1)
            return _table_from_bytes(blob)
        try:
            table = parquet.read_row_groups([group], columns=columns, use_threads=True)
        except Exception as error:
            raise FetchError(f"could not read row group {group} of {url}: {error}") from error
        if cache_key is not None and self.cache is not None:
            self.cache.put(cache_key, _table_to_bytes(table))
        return table

    def _each_group(self, url: str, parquet: pq.ParquetFile, jobs: list[Any], work: Any,
                    stats: ReadStats) -> list[Any]:
        """Run work(parquet_file, job) for each job, up to self.workers at once, each worker on its own handle.

        The footer already parsed in `parquet` is reused, so no worker downloads it again. Results keep the
        order of jobs.
        """
        if not jobs:
            return []
        size = self._sizes.get(url)

        def run(job: Any) -> Any:
            own, handle = self._open(url, stats, metadata=parquet.metadata, size=size)
            try:
                return work(own, job)
            finally:
                handle.close()

        if self.workers == 1 or len(jobs) == 1:
            return [run(job) for job in jobs]
        with futures.ThreadPoolExecutor(max_workers=min(self.workers, len(jobs))) as pool:
            return list(pool.map(run, jobs))

    # ---- queries

    def read_cone(self, ra: float, dec: float, radius_arcsec: float, columns: list[str]) -> pa.Table:
        """Rows within radius_arcsec of (ra, dec) (degrees), with the requested columns plus the position columns
        and a float64 column "_separation_arcsec".

        Near pixel edges a cone can touch more than one order-9 pixel or partition; the cells searched are those
        touched by the centre and by four points one radius away (as the first-light prototype did), plus the
        conservative HEALPix cone cover, and every partition overlapping them is read. Row groups are pruned
        with the `_healpix_9` / `_healpix_29` statistics; if neither has statistics, all row groups of the
        partition are scanned. Nested struct sub-fields can be requested as "struct.field".
        """
        return self.read_positions([("_cone", ra, dec)], radius_arcsec, columns)["_cone"]

    def read_positions(self, positions: list[tuple[str, float, float]], radius_arcsec: float,
                       columns: list[str]) -> dict[str, pa.Table]:
        """Cone searches around many (name, ra, dec) targets at once; returns {name: table} as for read_cone.

        Targets are grouped by partition; each partition's footer is read once, the position columns of each
        candidate row group are read once and matched against every target that can fall in it, and the
        requested columns of a row group are downloaded at most once, and only if some target matched there.
        A row within the radius of two targets appears in both tables. Costs are logged at INFO and kept in
        `self.last_stats`.
        """
        names = [p[0] for p in positions]
        if len(set(names)) != len(names):
            raise ValueError("target names must be unique")
        started = time.time()
        stats = ReadStats()
        wanted = _unique([*columns, self.ra_column, self.dec_column])
        plan: dict[Cell, list[tuple[str, float, float, set[Cell]]]] = {}
        for name, ra, dec in positions:
            order, pixels = cone_cells(ra, dec, radius_arcsec)
            cells = {(order, p) for p in pixels}
            for part in self.partitions_for_cells(order, pixels):
                plan.setdefault(part, []).append((name, ra, dec, cells))
        found: dict[str, list[pa.Table]] = {name: [] for name in names}
        for part, targets in sorted(plan.items()):
            for url in self.partition_files(*part):
                self._match_file(url, targets, radius_arcsec, wanted, stats, found)
        stats.seconds = time.time() - started
        self.last_stats = stats
        log.info("HATS %s: %d target(s), %s", self.base_url.rsplit("/", 1)[-1], len(positions), stats.summary())
        out = {}
        for name in names:
            tables = found[name]
            out[name] = pa.concat_tables(tables) if tables else self.empty_table(wanted)
        return out

    def _match_file(self, url: str, targets: list[tuple[str, float, float, set[Cell]]], radius_arcsec: float,
                    wanted: list[str], stats: ReadStats, found: dict[str, list[pa.Table]]) -> None:
        """Stage 1 (positions of candidate row groups, matched to targets), then stage 2 (full columns of the
        row groups that hold a match, in parallel). Appends each target's matched rows to found[name]."""
        parquet, handle = self._open(url, stats)
        jobs = []
        try:
            groups_for = {t[0]: set(prune_row_groups(parquet.metadata, t[3])) for t in targets}
            candidate = sorted(set().union(*groups_for.values()))
            stats.add(row_groups_total=parquet.metadata.num_row_groups, row_groups_scanned=len(candidate))
            position_columns = [self.ra_column, self.dec_column]
            for group in candidate:
                pos = self._read_group(parquet, url, group, position_columns, stats, f"hats-pos:{url}#rg={group}")
                ra_all = pos.column(self.ra_column).to_numpy(zero_copy_only=False)
                dec_all = pos.column(self.dec_column).to_numpy(zero_copy_only=False)
                hits = []
                for name, ra, dec, _ in (t for t in targets if group in groups_for[t[0]]):
                    sep = angular_separation_arcsec(ra, dec, ra_all, dec_all)
                    index = np.flatnonzero(sep <= radius_arcsec)
                    if index.size:
                        hits.append((name, ra, dec, index, sep[index]))
                if hits:
                    jobs.append((group, hits))
        finally:
            handle.close()

        def work(own: pq.ParquetFile, job: tuple[int, list]) -> list[pa.Table]:
            group, hits = job
            key = None
            if self.cache is not None:
                tag = ";".join(f"{r!r},{d!r}" for _, r, d, _, _ in hits)
                key = f"hats-rows:{url}#rg={group}#cols={','.join(wanted)}#r={radius_arcsec!r}#{tag}"
            return self._matched_rows(own, url, group, wanted, hits, stats, key)

        for (group, hits), rows in zip(jobs, self._each_group(url, parquet, jobs, work, stats)):
            for (name, *_), table in zip(hits, rows):
                found[name].append(table)

    def _matched_rows(self, parquet: pq.ParquetFile, url: str, group: int, wanted: list[str],
                      hits: list[tuple[str, float, float, np.ndarray, np.ndarray]], stats: ReadStats,
                      key: str | None) -> list[pa.Table]:
        """The matched rows of one row group for each hit, reading the full columns once (or from the cache)."""
        if key is not None and self.cache is not None and (blob := self.cache.get(key)) is not None:
            stats.add(cache_hits=1)
            combined = _table_from_bytes(blob)
        else:
            full = self._read_group(parquet, url, group, wanted, stats)
            stats.add(row_groups_read=1)
            pieces = []
            for n, (_, _, _, index, sep) in enumerate(hits):
                piece = full.take(pa.array(index)).append_column(SEPARATION_COLUMN, pa.array(sep))
                pieces.append(piece.append_column("_hit", pa.array(np.full(index.size, n, dtype=np.int32))))
            del full
            combined = pa.concat_tables(pieces)
            if key is not None and self.cache is not None:
                self.cache.put(key, _table_to_bytes(combined))
        hit = combined.column("_hit")
        combined = combined.drop_columns(["_hit"])
        return [combined.filter(pc.equal(hit, n)) for n in range(len(hits))]

    def read_pixel(self, order: int, pixel: int, columns: list[str], row_filter: pc.Expression | None = None,
                   bounds: dict[str, tuple[Any, Any]] | None = None) -> pa.Table:
        """All rows inside the HEALPix cell (order, pixel), with the requested columns.

        Works whether the cell is finer than, equal to, or coarser than the catalogue's partitions: row groups
        are pruned with the pixel statistics and, when a partition is coarser than the cell, rows are cut to the
        cell with `_healpix_29`. row_filter is a pyarrow compute expression applied to each row group after
        reading; every column it uses must be in `columns`. bounds ({column: (low, high)}) prunes row groups
        whose statistics cannot satisfy it (for example {"phot_g_mean_mag": (None, 19.0)}). Row groups are read
        and filtered as they arrive (self.workers at a time), so memory stays near that many row groups.
        """
        started = time.time()
        stats = ReadStats()
        cell = (order, pixel)
        parts = self.partitions_for_cells(order, [pixel])
        cut = any(o < order for o, _ in parts)
        read_columns = list(columns)
        added = cut and "_healpix_29" not in read_columns
        if added:
            read_columns.append("_healpix_29")
        lo29, hi29 = cell_range(order, pixel, MAX_ORDER)

        def work(url: str, own: pq.ParquetFile, group: int) -> pa.Table:
            key = None
            if self.cache is not None:
                key = f"hats-pixel:{url}#rg={group}#cell={order}/{pixel}#cols={','.join(read_columns)}#f={row_filter}"
                if (blob := self.cache.get(key)) is not None:
                    stats.add(cache_hits=1)
                    return _table_from_bytes(blob)
            table = self._read_group(own, url, group, read_columns, stats)
            stats.add(row_groups_read=1)
            if cut:
                h = table.column("_healpix_29")
                table = table.filter(pc.and_(pc.greater_equal(h, lo29), pc.less_equal(h, hi29)))
            if row_filter is not None:
                table = table.filter(row_filter)
            if key is not None:
                self.cache.put(key, _table_to_bytes(table))
            return table

        tables = []
        for part in parts:
            for url in self.partition_files(*part):
                parquet, handle = self._open(url, stats)
                handle.close()
                groups = prune_row_groups(parquet.metadata, [cell], bounds)
                stats.add(row_groups_total=parquet.metadata.num_row_groups, row_groups_scanned=len(groups))
                tables += self._each_group(url, parquet, groups, lambda own, g, u=url: work(u, own, g), stats)
        stats.seconds = time.time() - started
        self.last_stats = stats
        log.info("HATS %s pixel %d/%d: %s", self.base_url.rsplit("/", 1)[-1], order, pixel, stats.summary())
        if not tables:
            return self.empty_table(read_columns).drop_columns([SEPARATION_COLUMN])
        out = pa.concat_tables(tables)
        return out.drop_columns(["_healpix_29"]) if added else out


# ----------------------------------------------------------------------------------------------- parsers


def parse_partition_info(text: str) -> list[tuple[int, int, int]]:
    """(order, pixel, rows) tuples from the text of a HATS partition_info.csv (rows is -1 if not given)."""
    reader = csv.DictReader(io.StringIO(text))
    fields = {name.strip(): name for name in (reader.fieldnames or [])}
    if "Norder" not in fields or "Npix" not in fields:
        raise ValueError(f"partition_info.csv lacks Norder/Npix columns (has {reader.fieldnames})")
    rows_field = fields.get("num_rows")
    out = []
    for row in reader:
        rows = int(row[rows_field]) if rows_field and row.get(rows_field) not in (None, "") else -1
        out.append((int(row[fields["Norder"]]), int(row[fields["Npix"]]), rows))
    return out


def partitions_for_cells(index: dict[int, np.ndarray], order: int, pixels: Iterable[int]) -> list[Cell]:
    """Partitions overlapping any cell, given {partition order: sorted pixel array}: coarser partitions that
    contain a cell, and finer partitions inside one."""
    found: set[Cell] = set()
    for pixel in pixels:
        for part_order, part_pixels in index.items():
            lo, hi = cell_range(order, int(pixel), part_order)
            start = np.searchsorted(part_pixels, lo, side="left")
            stop = np.searchsorted(part_pixels, hi, side="right")
            found.update((part_order, int(p)) for p in part_pixels[start:stop])
    return sorted(found)


def parse_s3_listing(data: bytes) -> tuple[list[str], str | None]:
    """Object keys and the continuation token (None when complete) from an S3 ListObjectsV2 XML response."""
    root = ElementTree.fromstring(data)
    ns = root.tag[: root.tag.index("}") + 1] if root.tag.startswith("{") else ""
    keys = [el.text or "" for el in root.iter(f"{ns}Key")]
    truncated = (root.findtext(f"{ns}IsTruncated") or "false").lower() == "true"
    token = root.findtext(f"{ns}NextContinuationToken") if truncated else None
    return keys, token
