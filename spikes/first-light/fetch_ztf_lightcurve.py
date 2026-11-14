#!/usr/bin/env python3
"""Pull one object's ZTF DR24 light curves from the public cloud copy.

The light curves live on Amazon S3 as a HATS catalogue (Parquet files partitioned by
HEALPix sky pixel). This script finds the partition that covers the target, reads the
Parquet footer, and downloads only the row groups whose sky-pixel range includes the
target, using HTTP range requests. No account or API key is needed.

    python fetch_ztf_lightcurve.py pg1302 196.38754 -10.555389
"""
from __future__ import annotations

import argparse
import io
import json
import pathlib
import time
import urllib.request

import fsspec
import hats.pixel_math.healpix_shim as healpix
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BASE = "https://ipac-irsa-ztf.s3.amazonaws.com/ztf/enhanced/dr24/lc/hats/ztf_dr24_lc-hats"
CITATION = "ZTF DR24 light curves, IRSA, doi:10.26131/IRSA598"
BANDS = {1: "g", 2: "r", 3: "i"}
COLUMNS = ["objectid", "filterid", "objra", "objdec", "nepochs", "lightcurve"]


def nested_pixel(order: int, ra: float, dec: float) -> int:
    return int(np.atleast_1d(healpix.radec2pix(order, [ra], [dec]))[0])


def wanted_pixels(ra: float, dec: float, radius_arcsec: float) -> set[int]:
    """Order-9 pixels touched by the target and by four points one search radius away."""
    step = radius_arcsec / 3600.0
    stretch = max(np.cos(np.radians(dec)), 1e-6)
    points = [(ra, dec), (ra + step / stretch, dec), (ra - step / stretch, dec), (ra, dec + step), (ra, dec - step)]
    return {nested_pixel(9, r, d) for r, d in points}


def partitions_for(pixels: set[int]) -> set[tuple[int, int]]:
    """Map order-9 pixels to the (Norder, Npix) partitions that hold them."""
    with urllib.request.urlopen(f"{BASE}/partition_info.csv", timeout=60) as response:
        info = pd.read_csv(io.BytesIO(response.read()))
    existing = set(zip(info.Norder, info.Npix))
    found = set()
    for pixel in pixels:
        for order in sorted(info.Norder.unique()):
            candidate = (int(order), pixel >> (2 * (9 - int(order))))
            if candidate in existing:
                found.add(candidate)
    return found


def read_partition(order: int, npix: int, pixels: set[int]) -> pd.DataFrame:
    url = f"{BASE}/dataset/Norder={order}/Dir={(npix // 10000) * 10000}/Npix={npix}/part0.snappy.parquet"
    filesystem = fsspec.filesystem("https", client_kwargs={"trust_env": True})
    with filesystem.open(url, block_size=4 * 1024 * 1024) as handle:
        parquet = pq.ParquetFile(handle)
        column = parquet.schema_arrow.get_field_index("_healpix_9")
        groups = []
        for index in range(parquet.metadata.num_row_groups):
            stats = parquet.metadata.row_group(index).column(column).statistics
            if stats is None or any(stats.min <= pixel <= stats.max for pixel in pixels):
                groups.append(index)
        print(f"  Norder={order} Npix={npix}: reading {len(groups)} of {parquet.metadata.num_row_groups} row groups")
        return parquet.read_row_groups(groups, columns=COLUMNS).to_pandas()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("name")
    parser.add_argument("ra", type=float, help="right ascension, degrees (J2000)")
    parser.add_argument("dec", type=float, help="declination, degrees (J2000)")
    parser.add_argument("--radius", type=float, default=2.0, help="match radius, arcseconds")
    parser.add_argument("--out", type=pathlib.Path, default=pathlib.Path(__file__).parent / "data")
    args = parser.parse_args()

    started = time.time()
    pixels = wanted_pixels(args.ra, args.dec, args.radius)
    frames = [read_partition(order, npix, pixels) for order, npix in sorted(partitions_for(pixels))]
    table = pd.concat(frames, ignore_index=True)
    separation = np.hypot((table.objra - args.ra) * np.cos(np.radians(args.dec)), table.objdec - args.dec) * 3600.0
    matches = table[separation < args.radius]

    curves = []
    for _, row in matches.iterrows():
        curve = {key: np.asarray(row["lightcurve"][key]).tolist() for key in ("hmjd", "mag", "magerr", "catflags")}
        good = int(np.sum(np.asarray(curve["catflags"]) == 0))
        span = (max(curve["hmjd"]) - min(curve["hmjd"])) / 365.25
        print(f"  {BANDS[int(row.filterid)]} band: object {int(row.objectid)}, {int(row.nepochs)} epochs "
              f"({good} unflagged) over {span:.2f} years")
        curves.append({"objectid": int(row.objectid), "band": BANDS[int(row.filterid)], **curve})

    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / f"{args.name}_ztf_dr24.json"
    target.write_text(json.dumps({"name": args.name, "ra": args.ra, "dec": args.dec, "source": CITATION,
                                  "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                  "curves": curves}))
    print(f"wrote {target} ({len(curves)} light curves) in {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
