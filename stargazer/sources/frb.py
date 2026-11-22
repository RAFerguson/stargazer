"""Arrival times of bursts from repeating fast radio bursts, from burst tables deposited on Zenodo.

Written against the deposits as they stood on 22 Nov 2026. Of the candidates, CHIME/FRB's own site could not be
used: every page and path on www.chime-frb.ca answers with a redirect to chime-frb.ca, which this machine's proxy
refuses. Three Zenodo deposits work (files are fetched from zenodo.org/api/records/<id>/files/<name>/content, the
path Zenodo's robots.txt allows; requests to Zenodo are spaced by its 10 s crawl delay):

``FRB 20220912A``, Nançay Radio Telescope (Konijn et al. 2024), record 13889738, CC-BY-4.0.
    696 bursts in 68 scans (60.99 h) between MJD 59867 and 60054, 1.228-1.740 GHz, with the scan log, so with
    exposure. Times here are the topocentric (UTC) file start of each candidate's filterbank plus Heimdall's
    candidate time, both read from the candidate file name in ``NRT_bursts_with_corrected_drift.csv``. They are not
    dedispersed to infinite frequency: they refer to Heimdall's reference frequency at each burst's detection DM
    (believed to be the top of the band, 1.740 GHz; the deposit does not say). With ``barycentre=True`` times and
    scans are moved to the Solar System barycentre (TDB) with astropy's built-in ephemeris and the EVN position.
    The deposit's own "Time of arrival (MJD)" column (in ``Supplementary Material (online).csv``) is rounded to
    1e-5 d (0.86 s) and is mixed: 681 rows are topocentric, 15 rows (the paper's Table 2) are barycentric TDB at
    infinite frequency, and for 6 of those 15 the two files disagree by 1.5 to 61 minutes (B384 and B411 appear
    swapped). Every row is classified in ``meta["supplementary_toa_check"]``; the paper's printed Table 1 (arXiv
    v1) lists file start times rather than burst times and is not used.

``FRB 20180916B``, compilation by Pastor-Marazuela et al. (2021), record 4559593, Apache-2.0.
    ``files/r3all_data.json`` in the reproduction package: 151 burst times from nine instruments (CHIME/FRB 56,
    Apertif 52, GMRT 650 MHz 15, LOFAR 9, GBT 300 MHz 7, Effelsberg 4, GMRT 400 MHz 4, SRT P band 3, VLA 1) between
    MJD 58377 and 59144, S/N per burst, and observation start times and durations for 17 instrument set-ups
    including several that saw nothing. The Apertif and LOFAR times match the paper's barycentric tables; the
    others are as the authors compiled them from the literature, and the package does not state their time
    standard or reference frequency. The windows are incomplete (17 CHIME, 3 Apertif and the VLA burst lie
    outside every listed window of their instrument) and their starts are topocentric, so they are a guide to
    exposure, not a complete record. One CHIME time appears twice and is kept once.

``FRB 20220529``, FAST (Li et al. 2026, Science 391, 280), record 17446358, CC-BY-4.0.
    Supplementary Table S5: 643 bursts with S/N >= 10 and a rotation measure, outside the "RM flare", between MJD
    59758 and 60557, barycentric TDB at infinite frequency, to 1e-7 d (8.6 ms). A selected subset (bursts during
    the flare are not in it) and with no exposure information. Two pairs of rows share one MJD (different S/N and
    RM: components or distinct bursts closer than the rounding); both are kept and listed in meta.

`fetch_bursts` returns an `EventList` with time in seconds from the first burst; the first burst's MJD, the time
standard and everything per burst (ids, DM, S/N, instrument) are in ``meta``, aligned with the sorted times.
Intervals are formed from MJD differences times 86400 s; no leap second falls inside any of these spans.
"""
from __future__ import annotations

import csv
import io
import json
import logging
import re
import time as _time
import urllib.parse
import zipfile
from typing import Any

import numpy as np

from ..core import EventList, Provenance, content_hash
from .cache import Cache, cache_key, fetch_bytes, utc_now

log = logging.getLogger(__name__)

ZENODO_FILE = "https://zenodo.org/api/records/{record}/files/{name}/content"
MIN_INTERVAL_S = {"zenodo.org": 10.0}          # Zenodo's robots.txt crawl delay; 1 s for any other host
SECONDS_PER_DAY = 86400.0

NRT_CITATION = (
    "Konijn, D. C., Hewitt, D. M., Hessels, J. W. T., et al. 2024, \"A Nançay Radio Telescope study of the "
    "hyperactive repeating FRB 20220912A\", MNRAS 534, 3331-3348, doi:10.1093/mnras/stae2296, arXiv:2407.10155. "
    "Data: Konijn & Hewitt, reproduction package, Zenodo, doi:10.5281/zenodo.13889738 (CC-BY-4.0). The Nançay "
    "Radio Observatory is operated by the Paris Observatory, associated with the French CNRS."
)
R3_CITATION = (
    "Pastor-Marazuela, I., Connor, L., van Leeuwen, J., et al. 2021, \"Chromatic periodic activity down to 120 "
    "megahertz in a fast radio burst\", Nature 596, 505-508, doi:10.1038/s41586-021-03724-8, arXiv:2012.08348. "
    "Data: reproduction package (files/r3all_data.json), Zenodo, doi:10.5281/zenodo.4559593 (Apache-2.0), which "
    "compiles burst times and observing windows of FRB 20180916B from Apertif, LOFAR and published observations by "
    "CHIME/FRB, Effelsberg, GBT, uGMRT, SRT, VLA, DSN and others; cite the original papers for those."
)
FAST_20220529_CITATION = (
    "Li, Y., Zhang, S. B., Yang, Y. P., et al. 2026, \"A sudden change and recovery in the magnetic environment "
    "around a repeating fast radio burst\", Science 391, 280-284, doi:10.1126/science.adq3225, arXiv:2503.04727. "
    "Data: Supplementary Table S5, Zenodo, doi:10.5281/zenodo.17446358 (CC-BY-4.0). Based on observations with "
    "the Five-hundred-meter Aperture Spherical radio Telescope (FAST)."
)

# Nançay: ITRF position of the NRT (tempo2 observatory "NANCAY"), and the EVN position of FRB 20220912A
# (Hewitt et al. 2024, as quoted by Konijn et al. 2024).
NRT_ITRF_M = (4324165.81, 165927.11, 4670132.83)
FRB20220912A_RA, FRB20220912A_DEC = "23h09m04.8989s", "+48d42m23.9078s"
NRT_SOURCE_TAG = "2309+4842"

DATASETS: dict[str, dict[str, Any]] = {
    "FRB 20220912A": {
        "registry_id": "frb_nrt_20220912a", "record": 13889738, "doi": "10.5281/zenodo.13889738",
        "files": {"bursts": "NRT_bursts_with_corrected_drift.csv",
                  "supplementary": "Supplementary Material (online).csv", "scans": "Observation times.txt"},
        "telescope": "Nançay Radio Telescope (NUPPI, 1.228-1.740 GHz)", "licence": "CC-BY-4.0",
        "n_bursts": 696, "mjd_range": (59867.90, 60054.39), "exposure": True, "citation": NRT_CITATION,
        "time_standard": "topocentric UTC at Nançay, Heimdall reference frequency (not dedispersed); "
                         "barycentric TDB with barycentre=True",
    },
    "FRB 20180916B": {
        "registry_id": "frb_r3_compilation", "record": 4559593, "doi": "10.5281/zenodo.4559593",
        "files": {"package": "Apertif_LOFAR_R3.zip"}, "member": "files/r3all_data.json",
        "telescope": "compilation: CHIME/FRB, Apertif, LOFAR, uGMRT, GBT, Effelsberg, SRT, VLA",
        "licence": "Apache-2.0", "n_bursts": 150, "mjd_range": (58377.43, 59143.82), "exposure": True,
        "citation": R3_CITATION,
        "time_standard": "MJD as compiled (Apertif and LOFAR barycentric; others as published, standard and "
                         "reference frequency not stated in the package)",
    },
    "FRB 20220529": {
        "registry_id": "frb_fast_20220529", "record": 17446358, "doi": "10.5281/zenodo.17446358",
        "files": {"table": "FRB20220529_TableS5.csv"},
        "telescope": "FAST (19-beam receiver, 1.0-1.5 GHz)", "licence": "CC-BY-4.0",
        "n_bursts": 643, "mjd_range": (59758.93, 60556.87), "exposure": False, "citation": FAST_20220529_CITATION,
        "time_standard": "barycentric TDB, extrapolated to infinite frequency",
    },
}


def canonical_name(source_name: str) -> str:
    """The DATASETS key for a name such as "FRB 20220912A", "FRB20220912A" or "20220912a"."""
    m = re.fullmatch(r"\s*(?:FRB)?\s*_?(\d{8})([A-Za-z]{0,3})\s*", source_name, flags=re.I)
    if m:
        name = f"FRB {m[1]}{m[2].upper()}"
        if name in DATASETS:
            return name
    raise KeyError(f"no burst list for {source_name!r}; available: {sorted(DATASETS)}")


def zenodo_url(record: int, name: str) -> str:
    return ZENODO_FILE.format(record=record, name=urllib.parse.quote(name))


def list_sources() -> list[dict[str, Any]]:
    """The burst lists this module can fetch: name, number of bursts, MJD range, telescope, where from, licence.

    Counts are those of the deposits as read on 22 Nov 2026 (FRB 20180916B after removing one duplicated time).
    """
    out = []
    for name, d in DATASETS.items():
        out.append({"name": name, "n_bursts": d["n_bursts"], "mjd_range": d["mjd_range"],
                    "span_days": round(d["mjd_range"][1] - d["mjd_range"][0], 2), "telescope": d["telescope"],
                    "where": f"https://doi.org/{d['doi']}", "urls": [zenodo_url(d["record"], f)
                                                                     for f in d["files"].values()],
                    "licence": d["licence"], "time_standard": d["time_standard"], "has_exposure": d["exposure"],
                    "registry_id": d["registry_id"], "citation": d["citation"]})
    return out


# ------------------------------------------------------------------------------------------------- download

_last_request: dict[str, float] = {}


def _get(url: str, cache: Cache | None) -> tuple[bytes, str]:
    """fetch_bytes, waiting first if the host was asked less than its minimum interval ago (cache hits never wait)."""
    if cache is None or not cache.has(cache_key(url)):
        host = urllib.parse.urlsplit(url).hostname or ""
        wait = MIN_INTERVAL_S.get(host, 1.0) - (_time.monotonic() - _last_request.get(host, -1e9))
        if wait > 0:
            _time.sleep(wait)
        _last_request[host] = _time.monotonic()
    return fetch_bytes(url, cache=cache, timeout=120.0)


# ------------------------------------------------------------------------------------------------- parsing


def _split_mjd(text: str) -> tuple[int, float]:
    """An MJD string as (whole days, fraction of a day), so that no precision is lost to one float."""
    text = text.strip()
    if not re.fullmatch(r"\d+(\.\d*)?", text):
        raise ValueError(f"not an MJD: {text!r}")
    whole, _, frac = text.partition(".")
    return int(whole), float("0." + (frac or "0"))


def _seconds_since(days: np.ndarray, frac: np.ndarray, day0: int, frac0: float) -> np.ndarray:
    return (np.asarray(days, dtype=np.int64) - day0) * SECONDS_PER_DAY + (np.asarray(frac) - frac0) * SECONDS_PER_DAY


def _float(text: str) -> float:
    text = text.strip()
    return float(text) if text else float("nan")


_CAND = re.compile(r"tstart_(\d+\.\d+)_tcand_(\d+(?:\.\d+)?)_dm_(\d+(?:\.\d+)?)_snr_(\d+(?:\.\d+)?)")


def parse_nrt_bursts(text: str) -> list[dict[str, Any]]:
    """Rows of NRT_bursts_with_corrected_drift.csv, with the candidate file name split into its parts.

    Each row: burst_id, file_start (MJD string, topocentric UTC), tcand_s (Heimdall time within the file),
    heimdall_dm, snr, fluence_jyms, width_ms, freq_extent_mhz, centre_freq_ghz.
    """
    reader = csv.DictReader(io.StringIO(text))
    need = {"BurstID", "Candidate name", "Heimdall DM (pc cm^-3)", "S/N", "Fluence (Jy ms)"}
    if not reader.fieldnames or not need <= set(reader.fieldnames):
        raise ValueError(f"NRT burst table lacks columns {need - set(reader.fieldnames or [])}")
    rows = []
    for r in reader:
        m = _CAND.search(r["Candidate name"])
        if not m:
            raise ValueError(f"cannot read the candidate name of {r['BurstID']}: {r['Candidate name']!r}")
        rows.append({"burst_id": r["BurstID"].strip(), "file_start": m[1], "tcand_s": float(m[2]),
                     "heimdall_dm": _float(r["Heimdall DM (pc cm^-3)"]), "snr": _float(r["S/N"]),
                     "fluence_jyms": _float(r["Fluence (Jy ms)"]), "width_ms": _float(r.get("Temporal width (ms)", "")),
                     "freq_extent_mhz": _float(r.get("Frequency extent (MHz)", "")),
                     "centre_freq_ghz": _float(r.get("Central observed frequency (GHz)", ""))})
    if not rows:
        raise ValueError("NRT burst table has no rows")
    return rows


def parse_nrt_supplementary(text: str) -> dict[str, str]:
    """Burst id -> the "Time of arrival (MJD)" string of Supplementary Material (online).csv."""
    reader = csv.DictReader(io.StringIO(text))
    if not reader.fieldnames or "Time of arrival (MJD)" not in reader.fieldnames:
        raise ValueError("NRT supplementary table lacks 'Time of arrival (MJD)'")
    return {r["Burst ID"].strip(): r["Time of arrival (MJD)"].strip() for r in reader}


def parse_nrt_scans(text: str, tag: str = NRT_SOURCE_TAG) -> list[tuple[str, str, float]]:
    """(scan name, start MJD string, duration s) for the NUPPI scans of one source in "Observation times.txt".

    The log covers several sources; a scan belongs to FRB 20220912A when its name carries "2309+4842".
    Start times are topocentric (UTC).
    """
    scans = []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) < 3 or not fields[0].startswith("nuppi_"):
            continue
        if fields[0].split("_")[2:3] != [tag]:
            continue
        _split_mjd(fields[1])
        scans.append((fields[0], fields[1], float(fields[2])))
    scans.sort(key=lambda s: _split_mjd(s[1]))
    return scans


def _barycentre_nrt(days: np.ndarray, frac: np.ndarray):
    """Topocentric UTC MJDs at Nançay -> astropy Time in TDB at the barycentre (built-in ephemeris, no downloads)."""
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    from astropy.utils import iers

    site = EarthLocation.from_geocentric(*NRT_ITRF_M, unit=u.m)
    target = SkyCoord(FRB20220912A_RA, FRB20220912A_DEC, frame="icrs")
    with iers.conf.set_temp("auto_download", False):
        t = Time(np.asarray(days, float), np.asarray(frac, float), format="mjd", scale="utc", location=site)
        return t.tdb + t.light_travel_time(target, ephemeris="builtin")


def _mjd_text(t) -> str:
    return str(t.to_value("mjd", subfmt="str"))


def build_nrt(bursts_text: str, supplementary_text: str | None, scans_text: str, *, barycentre: bool = False) -> dict:
    """Times, amplitudes, exposure and meta of the FRB 20220912A NRT list (see the module docstring)."""
    rows = parse_nrt_bursts(bursts_text)
    scans = parse_nrt_scans(scans_text)
    if not scans:
        raise ValueError("no FRB 20220912A scans in the observation log")
    start = [_split_mjd(r["file_start"]) for r in rows]
    days = np.array([s[0] for s in start], dtype=np.int64)
    frac = np.array([s[1] for s in start]) + np.array([r["tcand_s"] for r in rows]) / SECONDS_PER_DAY
    days, frac = days + np.floor(frac).astype(np.int64), frac - np.floor(frac)
    order = np.lexsort((frac, days))
    rows = [rows[i] for i in order]
    days, frac = days[order], frac[order]

    scan_start = [_split_mjd(s[1]) for s in scans]
    sdays = np.array([s[0] for s in scan_start], dtype=np.int64)
    sfrac = np.array([s[1] for s in scan_start])
    sdur = np.array([s[2] for s in scans])

    bary = _barycentre_nrt(days, frac)
    if barycentre:
        bary_scan_start = _barycentre_nrt(sdays, sfrac)
        bary_scan_stop = _barycentre_nrt(sdays, sfrac + sdur / SECONDS_PER_DAY)
        t = (bary - bary[0]).to_value("s")
        exposure = [(float(a), float(b)) for a, b in zip((bary_scan_start - bary[0]).to_value("s"),
                                                         (bary_scan_stop - bary[0]).to_value("s"))]
        epoch = _mjd_text(bary[0])
        standard = ("barycentric TDB (astropy built-in ephemeris, NRT ITRF position, EVN source position); "
                    "not dedispersed: Heimdall reference frequency at each burst's detection DM")
    else:
        t = _seconds_since(days, frac, int(days[0]), float(frac[0]))
        s0 = _seconds_since(sdays, sfrac, int(days[0]), float(frac[0]))
        exposure = [(float(a), float(a + d)) for a, d in zip(s0, sdur)]
        epoch = f"{int(days[0])}{f'{frac[0]:.12f}'[1:]}"
        standard = ("topocentric UTC at the Nançay Radio Telescope (NUPPI file start + Heimdall candidate time); "
                    "not dedispersed: Heimdall reference frequency at each burst's detection DM")

    check = {"agrees_topocentric": [], "supplementary_is_barycentric": [], "conflict": [], "missing": []}
    if supplementary_text is not None:
        supp = parse_nrt_supplementary(supplementary_text)
        for r, d, f, b in zip(rows, days, frac, bary):
            s = supp.get(r["burst_id"])
            if s is None:
                check["missing"].append(r["burst_id"])
                continue
            sd, sf = _split_mjd(s)
            if abs(((sd - d) + (sf - f)) * SECONDS_PER_DAY) < 1.0:
                check["agrees_topocentric"].append(r["burst_id"])
            elif abs((sd + sf - b.mjd) * SECONDS_PER_DAY) < 1.5:
                check["supplementary_is_barycentric"].append(r["burst_id"])
            else:
                check["conflict"].append({"burst_id": r["burst_id"], "supplementary_mjd": s,
                                          "candidate_mjd_topocentric": f"{int(d)}{f'{f:.9f}'[1:]}",
                                          "difference_s": round(((sd - d) + (sf - f)) * SECONDS_PER_DAY, 2)})
    inside = np.zeros(len(t), dtype=bool)
    for a, b in exposure:
        inside |= (t >= a - 1e-3) & (t <= b + 1e-3)

    meta = {
        "source_name": "FRB 20220912A", "telescope": DATASETS["FRB 20220912A"]["telescope"],
        "epoch_mjd": epoch, "time_standard": standard, "time_frame": "barycentric" if barycentre else "topocentric",
        "dispersion_reference": "Heimdall reference frequency (believed top of band, 1740 MHz; not stated), "
                                "per-burst detection DM in meta['dm']",
        "band_mhz": [1228, 1740], "sampling_us": 16,
        "amplitude": "fluence (Jy ms)", "burst_id": [r["burst_id"] for r in rows],
        "snr": [r["snr"] for r in rows], "dm": [r["heimdall_dm"] for r in rows],
        "width_ms": [r["width_ms"] for r in rows], "centre_freq_ghz": [r["centre_freq_ghz"] for r in rows],
        "scans": [{"name": n, "start_mjd_utc_topocentric": s, "duration_s": d} for n, s, d in scans],
        "exposure_hours": float(sdur.sum() / 3600.0), "n_scans": len(scans),
        "n_bursts_outside_exposure": int((~inside).sum()),
        "supplementary_toa_check": check,
        "site_itrf_m": list(NRT_ITRF_M), "source_position_icrs": [FRB20220912A_RA, FRB20220912A_DEC],
    }
    return {"time": np.asarray(t, float), "amplitude": np.array([r["fluence_jyms"] for r in rows]),
            "exposure": exposure, "meta": meta}


def parse_r3_json(payload: dict[str, Any]) -> dict[str, Any]:
    """Check the layout of r3all_data.json: bursts and snr per instrument, obs_startmjds and obs_duration."""
    for key in ("bursts", "snr", "obs_startmjds", "obs_duration"):
        if not isinstance(payload.get(key), dict):
            raise ValueError(f"R3 data lacks {key!r} (has {sorted(payload)})")
    for inst, times in payload["bursts"].items():
        if len(payload["snr"].get(inst, [])) != len(times):
            raise ValueError(f"R3 data: {inst} has {len(times)} bursts but {len(payload['snr'].get(inst, []))} S/N")
    for inst, starts in payload["obs_startmjds"].items():
        if len(payload["obs_duration"].get(inst, [])) != len(starts):
            raise ValueError(f"R3 data: {inst} start and duration lists differ in length")
    return payload


def read_r3_zip(data: bytes, member: str = "files/r3all_data.json") -> dict[str, Any]:
    """The R3 compilation from the reproduction package's zip (only that one JSON member is read)."""
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        info = archive.getinfo(member)
        if info.file_size > 5_000_000:
            raise ValueError(f"{member} is unexpectedly large ({info.file_size} bytes)")
        return parse_r3_json(json.loads(archive.read(member)))


def build_r3(payload: dict[str, Any], instrument: str | None = None) -> dict:
    """Times, S/N, exposure and meta of FRB 20180916B from the compilation, for one instrument or all."""
    insts = sorted(payload["bursts"]) if instrument is None else [instrument]
    if instrument is not None and instrument not in payload["bursts"]:
        raise KeyError(f"no FRB 20180916B bursts from {instrument!r}; instruments: {sorted(payload['bursts'])}")
    events, dropped = [], []
    for inst in insts:
        seen: set[float] = set()
        for mjd, snr in zip(payload["bursts"][inst], payload["snr"][inst]):
            mjd = float(mjd)
            if mjd in seen:
                dropped.append({"instrument": inst, "mjd": mjd})
                continue
            seen.add(mjd)
            # GMRT 650 MHz S/N values in the file are all 1: placeholders, not measurements.
            amp = float("nan") if (inst == "GMRT650" or snr is None) else float(snr)
            events.append((mjd, inst, amp))
    if not events:
        raise ValueError("no bursts")
    events.sort()
    mjd = np.array([e[0] for e in events])
    day0, frac0 = int(np.floor(mjd[0])), float(mjd[0] - np.floor(mjd[0]))
    t = _seconds_since(np.floor(mjd).astype(np.int64), mjd - np.floor(mjd), day0, frac0)

    window_insts = sorted(payload["obs_startmjds"]) if instrument is None else (
        [instrument] if instrument in payload["obs_startmjds"] else [])
    windows = []
    per_inst = {}
    for inst in window_insts:
        starts = np.array(payload["obs_startmjds"][inst], float)
        dur = np.array(payload["obs_duration"][inst], float)
        per_inst[inst] = {"n_windows": int(starts.size), "hours": float(dur.sum() / 3600.0)}
        s = _seconds_since(np.floor(starts).astype(np.int64), starts - np.floor(starts), day0, frac0)
        windows += [(float(a), float(a + d)) for a, d in zip(s, dur)]
    windows.sort()
    merged: list[list[float]] = []
    for a, b in windows:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    exposure = [(a, b) for a, b in merged] or None
    inside = np.zeros(len(t), dtype=bool)
    for a, b in exposure or []:
        inside |= (t >= a - 864.0) & (t <= b + 864.0)
    close = [{"instruments": [events[i][1], events[i + 1][1]], "mjd": events[i][0],
              "separation_s": float(t[i + 1] - t[i])} for i in range(len(t) - 1) if t[i + 1] - t[i] < 1.0]
    meta = {
        "source_name": "FRB 20180916B", "instrument_selected": instrument,
        "epoch_mjd": repr(float(mjd[0])), "time_standard": DATASETS["FRB 20180916B"]["time_standard"],
        "time_frame": "barycentric (as compiled)",
        "dispersion_reference": "not stated in the package; as published by each instrument's team",
        "amplitude": "detection S/N (instrument-dependent; GMRT650 placeholders set to NaN)",
        "instrument": [e[1] for e in events],
        "n_bursts_by_instrument": {i: sum(e[1] == i for e in events) for i in insts},
        "exposure_by_instrument": per_inst,
        "exposure_note": "windows are the compilation's obs_startmjds (topocentric starts) and obs_duration; they "
                         "are incomplete and not barycentred (offsets up to ~8 min)",
        "n_bursts_outside_exposure": int((~inside).sum()) if exposure else None,
        "freq_range_mhz": {k: [payload.get("freq_min", {}).get(k), payload.get("freq_max", {}).get(k)]
                           for k in payload.get("freq_min", {})},
        "duplicates_dropped": dropped, "close_pairs_under_1s": close,
    }
    return {"time": t, "amplitude": np.array([e[2] for e in events]), "exposure": exposure, "meta": meta}


def parse_fast_20220529(text: str) -> list[dict[str, Any]]:
    """Rows of FRB20220529_TableS5.csv; '#' lines are comments and '#MJD,...' is the header."""
    header = None
    rows = []
    for line in text.splitlines():
        if line.startswith("#MJD"):
            header = next(csv.reader([line[1:]]))
            continue
        if not line.strip() or line.startswith("#"):
            continue
        if header is None:
            raise ValueError("FAST table has data before its '#MJD' header")
        fields = next(csv.reader([line]))
        if len(fields) != len(header):
            raise ValueError(f"FAST table row has {len(fields)} fields, header has {len(header)}: {line!r}")
        r = dict(zip(header, fields))
        rows.append({"mjd": r["MJD"].strip(), "snr": _float(r["S/N"]), "dm_det": _float(r["DM_det (pc cm^-3)"]),
                     "f_low_mhz": _float(r["f_low (MHz)"]), "f_high_mhz": _float(r["f_high (MHz)"]),
                     "dm_stru": r.get("DM_stru (pc cm^-3)", "").strip(), "rm": r.get("RM (rad m^-2)", "").strip()})
    if header is None or not rows:
        raise ValueError("FAST table has no '#MJD' header or no rows")
    return rows


def build_fast_20220529(text: str) -> dict:
    rows = parse_fast_20220529(text)
    split = [_split_mjd(r["mjd"]) for r in rows]
    days = np.array([s[0] for s in split], dtype=np.int64)
    frac = np.array([s[1] for s in split])
    order = np.lexsort((frac, days))
    rows, days, frac = [rows[i] for i in order], days[order], frac[order]
    t = _seconds_since(days, frac, int(days[0]), float(frac[0]))
    same = [rows[i]["mjd"] for i in range(len(rows) - 1) if rows[i]["mjd"] == rows[i + 1]["mjd"]]
    meta = {
        "source_name": "FRB 20220529", "same_mjd_rows": same, "telescope": DATASETS["FRB 20220529"]["telescope"],
        "epoch_mjd": rows[0]["mjd"], "time_standard": DATASETS["FRB 20220529"]["time_standard"],
        "time_frame": "barycentric", "dispersion_reference": "infinite frequency",
        "time_resolution_s": 1e-7 * SECONDS_PER_DAY,
        "selection": "S/N >= 10 bursts with an RM measurement outside the RM flare (Table S5); not complete",
        "amplitude": "S/N (maximum over DM trials)",
        "dm_det": [r["dm_det"] for r in rows], "f_low_mhz": [r["f_low_mhz"] for r in rows],
        "f_high_mhz": [r["f_high_mhz"] for r in rows], "rm": [r["rm"] for r in rows],
        "exposure_note": "the deposit gives no observing log; exposure is unknown",
    }
    return {"time": t, "amplitude": np.array([r["snr"] for r in rows]), "exposure": None, "meta": meta}


# ------------------------------------------------------------------------------------------------- public


def fetch_bursts(source_name: str, *, cache: Cache | None = None, instrument: str | None = None,
                 barycentre: bool = False) -> EventList:
    """Burst arrival times of one repeating FRB as an EventList (time in s from the first burst).

    source_name is one of `list_sources()` names (spacing and the "FRB" prefix are optional). instrument picks
    one instrument of the FRB 20180916B compilation (e.g. "CHIME", "Apertif"); barycentre=True moves the
    topocentric FRB 20220912A times and scans to the barycentre (the other two lists are barycentric already, and
    the flag is then recorded but changes nothing). See the module docstring for time standards and caveats.
    """
    name = canonical_name(source_name)
    d = DATASETS[name]
    if instrument is not None and name != "FRB 20180916B":
        raise ValueError("instrument applies only to the FRB 20180916B compilation")
    urls = {key: zenodo_url(d["record"], fname) for key, fname in d["files"].items()}
    raw, sha = {}, {}
    for key, url in urls.items():
        raw[key], sha[key] = _get(url, cache)
    if name == "FRB 20220912A":
        built = build_nrt(raw["bursts"].decode("utf-8"), raw["supplementary"].decode("utf-8"),
                          raw["scans"].decode("utf-8", "replace"), barycentre=barycentre)
    elif name == "FRB 20180916B":
        built = build_r3(read_r3_zip(raw["package"], d["member"]), instrument)
    else:
        built = build_fast_20220529(raw["table"].decode("utf-8"))
    meta = built["meta"]
    meta.update({"n_bursts": len(built["time"]), "span_days": float(built["time"][-1] / SECONDS_PER_DAY),
                 "licence": d["licence"], "doi": d["doi"]})
    events = EventList(built["time"], amplitude=built["amplitude"], exposure=built["exposure"],
                       source_id=f"frb:{name.removeprefix('FRB ')}", time_unit="s", meta=meta)
    events.provenance = Provenance(
        d["registry_id"],
        {"source_name": name, "instrument": instrument, "barycentre": barycentre, "urls": urls, "file_sha256": sha},
        utc_now(), content_hash({"time": events.time, "amplitude": events.amplitude, "exposure": events.exposure}),
        d["citation"])
    return events
