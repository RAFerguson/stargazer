"""Gamma-ray light curves from the Fermi-LAT Light Curve Repository (LCR) at the Fermi Science Support Center.

The LCR's JSON endpoint is undocumented; this module was written against live responses of 22 Nov 2026.

Source names. The repository holds the ~1,525 variable sources of the 4FGL-DR2 catalogue and names each by its
4FGL catalogue designation with the space, for example "4FGL J1555.7+1111" (the blazar PG 1553+113). Any other
name, including the common name, or "4FGL_J1555.7+1111", returns HTTP 200 with every array empty, which this
module reports as a LookupError.

Response. A JSON object of parallel arrays keyed by Fermi mission elapsed time (MET, seconds since
2001-01-01T00:00:00 UTC): "ts" ([MET, TS]) and "fit_convergence" ([MET, code], 0 = converged) and "fit_tolerance"
for every bin; "flux" ([MET, photon flux]) and "flux_error" ([MET, lower, upper], the 1-sigma interval of the
flux, not offsets) for bins with TS >= ts_min; "flux_upper_limits" ([MET, 95% upper limit]) for the others;
plus photon_index, photon_index_interval, dlogl, EG, GAL and bin_id.

Time stamps are bin centres: weekly stamps are 604800 s apart starting at MET 239889601, which is half a week
after the repository's start of MET 239587201 (2008-08-05); the monthly series starts at 239587201 + 15 days
and the 3-day ("daily") one at 239587201 + 1.5 days.
"""
from __future__ import annotations

import json
import logging
import urllib.parse
from typing import Any

import numpy as np

from ..core import Provenance, Series, content_hash
from .cache import Cache, fetch_bytes, utc_now
from .registry import FERMI_CITATION

log = logging.getLogger(__name__)

SOURCE_ID = "fermi_lcr"
ENDPOINT = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/queryDB.php"
MJDREF_TT = 51910.0 + 7.428703703703703e-4     # Fermi MET zero, 2001-01-01T00:00:00 UTC, as an MJD in TT
TIME_STANDARD = "MJD_TT (Fermi MET / 86400 + 51910.000742870); bin centres"
CADENCES = {"daily": "daily", "3day": "daily", "weekly": "weekly", "monthly": "monthly"}
BIN_DAYS = {"daily": 3.0, "weekly": 7.0, "monthly": 30.0}


def met_to_mjd(met: np.ndarray | float) -> np.ndarray:
    """Fermi mission elapsed time (s) to MJD in the TT scale."""
    return MJDREF_TT + np.asarray(met, dtype=float) / 86400.0


def lcr_url(source_name: str, cadence: str = "weekly") -> str:
    """The LCR request for a source's photon-flux, fixed-index light curve with TS >= 4 detections."""
    if cadence not in CADENCES:
        raise ValueError(f"cadence must be one of {sorted(CADENCES)}, not {cadence!r}")
    params = {"typeOfRequest": "lightCurveData", "source_name": source_name, "cadence": CADENCES[cadence],
              "flux_type": "photon", "index_type": "fixed", "ts_min": "4"}
    return f"{ENDPOINT}?{urllib.parse.urlencode(params)}"


def _by_met(pairs: list[list[float]]) -> dict[float, list[float]]:
    return {float(p[0]): [float(x) for x in p[1:]] for p in pairs}


def parse_lcr_json(payload: dict[str, Any]) -> dict[str, Any]:
    """Detections from an LCR lightCurveData response.

    Returns a dict: met, flux, error (arrays; flux and error in photons cm^-2 s^-1, 100 MeV-100 GeV; error is
    half the width of the reported 1-sigma interval), n_bins, n_upper_limits (bins with only an upper limit),
    n_not_converged (detection bins whose fit_convergence is non-zero), all dropped from the arrays.
    """
    for key in ("ts", "flux", "flux_error", "flux_upper_limits", "fit_convergence"):
        if key not in payload:
            raise ValueError(f"LCR response lacks {key!r} (has {sorted(payload)})")
    flux = _by_met(payload["flux"])
    interval = _by_met(payload["flux_error"])
    limits = _by_met(payload["flux_upper_limits"])
    convergence = {met: int(v[0]) for met, v in _by_met(payload["fit_convergence"]).items()}
    bins = sorted(_by_met(payload["ts"]))
    met, value, error = [], [], []
    n_limits = n_bad = 0
    for t in bins:
        if t in limits or t not in flux:
            n_limits += 1
            continue
        if convergence.get(t, 0) != 0:
            n_bad += 1
            continue
        f = flux[t][0]
        if t not in interval:
            raise ValueError(f"LCR response has a flux but no flux_error at MET {t}")
        lower, upper = interval[t][:2]
        if not lower <= f <= upper:
            raise ValueError(f"LCR flux_error at MET {t} is not an interval around the flux: {lower}, {f}, {upper}")
        met.append(t)
        value.append(f)
        error.append((upper - lower) / 2.0)
    return {"met": np.array(met), "flux": np.array(value), "error": np.array(error), "n_bins": len(bins),
            "n_upper_limits": n_limits, "n_not_converged": n_bad}


def fetch_lightcurve(source_name: str, *, cadence: str = "weekly", cache: Cache | None = None) -> Series:
    """A Fermi-LAT photon-flux light curve from the LCR (see the module docstring for names and time stamps).

    cadence is "daily" (the repository's 3-day bins; "3day" is accepted too), "weekly" or "monthly". Returns a
    Series of kind "flux" in photons cm^-2 s^-1 (100 MeV-100 GeV, fixed photon index), time in MJD (TT) at bin
    centres, with upper-limit bins and unconverged fits dropped and counted in meta. With a cache the first
    response is kept for good, although the repository itself is updated continually.
    """
    url = lcr_url(source_name, cadence)
    data, _ = fetch_bytes(url, cache=cache, timeout=120.0)
    try:
        payload = json.loads(data)
    except json.JSONDecodeError as error:
        raise ValueError(f"LCR returned non-JSON from {url}: {data[:200]!r}") from error
    parsed = parse_lcr_json(payload)
    if parsed["n_bins"] == 0:
        raise LookupError(f"the Fermi LCR has no light curve named {source_name!r} (use the 4FGL designation, "
                          f"e.g. '4FGL J1555.7+1111'); request was {url}")
    meta = {"source_name": source_name, "cadence": CADENCES[cadence], "bin_days": BIN_DAYS[CADENCES[cadence]],
            "n_bins": parsed["n_bins"], "n_upper_limits_dropped": parsed["n_upper_limits"],
            "n_not_converged_dropped": parsed["n_not_converged"], "flux_unit": "ph cm^-2 s^-1",
            "energy_range_mev": [100, 100000], "index_type": "fixed", "ts_min": 4}
    series = Series(met_to_mjd(parsed["met"]), parsed["flux"], parsed["error"], band="LAT 0.1-100 GeV",
                    source_id=f"4fgl:{source_name.removeprefix('4FGL ').strip()}", kind="flux",
                    time_standard=TIME_STANDARD, meta=meta)
    series.provenance = Provenance(SOURCE_ID, {"source_name": source_name, "cadence": cadence, "url": url},
                                   utc_now(), content_hash({"time": series.time, "value": series.value,
                                                            "error": series.error}), FERMI_CITATION)
    return series
