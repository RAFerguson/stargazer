"""Tests of the repeating-FRB burst-list connector (stargazer.sources.frb).

Unmarked tests run offline on trimmed copies of the real deposits in tests/fixtures/frb (see ATTRIBUTION.txt
there for origin and licences), fed to `fetch_bursts` through a pre-filled Cache so the full path runs without
the network. The interval-test check uses a synthetic table only; no search is run on real burst lists.
Tests marked `network` read Zenodo (deselected by default; run with `pytest -m network`).
"""
from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from stargazer.core import EventList
from stargazer.detect.intervals import clock_test, to_bits
from stargazer.sources import Cache, frb

FIXTURES = Path(__file__).parent / "fixtures" / "frb"
NRT_FILES = {"bursts": "nrt_bursts_excerpt.csv", "supplementary": "nrt_supplementary_excerpt.csv",
             "scans": "nrt_observation_times_excerpt.txt"}


def _seed(cache: Cache, name: str, files: dict[str, bytes]) -> None:
    d = frb.DATASETS[name]
    for key, data in files.items():
        cache.put(frb.zenodo_url(d["record"], d["files"][key]), data)


def _nrt_cache(tmp_path) -> Cache:
    cache = Cache(tmp_path / "cache")
    _seed(cache, "FRB 20220912A", {k: (FIXTURES / v).read_bytes() for k, v in NRT_FILES.items()})
    return cache


# ------------------------------------------------------------------------------------------------ offline


def test_names_and_listing():
    assert frb.canonical_name("FRB20220912A") == "FRB 20220912A"
    assert frb.canonical_name("20180916b") == "FRB 20180916B"
    with pytest.raises(KeyError):
        frb.canonical_name("FRB 20121102A")
    listing = {s["name"]: s for s in frb.list_sources()}
    assert set(listing) == {"FRB 20220912A", "FRB 20180916B", "FRB 20220529"}
    assert listing["FRB 20220912A"]["n_bursts"] == 696 and listing["FRB 20220912A"]["has_exposure"]
    assert all(s["urls"] and s["urls"][0].startswith("https://zenodo.org/api/records/") for s in listing.values())
    assert "%28online%29" in listing["FRB 20220912A"]["urls"][1]


def test_nrt_times_scans_and_toa_check(tmp_path):
    events = frb.fetch_bursts("FRB 20220912A", cache=_nrt_cache(tmp_path))
    assert isinstance(events, EventList) and len(events) == 25 and events.time_unit == "s"
    m = events.meta
    # B0 is the first burst: file start 59867.901782407411 + 17.0225 s; B1 follows 0.092704 s later.
    assert events.time[0] == 0.0 and m["burst_id"][:2] == ["B0", "B1"]
    assert events.time[1] == pytest.approx(17.115204 - 17.0225, abs=1e-6)
    assert float(m["epoch_mjd"]) == pytest.approx(59867.901782407411 + 17.0225 / 86400, abs=1e-11)
    assert m["time_frame"] == "topocentric" and "Heimdall" in m["dispersion_reference"]
    assert np.all(np.diff(events.time) >= 0) and len(m["snr"]) == len(m["dm"]) == len(events)
    assert events.amplitude[m["burst_id"].index("B10")] == pytest.approx(13.58)       # fluence, Jy ms
    # Exposure: the five FRB 20220912A scans of the excerpt (other sources ignored), every burst inside one.
    assert m["n_scans"] == 5 and len(events.exposure) == 5 and m["n_bursts_outside_exposure"] == 0
    assert events.exposure[0][0] < 0 < events.exposure[0][1]
    assert events.exposure[0][1] - events.exposure[0][0] == pytest.approx(1298.7)
    check = m["supplementary_toa_check"]
    assert sorted(check["supplementary_is_barycentric"]) == ["B318", "B38"]
    assert sorted(c["burst_id"] for c in check["conflict"]) == ["B384", "B411", "B50"]
    assert len(check["agrees_topocentric"]) == 20 and not check["missing"]
    prov = events.provenance
    assert prov.source == "frb_nrt_20220912a" and "Konijn" in prov.citation and len(prov.content_sha256) == 64
    assert set(prov.query["file_sha256"]) == set(NRT_FILES)


def test_nrt_barycentre_matches_paper_table(tmp_path):
    events = frb.fetch_bursts("FRB 20220912A", cache=_nrt_cache(tmp_path), barycentre=True)
    m = events.meta
    assert m["time_frame"] == "barycentric" and "TDB" in m["time_standard"]
    epoch = float(m["epoch_mjd"])
    # Konijn et al. Table 2: barycentric TDB at infinite frequency. Our times are not dedispersed, so they differ
    # by the dispersion delay at the Heimdall reference (~0.3-0.6 s at DM ~220 and 1.2-1.7 GHz).
    for burst, paper_mjd in (("B38", 59870.88721446), ("B318", 59879.83220702)):
        t = epoch + events.time[m["burst_id"].index(burst)] / 86400
        assert 0.0 < (t - paper_mjd) * 86400 < 1.0
    topo = frb.fetch_bursts("FRB 20220912A", cache=_nrt_cache(tmp_path / "b"))
    shift = (float(m["epoch_mjd"]) - float(topo.meta["epoch_mjd"])) * 86400
    assert 300 < shift < 500                                       # Roemer delay plus TT-UTC for this epoch
    assert events.time[-1] - topo.time[-1] == pytest.approx(0.0, abs=30.0)   # nine days later, delay changed
    assert m["n_bursts_outside_exposure"] == 0


def test_r3_compilation(tmp_path):
    cache = Cache(tmp_path / "cache")
    _seed(cache, "FRB 20180916B", {"package": (FIXTURES / "r3_package_excerpt.zip").read_bytes()})
    events = frb.fetch_bursts("FRB 20180916B", cache=cache)
    m = events.meta
    assert len(events) == 150 and m["n_bursts_by_instrument"]["CHIME"] == 55 and m["n_bursts_by_instrument"]["Apertif"] == 52
    assert m["duplicates_dropped"] == [{"instrument": "CHIME", "mjd": 58835.17323650127}]
    assert len(m["close_pairs_under_1s"]) == 2
    assert float(m["epoch_mjd"]) == pytest.approx(58377.427312524305)
    assert events.time[-1] / 86400 == pytest.approx(59143.81778929 - 58377.427312524305, abs=1e-9)
    gmrt = np.array(m["instrument"]) == "GMRT650"
    assert gmrt.sum() == 15 and np.all(np.isnan(events.amplitude[gmrt]))
    assert np.all(np.isfinite(events.amplitude[~gmrt]))
    assert events.exposure and all(b > a for a, b in events.exposure)
    assert "DSS-63-2250" in m["exposure_by_instrument"]                # windows of instruments with no bursts
    apertif = frb.fetch_bursts("FRB20180916B", cache=cache, instrument="Apertif")
    assert len(apertif) == 52 and set(apertif.meta["instrument"]) == {"Apertif"}
    assert apertif.meta["n_bursts_outside_exposure"] == 3
    assert apertif.provenance.query["instrument"] == "Apertif" and "Pastor-Marazuela" in apertif.provenance.citation
    with pytest.raises(KeyError):
        frb.fetch_bursts("FRB 20180916B", cache=cache, instrument="Arecibo")
    with pytest.raises(ValueError):
        frb.fetch_bursts("FRB 20220529", instrument="CHIME")


def test_fast_20220529_table(tmp_path):
    cache = Cache(tmp_path / "cache")
    _seed(cache, "FRB 20220529", {"table": (FIXTURES / "fast_frb20220529_table_s5_excerpt.csv").read_bytes()})
    events = frb.fetch_bursts("FRB 20220529", cache=cache)
    m = events.meta
    assert len(events) == 33 and events.exposure is None and m["time_frame"] == "barycentric"
    assert m["epoch_mjd"] == "59758.9282742" and m["same_mjd_rows"] == ["59808.8361080"]
    assert np.sum(np.diff(events.time) == 0) == 1
    assert events.amplitude[0] == 27 and len(m["rm"]) == len(events)


def test_malformed_inputs():
    with pytest.raises(ValueError):
        frb.parse_nrt_bursts("BurstID,Candidate name\nB0,nonsense.h5\n")
    with pytest.raises(ValueError):
        frb.parse_fast_20220529("59758.9,27\n")
    with pytest.raises(ValueError):
        frb.parse_r3_json({"bursts": {"X": [1.0]}, "snr": {"X": []}, "obs_startmjds": {}, "obs_duration": {}})
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("files/other.json", json.dumps({}))
    with pytest.raises(KeyError):
        frb.read_r3_zip(buf.getvalue())


def test_eventlist_feeds_interval_tests_synthetic(tmp_path):
    """Shape and units only: a synthetic clock (0.7 s, 40% of ticks empty) written in the FAST table's format."""
    rng = np.random.default_rng(7)
    ticks = np.flatnonzero(rng.random(200) < 0.6)
    mjd = 60000.0 + (ticks * 0.7 + rng.normal(0, 0.003, ticks.size)) / 86400
    header = "# synthetic\n#MJD,S/N,DM_det (pc cm^-3),f_low (MHz),f_high (MHz),DM_stru (pc cm^-3),RM (rad m^-2),f_L (%)\n"
    body = "".join(f"{x:.7f},20,250.0,1000,1500,250 ± 1,0 ± 1,90 ± 1\n" for x in mjd)
    cache = Cache(tmp_path / "cache")
    _seed(cache, "FRB 20220529", {"table": (header + body).encode()})
    events = frb.fetch_bursts("FRB 20220529", cache=cache)
    assert events.time_unit == "s" and events.time[0] == 0.0
    assert events.time[-1] == pytest.approx((ticks[-1] - ticks[0]) * 0.7, abs=0.05)
    result = clock_test(events, np.random.default_rng(1), n_surrogates=40, min_period=0.3, max_period=1.5)
    assert result.period == pytest.approx(0.7, rel=1e-3)
    bits = to_bits(events, result.period, result.phase)
    assert len(bits) >= ticks[-1] - ticks[0] and bits.bits.sum() == pytest.approx(ticks.size, abs=3)


# ------------------------------------------------------------------------------------------------ network


@pytest.mark.network
def test_live_nrt(tmp_path):
    events = frb.fetch_bursts("FRB 20220912A", cache=Cache(tmp_path / "cache"))
    m = events.meta
    assert len(events) == 696 and m["n_scans"] == 68 and m["exposure_hours"] == pytest.approx(60.99, abs=0.01)
    assert m["n_bursts_outside_exposure"] == 0 and len(m["supplementary_toa_check"]["conflict"]) == 6
    assert len(m["supplementary_toa_check"]["agrees_topocentric"]) == 681


@pytest.mark.network
def test_live_r3_and_fast(tmp_path):
    cache = Cache(tmp_path / "cache")
    r3 = frb.fetch_bursts("FRB 20180916B", cache=cache)
    assert len(r3) == 150 and r3.meta["exposure_by_instrument"]["CHIME"]["n_windows"] == 824
    fast = frb.fetch_bursts("FRB 20220529", cache=cache)
    assert len(fast) == 643 and fast.meta["epoch_mjd"].startswith("59758.")
