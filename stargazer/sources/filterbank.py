"""SIGPROC filterbank (.fil) files: the simplest format radio SETI dynamic spectra come in.

A file is a keyword header followed by raw samples. The header starts with the string HEADER_START and ends
with HEADER_END; every string is written as a 4-byte little-endian length followed by its characters, and each
keyword is followed by its value: a 4-byte integer, an 8-byte double or another string, according to the
keyword. The data are time-major: for every time sample, for every IF (polarisation), nchans values of nbits
bits (8 and 16 bits unsigned integers, 32 bits IEEE floats). fch1 is the centre of the first channel in MHz,
foff the channel step in MHz (usually negative), tsamp the integration time in seconds and tstart the start
of the scan as a Modified Julian Date.

Only 8, 16 and 32-bit data are handled (1, 2 and 4-bit packed data are refused), and only little-endian files,
which is what SIGPROC, rawspec and blimpy write. The header keyword FREQUENCY_START (a channel-by-channel
frequency table) is not supported. Breakthrough Listen's HDF5 product is not handled here (it needs h5py).
"""
from __future__ import annotations

import os
import struct
from pathlib import Path
from typing import Any

import numpy as np

from ..detect.narrowband import Spectrogram

INT_KEYS = {"telescope_id", "machine_id", "data_type", "barycentric", "pulsarcentric", "nbits", "nsamples",
            "nchans", "nifs", "nbeams", "ibeam"}
DOUBLE_KEYS = {"tstart", "tsamp", "fch1", "foff", "refdm", "az_start", "za_start", "src_raj", "src_dej",
               "period"}
STRING_KEYS = {"source_name", "rawdatafile"}
CHAR_KEYS = {"signed"}
WRITE_ORDER = ["telescope_id", "machine_id", "data_type", "rawdatafile", "source_name", "barycentric",
               "pulsarcentric", "az_start", "za_start", "src_raj", "src_dej", "tstart", "tsamp", "nbits", "fch1",
               "foff", "nchans", "nifs", "nbeams", "ibeam", "refdm", "period"]


class FilterbankError(ValueError):
    """The file is not a filterbank file this reader can handle."""


def _read_string(fh) -> str:
    raw = fh.read(4)
    if len(raw) < 4:
        raise FilterbankError("file ended inside the header")
    (length,) = struct.unpack("<i", raw)
    if not 0 < length < 256:
        raise FilterbankError(f"implausible header string length {length}: not a little-endian SIGPROC file")
    text = fh.read(length)
    if len(text) < length:
        raise FilterbankError("file ended inside the header")
    return text.decode("ascii", errors="replace")


def read_header(path: str | os.PathLike) -> tuple[dict[str, Any], int]:
    """The header keywords of a filterbank file and the byte offset at which its data start."""
    with open(path, "rb") as fh:
        if _read_string(fh) != "HEADER_START":
            raise FilterbankError(f"{path}: does not start with HEADER_START")
        header: dict[str, Any] = {}
        while True:
            key = _read_string(fh)
            if key == "HEADER_END":
                break
            if key in INT_KEYS:
                header[key] = struct.unpack("<i", fh.read(4))[0]
            elif key in DOUBLE_KEYS:
                header[key] = struct.unpack("<d", fh.read(8))[0]
            elif key in STRING_KEYS:
                header[key] = _read_string(fh)
            elif key in CHAR_KEYS:
                header[key] = struct.unpack("<b", fh.read(1))[0]
            else:
                raise FilterbankError(f"{path}: unsupported header keyword {key!r}")
        offset = fh.tell()
    for key in ("nchans", "nbits", "tsamp", "fch1", "foff"):
        if key not in header:
            raise FilterbankError(f"{path}: header lacks {key}")
    header.setdefault("nifs", 1)
    return header, offset


def _dtype(nbits: int, signed: bool = False) -> np.dtype:
    if nbits == 8:
        return np.dtype("i1" if signed else "u1")
    if nbits == 16:
        return np.dtype("<i2" if signed else "<u2")
    if nbits == 32:
        return np.dtype("<f4")
    raise FilterbankError(f"nbits = {nbits} is not supported (8, 16 or 32 only)")


def _index_range(values: np.ndarray, lo: float, hi: float) -> tuple[int, int]:
    """First and one-past-last index of the values inside [lo, hi] (values monotonic)."""
    lo, hi = min(lo, hi), max(lo, hi)
    inside = np.flatnonzero((values >= lo) & (values <= hi))
    if inside.size == 0:
        raise ValueError(f"no samples between {lo} and {hi}")
    return int(inside[0]), int(inside[-1]) + 1


def read(path: str | os.PathLike, f_range: tuple[float, float] | None = None,
         t_range: tuple[float, float] | None = None, if_index: int = 0) -> Spectrogram:
    """Read a filterbank file, or a slice of it, as a Spectrogram.

    f_range is (low, high) in Hz: channels whose centre lies inside are kept. t_range is (start, stop) in
    seconds since the start of the scan: integrations that begin inside are kept. The data are memory-mapped,
    so only the slice asked for is read from disk. if_index picks one IF (polarisation) when nifs > 1.
    meta holds every header keyword, plus "file" and "channel_offset"/"sample_offset" (where the slice starts).
    The format does not record how many raw samples were averaged per cell, so meta has no "n_avg" and the
    search will estimate it unless the caller adds it.
    """
    header, offset = read_header(path)
    nchans, nifs = int(header["nchans"]), int(header["nifs"])
    dtype = _dtype(int(header["nbits"]), bool(header.get("signed", 0)))
    data_bytes = os.path.getsize(path) - offset
    nsamples = data_bytes // (nchans * nifs * dtype.itemsize)
    if header.get("nsamples") and header["nsamples"] < nsamples:
        nsamples = int(header["nsamples"])
    if nsamples < 1:
        raise FilterbankError(f"{path}: no data after the header")
    if not 0 <= if_index < nifs:
        raise ValueError(f"if_index {if_index} out of range for nifs = {nifs}")
    data = np.memmap(path, dtype=dtype, mode="r", offset=offset, shape=(nsamples, nifs, nchans))
    tsamp = float(header["tsamp"])
    fch1, foff = float(header["fch1"]) * 1e6, float(header["foff"]) * 1e6
    c0, c1 = (0, nchans) if f_range is None else _index_range(fch1 + foff * np.arange(nchans), *f_range)
    t0, t1 = (0, nsamples) if t_range is None else _index_range(tsamp * np.arange(nsamples), *t_range)
    power = np.array(data[t0:t1, if_index, c0:c1], dtype=float)
    del data
    meta: dict[str, Any] = dict(header)
    meta["file"] = str(Path(path))
    meta["channel_offset"], meta["sample_offset"] = c0, t0
    if "tstart" in header:
        meta["tstart"] = float(header["tstart"]) + t0 * tsamp / 86400.0
    return Spectrogram(power, tsamp, fch1 + c0 * foff, foff, str(header.get("source_name", "")), meta)


def _pack(key: str, value: Any) -> bytes:
    name = struct.pack("<i", len(key)) + key.encode("ascii")
    if key in INT_KEYS:
        return name + struct.pack("<i", int(value))
    if key in DOUBLE_KEYS:
        return name + struct.pack("<d", float(value))
    if key in STRING_KEYS:
        text = str(value).encode("ascii", errors="replace")
        return name + struct.pack("<i", len(text)) + text
    if key in CHAR_KEYS:
        return name + struct.pack("<b", int(value))
    raise FilterbankError(f"cannot write header keyword {key!r}")


def write(spec: Spectrogram, path: str | os.PathLike, nbits: int = 32) -> Path:
    """Write a Spectrogram as a single-IF filterbank file.

    Header keywords in spec.meta are kept where SIGPROC knows them (telescope_id, tstart, src_raj and so on);
    nchans, nbits, tsamp, fch1, foff and nifs always come from the spectrogram. With 8 or 16 bits every value
    must already be an integer within the unsigned range of that width: values are not rescaled silently,
    because the format has no keyword to record a scale factor. 32-bit data are written as float32, which
    rounds float64 power to about 7 significant digits.
    """
    dtype = _dtype(nbits)
    power = spec.power
    if nbits in (8, 16):
        top = 2 ** nbits - 1
        if np.any(power != np.round(power)) or power.min() < 0 or power.max() > top:
            raise ValueError(f"{nbits}-bit output needs integer values in [0, {top}]; rescale and round first")
    header: dict[str, Any] = {"telescope_id": 0, "machine_id": 0, "data_type": 1, "source_name": spec.source_id,
                              "tstart": float(spec.meta.get("tstart", 0.0))}
    for key in WRITE_ORDER:
        if key in spec.meta and key not in ("nbits", "nchans", "nifs", "tsamp", "fch1", "foff", "nsamples"):
            header[key] = spec.meta[key]
    header.update({"tsamp": spec.t_samp, "nbits": nbits, "fch1": spec.f_start / 1e6, "foff": spec.f_step / 1e6,
                   "nchans": spec.n_freq, "nifs": 1})
    if not header.get("source_name"):
        header["source_name"] = spec.source_id or "unknown"
    out = Path(path)
    with open(out, "wb") as fh:
        fh.write(struct.pack("<i", 12) + b"HEADER_START")
        for key in WRITE_ORDER:
            if key in header:
                fh.write(_pack(key, header[key]))
        fh.write(struct.pack("<i", 10) + b"HEADER_END")
        fh.write(np.ascontiguousarray(power, dtype=dtype).tobytes())
    return out
