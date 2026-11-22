"""A content-addressed download cache and the one HTTP helper every connector uses.

Objects are stored under their SHA-256 (`objects/ab/abcdef...`); a key (normally a URL plus byte range) maps to
an object through a small JSON sidecar (`keys/<sha256 of key>.json`) recording key, sha256, size and UTC time.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import logging
import os
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

from .registry import USER_AGENT

log = logging.getLogger(__name__)


class FetchError(RuntimeError):
    """A download failed; the message names the URL."""


def utc_now() -> str:
    """The current UTC time as an ISO 8601 string with a trailing Z, to the second."""
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".tmp-")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


class Cache:
    """Content-addressed byte store on local disk.

    root defaults to $STARGAZER_CACHE, else ~/.stargazer/cache. Safe to share between processes: every write is
    atomic, and identical content is stored once.
    """

    def __init__(self, root: str | Path | None = None) -> None:
        if root is None:
            root = os.environ.get("STARGAZER_CACHE") or Path.home() / ".stargazer" / "cache"
        self.root = Path(root).expanduser()

    def _sidecar(self, key: str) -> Path:
        return self.root / "keys" / f"{hashlib.sha256(key.encode()).hexdigest()}.json"

    def _object(self, sha256: str) -> Path:
        return self.root / "objects" / sha256[:2] / sha256

    def put(self, key: str, data: bytes) -> str:
        """Store data under key; returns the SHA-256 (hex) of the data."""
        sha = hashlib.sha256(data).hexdigest()
        target = self._object(sha)
        if not target.exists():
            _atomic_write(target, data)
        record = {"key": key, "sha256": sha, "size": len(data), "stored_utc": utc_now()}
        _atomic_write(self._sidecar(key), json.dumps(record, sort_keys=True).encode())
        return sha

    def info(self, key: str) -> dict | None:
        """The sidecar record for key (key, sha256, size, stored_utc), or None if key is not cached."""
        path = self._sidecar(key)
        if not path.exists():
            return None
        return json.loads(path.read_text())

    def get(self, key: str) -> bytes | None:
        """The bytes stored under key, or None if absent. Content that fails its hash check counts as absent."""
        record = self.info(key)
        if record is None:
            return None
        path = self._object(record["sha256"])
        if not path.exists():
            return None
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            log.warning("cache object for %s is corrupt; ignoring it", key)
            return None
        return data

    def has(self, key: str) -> bool:
        """Whether key has a stored, intact object."""
        record = self.info(key)
        return record is not None and self._object(record["sha256"]).exists()


def cache_key(url: str, byte_range: tuple[int, int] | None = None) -> str:
    """The cache key for a URL and an optional inclusive byte range."""
    return url if byte_range is None else f"{url}#bytes={byte_range[0]}-{byte_range[1]}"


def fetch_bytes(url: str, *, headers: dict[str, str] | None = None, byte_range: tuple[int, int] | None = None,
                cache: Cache | None = None, timeout: float = 60.0) -> tuple[bytes, str]:
    """Download url (optionally only the inclusive byte range (first, last)); returns (data, sha256 hex).

    Uses urllib, which honours the HTTPS_PROXY/HTTP_PROXY environment. With a cache, a hit is returned without
    any network access and a miss is stored. Raises FetchError naming the URL on any HTTP or network failure.
    """
    key = cache_key(url, byte_range)
    if cache is not None:
        data = cache.get(key)
        if data is not None:
            log.debug("cache hit %s", key)
            return data, hashlib.sha256(data).hexdigest()
    request_headers = {"User-Agent": USER_AGENT, **(headers or {})}
    if byte_range is not None:
        first, last = byte_range
        if first < 0 or last < first:
            raise ValueError(f"bad byte range {byte_range}")
        request_headers["Range"] = f"bytes={first}-{last}"
    request = urllib.request.Request(url, headers=request_headers)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            status = response.status
            data = response.read()
    except urllib.error.HTTPError as error:
        raise FetchError(f"HTTP {error.code} {error.reason} from {url}") from error
    except (urllib.error.URLError, OSError) as error:
        raise FetchError(f"could not fetch {url}: {getattr(error, 'reason', error)}") from error
    if byte_range is not None and status != 206:
        raise FetchError(f"{url} ignored the byte range {byte_range} (HTTP {status})")
    log.debug("fetched %d bytes from %s", len(data), key)
    if cache is not None:
        return data, cache.put(key, data)
    return data, hashlib.sha256(data).hexdigest()
