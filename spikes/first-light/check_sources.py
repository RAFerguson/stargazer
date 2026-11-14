#!/usr/bin/env python3
"""Probe each data source once and report whether this machine can reach it.

Standard library only, so it runs before anything is installed. A source counts as
reachable when its server answers at all (any HTTP status), as blocked when the outbound
proxy refuses the connection, and as unreachable on a timeout or network error.

    python check_sources.py              # table to stdout
    python check_sources.py --markdown   # Markdown table
"""
from __future__ import annotations

import argparse
import concurrent.futures as futures
import pathlib
import urllib.error
import urllib.request

HERE = pathlib.Path(__file__).parent


def load_sources(path: pathlib.Path) -> list[tuple[str, str, str]]:
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        host, what, route = (part.strip() for part in line.split("|"))
        rows.append((host, what, route))
    return rows


def probe(host: str, timeout: float) -> tuple[str, str]:
    request = urllib.request.Request(f"https://{host}/", method="GET", headers={"User-Agent": "source-check/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return "reachable", f"HTTP {response.status}"
    except urllib.error.HTTPError as error:  # the server itself answered
        return "reachable", f"HTTP {error.code}"
    except urllib.error.URLError as error:
        reason = str(error.reason)
        if "Tunnel connection failed" in reason or "407" in reason:
            return "blocked", reason
        return "unreachable", reason
    except Exception as error:  # timeouts, TLS failures
        return "unreachable", f"{type(error).__name__}: {error}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sources", type=pathlib.Path, default=HERE / "sources.txt")
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()

    sources = load_sources(args.sources)
    with futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda row: probe(row[0], args.timeout), sources))

    if args.markdown:
        print("| Host | Serves | Route | Status | Detail |")
        print("|---|---|---|---|---|")
        for (host, what, route), (status, detail) in zip(sources, results):
            print(f"| `{host}` | {what} | {route} | {status} | {detail} |")
    else:
        width = max(len(host) for host, _, _ in sources)
        for (host, what, route), (status, detail) in zip(sources, results):
            print(f"{status:<12}{route:<7}{host:<{width + 2}}{detail}")
    counts = {state: sum(1 for status, _ in results if status == state) for state in ("reachable", "blocked", "unreachable")}
    print()
    print(", ".join(f"{count} {state}" for state, count in counts.items()))


if __name__ == "__main__":
    main()
