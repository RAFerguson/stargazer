#!/usr/bin/env python3
"""Check every arXiv identifier cited in the research notes against arXiv.

For each identifier the script asks the arXiv API whether the paper exists and records its title,
first author and date. It also counts the vocabulary shared between the paper's title and abstract
and the note line that cites it; a low count flags the citation for a person to read. A low count
is not proof of error, and a high count is not proof the claim is right.

The script is patient and resumable: it keeps what it has already resolved in a JSON file, asks
for the rest in small batches, and backs off when arXiv answers 429.

    python arxiv_citation_check.py ../research_notes/<title>/ resolved.json citations.tsv
"""
from __future__ import annotations

import csv
import json
import pathlib
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ATOM = "{http://www.w3.org/2005/Atom}"
ID_PATTERNS = (re.compile(r"arxiv\.org/(?:abs|pdf|html)/(\d{4}\.\d{4,5})", re.I), re.compile(r"arXiv:(\d{4}\.\d{4,5})"))
STOP = set("which their there these those about would could should between through other using based "
           "paper study results result search searches first found shows shown from with that this have".split())
BATCH, PAUSE_S, MAX_MINUTES = 20, 12, 26


def words(text: str) -> set[str]:
    return {w for w in re.findall(r"[a-z]{5,}", text.lower()) if w not in STOP}


def fetch(ids: list[str]) -> dict | None:
    """Return metadata for the identifiers, or None when the service refuses or the reply is unreadable."""
    url = "https://export.arxiv.org/api/query?max_results=%d&id_list=%s" % (len(ids), ",".join(ids))
    run = subprocess.run(["curl", "-sS", "-L", "--max-time", "60", "-A", "citation-check/0.1", "-w", "\n%{http_code}", url],
                         capture_output=True, text=True)
    body, _, code = run.stdout.rpartition("\n")
    if code.strip() != "200":
        return None
    try:
        root = ET.fromstring(body)
    except ET.ParseError:
        return None
    found = {}
    for entry in root.findall(ATOM + "entry"):
        ident = re.sub(r"v\d+$", "", (entry.findtext(ATOM + "id") or "").rsplit("/abs/", 1)[-1])
        title = " ".join((entry.findtext(ATOM + "title") or "").split())
        if not title or title.lower() == "error":
            continue
        author = entry.find(ATOM + "author")
        found[ident] = {"title": title, "published": (entry.findtext(ATOM + "published") or "")[:10],
                        "first_author": author.findtext(ATOM + "name") if author is not None else "",
                        "abstract": " ".join((entry.findtext(ATOM + "summary") or "").split())}
    return found


def main() -> None:
    notes_dir, resolved_path, table_path = (pathlib.Path(arg) for arg in sys.argv[1:4])
    citations = []
    for path in sorted(notes_dir.glob("*.md")):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            for ident in sorted({m for pattern in ID_PATTERNS for m in pattern.findall(line)}):
                citations.append({"note": path.name, "line": number, "id": ident, "text": line.strip()})
    ids = sorted({c["id"] for c in citations})
    resolved = json.loads(resolved_path.read_text()) if resolved_path.exists() else {}
    todo = [i for i in ids if i not in resolved]
    deadline = time.time() + MAX_MINUTES * 60
    for start in range(0, len(todo), BATCH):
        batch, wait = todo[start:start + BATCH], 90
        while time.time() < deadline:
            got = fetch(batch)
            if got is not None:
                resolved.update(got)
                resolved_path.write_text(json.dumps(resolved, indent=1, sort_keys=True))
                break
            time.sleep(wait)
            wait = min(wait * 2, 360)
        time.sleep(PAUSE_S)

    by_id: dict[str, dict] = {}
    for c in citations:
        row = by_id.setdefault(c["id"], {"notes": set(), "min_overlap": None})
        row["notes"].add(c["note"].removesuffix(".md"))
        meta = resolved.get(c["id"])
        if meta and meta.get("abstract"):
            overlap = len(words(c["text"]) & words(meta["title"] + " " + meta["abstract"]))
            row["min_overlap"] = overlap if row["min_overlap"] is None else min(row["min_overlap"], overlap)
    with table_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["arxiv_id", "status", "published", "first_author", "title", "min_shared_words", "cited_in"])
        for ident in ids:
            meta = resolved.get(ident)
            writer.writerow([ident, "resolved" if meta else "not checked", meta["published"] if meta else "",
                             meta["first_author"] if meta else "", meta["title"] if meta else "",
                             "" if by_id[ident]["min_overlap"] is None else by_id[ident]["min_overlap"],
                             ", ".join(sorted(by_id[ident]["notes"]))])
    done = sum(1 for i in ids if i in resolved)
    print(f"{len(citations)} citations of {len(ids)} distinct identifiers; {done} resolved; {len(ids) - done} not checked")


if __name__ == "__main__":
    main()
