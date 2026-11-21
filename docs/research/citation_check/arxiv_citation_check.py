#!/usr/bin/env python3
"""Check every arXiv identifier cited in the research notes against arXiv.

For each identifier the script confirms that the paper exists and records its title, first author
and date. It also counts the vocabulary shared between the paper's title and abstract and each
note line that cites it, and keeps the smallest count; a low count flags the citation for a person
to read. A low count is not proof of error, and a high count is not proof the claim is right.

Two routes are available. The default asks the arXiv API in batches of 20 and backs off when it
answers 429. With --pages the script reads one abstract page every 15 seconds, the delay arXiv's
robots.txt asks for; use it when the API is rate-limiting. Either way the script is resumable:
what has been resolved is kept in a JSON file and never fetched twice.

    python arxiv_citation_check.py ../research_notes/<title>/ resolved.json citations.tsv [--pages]
"""
from __future__ import annotations

import csv
import html
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
AGENT = "citation-check/0.1"
MAX_MINUTES = 80


def words(text: str) -> set[str]:
    return {w for w in re.findall(r"[a-z]{5,}", text.lower()) if w not in STOP}


def curl(url: str) -> tuple[str, str]:
    run = subprocess.run(["curl", "-sS", "-L", "--max-time", "60", "-A", AGENT, "-w", "\n%{http_code}", url],
                         capture_output=True, text=True)
    body, _, code = run.stdout.rpartition("\n")
    return body, code.strip()


def from_api(ids: list[str]) -> dict | None:
    """Metadata for a batch, or None when the service refuses or the reply is unreadable."""
    body, code = curl("https://export.arxiv.org/api/query?max_results=%d&id_list=%s" % (len(ids), ",".join(ids)))
    if code != "200":
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


def from_page(ident: str) -> dict | str:
    """Metadata from one abstract page, or the HTTP status when the page could not be read."""
    page, code = curl(f"https://arxiv.org/abs/{ident}")
    tag = lambda name: [html.unescape(" ".join(m.split()))
                        for m in re.findall(r'<meta name="%s" content="(.*?)"\s*/?>' % name, page, re.S)]
    if code != "200" or not tag("citation_title"):
        return code
    return {"title": tag("citation_title")[0], "first_author": (tag("citation_author") or [""])[0],
            "published": (tag("citation_date") or [""])[0].replace("/", "-"), "abstract": (tag("citation_abstract") or [""])[0]}


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    use_pages = "--pages" in sys.argv
    notes_dir, resolved_path, table_path = (pathlib.Path(a) for a in args[:3])
    lines_citing: dict[str, list[str]] = {}
    cited_in: dict[str, set[str]] = {}
    for path in sorted(notes_dir.glob("*.md")):
        for line in path.read_text().splitlines():
            for ident in {m for pattern in ID_PATTERNS for m in pattern.findall(line)}:
                lines_citing.setdefault(ident, []).append(line)
                cited_in.setdefault(ident, set()).add(path.stem)
    ids = sorted(lines_citing)
    resolved = json.loads(resolved_path.read_text()) if resolved_path.exists() else {}

    def keep(ident: str, meta: dict) -> None:
        shared = min(len(words(line) & words(meta["title"] + " " + meta["abstract"])) for line in lines_citing[ident])
        resolved[ident] = {"title": meta["title"], "first_author": meta["first_author"], "published": meta["published"],
                           "min_shared_words": shared}
        resolved_path.write_text(json.dumps(resolved, indent=1, sort_keys=True))

    todo = [i for i in ids if i not in resolved]
    deadline = time.time() + MAX_MINUTES * 60
    if use_pages:
        for ident in todo:
            if time.time() > deadline:
                break
            got = from_page(ident)
            if got in ("403", "429"):
                print(f"stopped: arXiv answered {got}")
                break
            if isinstance(got, dict):
                keep(ident, got)
            time.sleep(15)
    else:
        for start in range(0, len(todo), 20):
            batch, wait = todo[start:start + 20], 90
            while time.time() < deadline:
                got = from_api(batch)
                if got is not None:
                    for ident, meta in got.items():
                        if ident in lines_citing:
                            keep(ident, meta)
                    break
                time.sleep(wait)
                wait = min(wait * 2, 360)
            time.sleep(12)

    with table_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["arxiv_id", "status", "published", "first_author", "title", "min_shared_words", "cited_in"])
        for ident in ids:
            meta = resolved.get(ident)
            writer.writerow([ident, "resolved" if meta else "not checked", meta["published"] if meta else "",
                             meta["first_author"] if meta else "", meta["title"] if meta else "",
                             meta["min_shared_words"] if meta else "", ", ".join(sorted(cited_in[ident]))])
    done = sum(1 for i in ids if i in resolved)
    print(f"{len(ids)} distinct identifiers; {done} resolved; {len(ids) - done} not checked")


if __name__ == "__main__":
    main()
