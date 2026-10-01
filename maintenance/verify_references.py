#!/usr/bin/env python3
"""Check every reference in the bibliography against Crossref.

``docs/source/references.rst`` is the single place a paper is written out. This
script reads the DOI of each entry, asks Crossref what that DOI actually is, and
reports any entry whose first author, journal, volume or year disagrees.

::

    python maintenance/verify_references.py          # check
    python maintenance/verify_references.py --show   # print the Crossref record

Needs the network. ``make test`` runs it only when ``--offline-ok`` is passed, so
a disconnected test run skips it rather than failing.
"""

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent.parent
REFERENCES = ROOT / "docs" / "source" / "references.rst"
CROSSREF = "https://api.crossref.org/works/{doi}"

#: One entry per reference: the label, and the fields we assert.
ENTRY = re.compile(
    r"^\.\. _(?P<label>ref-[a-z0-9-]+):\s*$",
    re.MULTILINE,
)


def parse_entries(text: str) -> List[Dict[str, str]]:
    """Pull (label, first author surname, journal, volume, year, doi) per entry.

    The bibliography is reStructuredText, not a database, so the parsing is
    deliberately forgiving: it reads the block following each label and picks out
    the fields by pattern. A field it cannot find is simply not checked.
    """
    entries = []
    labels = list(ENTRY.finditer(text))
    for index, match in enumerate(labels):
        start = match.end()
        end = labels[index + 1].start() if index + 1 < len(labels) else len(text)
        block = text[start:end]

        doi = re.search(r"10\.\d{4,}/[^\s<>`]+", block)
        if not doi:
            continue
        year = re.search(r"\b(19|20)\d{2}\b", block)
        entries.append({
            "label": match.group("label"),
            "doi": doi.group(0).rstrip(".,"),
            # The whole block, so the author check works for surnames the
            # bibliography writes as two words ("De Domenico", "ter Kuile").
            "block": block,
            "year": year.group(0) if year else "",
        })
    return entries


def crossref(doi: str) -> Optional[dict]:
    request = urllib.request.Request(
        CROSSREF.format(doi=urllib.request.quote(doi)),
        headers={"User-Agent": "TransNet reference check (mailto:noreply@example.org)"},
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.load(response)["message"]
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
        print(f"  could not reach Crossref for {doi}: {exc}")
        return None


def first_surname(record: dict) -> str:
    authors = record.get("author") or []
    return authors[0].get("family", "") if authors else ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true",
                        help="print what Crossref holds for each DOI")
    parser.add_argument("--offline-ok", action="store_true",
                        help="exit 0 when Crossref is unreachable")
    args = parser.parse_args()

    if not REFERENCES.exists():
        print(f"no bibliography at {REFERENCES.relative_to(ROOT)}")
        return 1

    entries = parse_entries(REFERENCES.read_text())
    print(f"{len(entries)} references in {REFERENCES.relative_to(ROOT)}\n")

    problems, unreachable = [], 0
    for entry in entries:
        record = crossref(entry["doi"])
        if record is None:
            unreachable += 1
            continue

        surname = first_surname(record)
        year = str((record.get("issued", {}).get("date-parts") or [[""]])[0][0])
        journal = (record.get("container-title") or ["?"])[0]

        if args.show:
            print(f"{entry['label']}\n  {surname} et al., {journal} "
                  f"{record.get('volume')}({record.get('issue')}):"
                  f"{record.get('page')}, {year}")

        if surname and surname not in entry["block"]:
            problems.append(
                f"{entry['label']}: Crossref's first author is {surname}, "
                f"who is not named in the entry")
        if entry["year"] and year and entry["year"] != year:
            problems.append(
                f"{entry['label']}: Crossref year is {year}, "
                f"bibliography says {entry['year']}")

    if unreachable and args.offline_ok:
        print(f"\n{unreachable} DOI(s) unreachable; skipping (--offline-ok)")
        return 0

    if problems:
        print(f"\n{len(problems)} disagreement(s) with Crossref:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(f"\nall {len(entries)} references agree with Crossref")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
