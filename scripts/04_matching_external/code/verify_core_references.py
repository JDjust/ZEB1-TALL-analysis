#!/usr/bin/env python3
"""Fetch Crossref identity metadata for core cited papers; no full-text claim check."""
from pathlib import Path
import json
from datetime import datetime, timezone
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/core_reference_identity.json"
DOIS = [
    "10.1016/j.bbadis.2018.05.013",
    "10.1038/s41591-022-02112-7",
    "10.1038/s43018-024-00863-5",
    "10.1038/s41467-025-65134-y",
    "10.1038/s41467-021-25960-2",
    "10.1038/s41586-024-07807-0",
    "10.1038/s41598-023-39152-z",
    "10.1186/s13059-021-02540-7",
]


def main():
    session = requests.Session()
    rows = []
    for doi in DOIS:
        response = session.get(f"https://api.crossref.org/works/{doi}", timeout=30)
        response.raise_for_status()
        m = response.json()["message"]
        if m["DOI"].lower() != doi.lower() or not m.get("title"):
            raise ValueError(f"identity mismatch: {doi}")
        pub = m.get("published-print") or m.get("published") or m.get("published-online")
        year = pub["date-parts"][0][0] if pub else None
        authors = [a.get("family", "") for a in m.get("author", [])]
        rows.append({"doi": doi, "title": m["title"][0],
                     "authors_family": authors, "journal": m.get("container-title", [""])[0],
                     "year": year, "volume": m.get("volume"), "issue": m.get("issue"),
                     "page": m.get("page"), "article_number": m.get("article-number"),
                     "crossref_url": f"https://api.crossref.org/works/{doi}",
                     "verification_level": "bibliographic identity only; supporting passages not inspected"})
    payload = {"checked_utc": datetime.now(timezone.utc).isoformat(), "records": rows}
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"verified {len(rows)} DOI-title identities -> {OUT}")


if __name__ == "__main__":
    main()
