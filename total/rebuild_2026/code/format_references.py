"""Render the first 60 DOI-identity-checked references for the manuscript."""
from pathlib import Path
import json, re

base = Path(__file__).resolve().parents[1]
records = json.loads((base / "data/manifests/references_verified.json").read_text(encoding="utf-8"))[:60]
if any("error" in r for r in records):
    raise SystemExit("Unverified reference cannot be formatted")

def clean(x):
    return re.sub(r"\s+", " ", str(x or "")).strip()

lines = ["# References", "", "Bibliographic DOI/title identity was checked against Crossref on 2026-09-25. Citation-to-claim suitability remains a manuscript-level judgment.", ""]
for i, r in enumerate(records, 1):
    authors = r["authors"]
    displayed = [clean((a.get("family", "") + " " + a.get("given", "")).strip()) for a in authors[:6]]
    auth = ", ".join(displayed) + (", et al." if len(authors) > 6 else ".")
    vol = clean(r.get("volume"))
    issue = clean(r.get("issue"))
    page = clean(r.get("page"))
    detail = (vol + ("(" + issue + ")" if issue else "") + (":" + page if page else "")).strip(":")
    lines.append(f"{i}. {auth} {clean(r['title'])}. *{clean(r['journal'])}*. {r['year']};{detail}. doi:{r['doi']}")
    lines.append("")
(base / "manuscript/References.md").write_text("\n".join(lines), encoding="utf-8")
print("Rendered", len(records), "references")
