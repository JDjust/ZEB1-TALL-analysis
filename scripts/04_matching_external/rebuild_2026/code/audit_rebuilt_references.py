"""Audit rebuilt numeric citations against DOI records without changing the manuscript."""
from __future__ import annotations

import concurrent.futures
import difflib
import html
import re
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "manuscript/ZEB1_ZEB2_TALL_MASTER_REBUILT.md"
OUT = ROOT / "manuscript/REFERENCE_DOI_AUDIT.tsv"


def clean(value: str) -> str:
    value = re.sub(r"<[^>]+>", " ", html.unescape(str(value)))
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def audit(item):
    number, line, cited = item
    doi_match = re.search(r"doi:(\S+)", line, re.I)
    doi = doi_match.group(1).rstrip(". ") if doi_match else ""
    title_match = re.search(r"\. ([^.]+)\. \*[^*]+\*\.", line)
    title = title_match.group(1) if title_match else ""
    row = dict(old_number=number, cited=cited, doi=doi, manuscript_title=title,
               doi_title="", title_similarity=None, doi_volume="", doi_issue="",
               doi_page="", doi_online_year="", lookup_status="")
    if not doi:
        row["lookup_status"] = "NO_DOI"
        return row
    try:
        resp = requests.get("https://api.crossref.org/works/" + doi,
                            headers={"User-Agent": "ZEB1-TALL-submission-audit/1.0"}, timeout=20)
        resp.raise_for_status()
        data = resp.json()["message"]
        row["doi_title"] = html.unescape(data.get("title", [""])[0])
        row["title_similarity"] = round(difflib.SequenceMatcher(
            None, clean(title), clean(row["doi_title"])).ratio(), 3)
        row["doi_volume"] = data.get("volume", "")
        row["doi_issue"] = data.get("issue", "")
        row["doi_page"] = data.get("page", "")
        parts = data.get("published-online", data.get("published", {})).get("date-parts", [[]])
        row["doi_online_year"] = parts[0][0] if parts and parts[0] else ""
        row["lookup_status"] = "MATCH" if row["title_similarity"] >= .85 else "REVIEW_TITLE"
    except Exception as exc:
        row["lookup_status"] = "LOOKUP_ERROR:" + type(exc).__name__
    return row


def main():
    body, refs = MAIN.read_text(encoding="utf-8").split("## References", 1)
    cited = set()
    for group in re.findall(r"\[([\d,\s–-]+)\]", body):
        for start, end in re.findall(r"(\d+)(?:[–-](\d+))?", group):
            cited.update(range(int(start), int(end) + 1)) if end else cited.add(int(start))
    entries = []
    for line in refs.splitlines():
        match = re.match(r"^(\d+)\. (.+)$", line)
        if match:
            number = int(match.group(1))
            entries.append((number, match.group(2), number in cited))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(audit, entries))
    table = pd.DataFrame(rows).sort_values("old_number")
    table.to_csv(OUT, sep="\t", index=False)
    print("references", len(table), "cited", int(table.cited.sum()))
    print(table.lookup_status.value_counts().to_string())
    print("audit", OUT)


if __name__ == "__main__":
    main()
