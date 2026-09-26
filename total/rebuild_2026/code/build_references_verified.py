"""Identity-check a manuscript candidate bibliography against Crossref DOI records.

Scientific support for each manuscript assertion still requires author reading;
this script verifies bibliographic identity and records the returned metadata.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib, json, re, time
import requests

BASE = Path(__file__).resolve().parents[1]
DEST = BASE / "data" / "manifests"
DEST.mkdir(parents=True, exist_ok=True)
DOIS = [
"10.1038/s41586-024-07807-0",
"10.1038/ng.3909",
"10.1126/science.aay3224",
"10.1126/sciimmunol.ade0182",
"10.1084/jem.20192360",
"10.1038/s43018-024-00863-5",
"10.1038/s41467-021-24044-5",
"10.1182/blood.2020010510",
"10.1016/j.molcel.2024.11.040",
"10.1038/s41598-023-39152-z",
"10.1172/jci61269",
"10.1038/nature10725",
"10.1016/S1470-2045(08)70314-0",
"10.1038/s41586-018-0436-0",
"10.1182/blood-2015-08-661702",
"10.1182/blood-2010-11-318873",
"10.1158/2159-8290.CD-21-0145",
"10.1126/science.1188989",
"10.1126/science.1188995",
"10.1038/s41590-018-0238-4",
"10.1016/j.ccr.2011.02.008",
"10.1016/j.ccr.2012.06.007",
"10.1038/leu.2016.82",
"10.1016/j.bbadis.2018.05.013",
"10.1182/blood-2012-09-458570",
"10.1038/s41388-022-02414-7",
"10.1038/leu.2015.162",
"10.1182/blood-2012-08-447839",
"10.1084/jem.20212383",
"10.1038/s41586-024-07944-6",
"10.1016/j.cels.2019.09.008",
"10.1084/jem.20150194",
"10.1182/blood-2009-08-231217",
"10.1038/s41467-021-25960-2",
"10.1038/s41467-020-19894-4",
"10.1038/s41467-022-35519-4",
"10.1093/bioinformatics/btab337",
"10.1038/nmeth.3999",
"10.1038/s41467-019-11950-y",
"10.1038/nmeth.4583",
"10.1038/s41588-020-0602-9",
"10.1200/JCO.20.00256",
"10.1200/JCO.21.02678",
"10.1038/s41375-018-0307-6",
"10.1200/JCO.2017.74.0449",
"10.1038/s41467-021-21346-6",
"10.1186/gb-2010-11-3-r25",
"10.1093/bioinformatics/btp616",
"10.1186/gb-2014-15-2-r29",
"10.1093/biomet/80.1.27",
"10.1016/0197-2456(86)90046-2",
"10.1002/sim.791",
"10.1038/nri1883",
"10.1038/ni.3514",
"10.1038/nature10279",
"10.1038/ncomms11171",
"10.4049/jimmunol.1301663",
"10.1182/blood-2012-11-465138",
"10.1182/blood-2013-03-491092",
"10.1101/gad.1897910",
"10.1038/nm.3665",
"10.1158/2159-8290.CD-14-0353",
"10.1038/nm.2651",
"10.1038/ng.2508",
"10.1038/s41467-025-65134-y",
"10.1038/s41467-025-65049-8",
]

def fetch(doi):
    for attempt in range(3):
        try:
            r = requests.get("https://api.crossref.org/works/" + doi,
                headers={"User-Agent": "ZEB-axis-manuscript/1.0 (academic reanalysis)"},
                timeout=30)
            r.raise_for_status()
            m = r.json()["message"]
            if m["DOI"].lower() != doi.lower():
                raise ValueError("DOI identity mismatch")
            authors = m.get("author", [])
            return {"doi": m["DOI"], "title": re.sub("<[^>]+>", "", m.get("title", [""])[0]),
                "authors": [{"given": a.get("given", ""), "family": a.get("family", "")}
                            for a in authors],
                "journal": m.get("container-title", [""])[0],
                "year": m.get("published", {}).get("date-parts", [[None]])[0][0],
                "volume": m.get("volume", ""), "issue": m.get("issue", ""),
                "page": m.get("page", m.get("article-number", "")),
                "url": "https://doi.org/" + m["DOI"],
                "verification": "Crossref DOI and title identity checked",
                "crossref_sha256": hashlib.sha256(r.content).hexdigest()}
        except Exception as exc:
            error = str(exc)
            time.sleep(attempt + 1)
    return {"doi": doi, "error": error}

if __name__ == "__main__":
    if len({x.lower() for x in DOIS}) != len(DOIS):
        raise SystemExit("Duplicate DOI in curated candidates")
    with ThreadPoolExecutor(max_workers=6) as pool:
        records = list(pool.map(fetch, DOIS))
    (DEST / "references_verified.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
    failures = [x for x in records if "error" in x]
    print("verified", len(records) - len(failures), "failed", len(failures))
    for x in failures:
        print(x["doi"], x["error"])
    for i, x in enumerate(records, 1):
        print(i, x["doi"], x.get("title", "ERROR"))
