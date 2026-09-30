"""Cache GEO accession citation links and PubMed identity records.

Accession links identify candidate source publications; they do not establish
that a particular analysis or claim is supported by that publication.
"""
from pathlib import Path
import hashlib
import json
import re
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/reference_records"
OUT.mkdir(parents=True, exist_ok=True)


def download(url, cache):
    if cache.exists():
        return cache.read_bytes()
    request = urllib.request.Request(url, headers={"User-Agent": "ZEB1-public-data-audit/1.0"})
    with urllib.request.urlopen(request, timeout=15) as response:
        raw = response.read()
    cache.write_bytes(raw)
    return raw


def plain(element, tag):
    node = element.find(tag)
    return "" if node is None else "".join(node.itertext()).strip()


def main():
    manuscript = "\n".join((ROOT / "code" / name).read_text(encoding="utf-8")
                           for name in ["build_manuscript.py", "build_supplementary.py"])
    accessions = sorted(set(re.findall(r"GSE\d+", manuscript)))
    registry = []
    for accession in accessions:
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=brief"
        record = {"accession": accession, "geo_url": url}
        try:
            raw = download(url, OUT / f"{accession}.soft.txt")
            text = raw.decode("utf-8", errors="replace")
            if f"^SERIES = {accession}" not in text:
                raise ValueError("Response does not identify requested GEO series")
            record.update(geo_title=re.findall(r"^!Series_title = (.+)$", text, re.M),
                          pmids=re.findall(r"^!Series_pubmed_id = (\d+)\s*$", text, re.M),
                          sha256=hashlib.sha256(raw).hexdigest(), status="retrieved")
        except Exception as exc:
            record.update(status="unresolved", error=str(exc), pmids=[])
        registry.append(record)
        print(accession, record["status"], record["pmids"], flush=True)
        (OUT / "geo_citation_links.json").write_text(json.dumps(registry, indent=2), encoding="utf-8")
        time.sleep(.4)
    pmids = sorted({pmid for r in registry for pmid in r["pmids"]})
    if not pmids:
        return
    query = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"})
    raw = download("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + query,
                   OUT / ("pubmed_" + hashlib.sha256(query.encode()).hexdigest()[:12] + ".xml"))
    records = []
    for article in ET.fromstring(raw).findall("PubmedArticle"):
        citation = article.find("MedlineCitation")
        content = citation.find("Article")
        ids = {x.attrib.get("IdType"): x.text for x in article.findall("PubmedData/ArticleIdList/ArticleId")}
        records.append({"pmid": plain(citation, "PMID"), "title": plain(content, "ArticleTitle"),
                        "doi": ids.get("doi", ""), "journal": plain(content, "Journal/Title"),
                        "year": plain(content, "Journal/JournalIssue/PubDate/Year") or
                                plain(content, "Journal/JournalIssue/PubDate/MedlineDate"),
                        "volume": plain(content, "Journal/JournalIssue/Volume"),
                        "pages": plain(content, "Pagination/MedlinePgn"),
                        "authors": [plain(a, "LastName") for a in content.findall("AuthorList/Author")],
                        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{plain(citation, 'PMID')}/"})
    (OUT / "pubmed_identities.json").write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"PubMed identities: {len(records)}/{len(pmids)}", flush=True)


if __name__ == "__main__":
    main()
