"""Extend the eight existing references with selected, cached PubMed identities."""
from pathlib import Path
import json
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SELECTED = {
    "GSE13159/MILE": ["20406941"], "TARGET": ["28671688"],
    "GSE107011": ["30726743"], "GSE142522": ["32667968"],
    "GSE195812": ["36367948"], "GSE206710/Park atlas": ["32079746"],
    "GSE62156": ["25301704"], "GSE26713": ["21481790"],
    "GSE110635/GSE110636/GSE110637": ["29954933", "33259601"],
    "GSE272023": ["39841000"], "GSE69954": ["26928953", "32345961"],
    "GSE188225": ["36920307"], "GSE49164": ["23926305"],
    "GSE186943": ["35851847"], "GSE70734": ["26108691"],
    "GSE243915/GSE280728": ["39719705", "39880094"],
    "GSE287751/GSE287757": ["41824555"],
    "GSE260697": ["38553571"],
    "GSE154675": ["40637766"],
}


def text(node, query):
    child = node.find(query)
    return "" if child is None else "".join(child.itertext()).strip()


def main():
    refs = json.loads((ROOT / "data/core_reference_identity.json").read_text(encoding="utf-8"))["records"]
    records = {}
    for source in (ROOT / "data/reference_records").glob("pubmed*.xml"):
        for article in ET.fromstring(source.read_bytes()).findall("PubmedArticle"):
            citation = article.find("MedlineCitation")
            content = citation.find("Article")
            pmid = text(citation, "PMID")
            ids = {x.attrib.get("IdType"): x.text for x in article.findall("PubmedData/ArticleIdList/ArticleId")}
            records[pmid] = {"pmid": pmid, "doi": ids.get("doi", ""), "title": text(content, "ArticleTitle"),
                             "authors_family": [text(a, "LastName") or text(a, "CollectiveName") for a in content.findall("AuthorList/Author")],
                             "journal": text(content, "Journal/Title"), "year": text(content, "Journal/JournalIssue/PubDate/Year"),
                             "volume": text(content, "Journal/JournalIssue/Volume"),
                             "page": text(content, "Pagination/MedlinePgn"), "article_number": "",
                             "source_record": str(source.relative_to(ROOT)).replace("\\", "/"),
                             "verification_level": "PubMed title/PMID/DOI identity; dataset-source linkage requires assay-level reading"}
    mapping = []
    for dataset, pmids in SELECTED.items():
        numbers = []
        for pmid in pmids:
            ref = records[pmid]
            if not ref["title"] or not ref["doi"]:
                raise ValueError(f"Incomplete identity: PMID {pmid}")
            index = next((i for i, old in enumerate(refs) if old["doi"].lower() == ref["doi"].lower()), None)
            if index is None:
                refs.append(ref)
                index = len(refs) - 1
            numbers.append(index + 1)
        mapping.append({"dataset": dataset, "reference_numbers": numbers, "pmids": pmids})
    geo_path = ROOT / "data/reference_records/GSE144035.soft.txt"
    geo_text = geo_path.read_text(encoding="utf-8")
    title = next(line.split(" = ", 1)[1].strip() for line in geo_text.splitlines() if line.startswith("!Series_title = "))
    refs.append({"record_type": "dataset", "accession": "GSE144035", "title": title,
                 "authors_family": ["Vo A", "Wong NC", "Shields BJ", "McCormack MP"],
                 "journal": "NCBI Gene Expression Omnibus", "year": "2020", "volume": "",
                 "page": "", "article_number": "GSE144035", "doi": "",
                 "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144035",
                 "source_record": str(geo_path.relative_to(ROOT)).replace("\\", "/"),
                 "verification_level": "GEO accession/title/contributor record; primary paper linkage unresolved"})
    mapping.append({"dataset": "GSE144035", "reference_numbers": [len(refs)], "pmids": []})
    output = {"records": refs, "dataset_citations": mapping,
              "unresolved_primary_publications": ["GSE144035"]}
    (ROOT / "data/reference_catalog.json").write_text(json.dumps(output, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"References: {len(refs)}; accession/source groups: {len(mapping)}")


if __name__ == "__main__":
    main()
