"""Build an accession/access/citation registry from cached GEO source records."""
from pathlib import Path
import csv
import json
import re

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "data/reference_records"


def main():
    links = json.loads((RECORDS / "geo_citation_links.json").read_text(encoding="utf-8"))
    papers = {r["pmid"]: r for r in json.loads((RECORDS / "pubmed_identities.json").read_text(encoding="utf-8"))}
    rows = []
    for record in links:
        accession = record["accession"]
        raw_path = RECORDS / f"{accession}.soft.txt"
        raw = raw_path.read_text(encoding="utf-8") if raw_path.exists() else ""
        def fields(name):
            return re.findall(rf"^!Series_{name} = (.+)$", raw, re.M)
        design = " ".join(fields("overall_design"))
        relations = fields("relation")
        if "raw data are not available" in design.lower():
            access = "Submitter states raw data unavailable due to patient privacy"
        elif "dbgap" in design.lower() or any("dbgap" in x.lower() for x in relations):
            access = "Controlled raw-data repository identified; authorization not inferred"
        elif any("SRA" in x for x in relations):
            access = "SRA relation listed; individual raw-file accessibility not verified"
        else:
            access = "Raw-data access unresolved from this record"
        notes = "GEO-linked citation; supporting passages require source-level review"
        if accession == "GSE227122":
            notes += "; GEO currently links a 2025 atlas; 2023 primary DOI 10.1038/s41598-023-39152-z is separately verified in core_reference_identity.json"
        if accession == "GSE248287":
            notes += "; immune/leukemia fractions sorted and pooled; observed proportions are not native burden"
        rows.append({"accession": accession, "title": " | ".join(record.get("geo_title", [])),
                     "source_url": record["geo_url"], "status": record["status"],
                     "linked_pmids": ";".join(record["pmids"]),
                     "linked_dois": ";".join(papers[p]["doi"] for p in record["pmids"] if p in papers),
                     "raw_access_evidence": access, "deposited_design": design,
                     "supplementary_urls": " | ".join(fields("supplementary_file")),
                     "relations": " | ".join(relations),
                     "geo_source_sha256": record.get("sha256", ""), "notes": notes})
    with (ROOT / "data/dataset_registry.tsv").open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Dataset registry: {len(rows)} accessions")


if __name__ == "__main__":
    main()
