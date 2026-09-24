"""Read the author's Strict OOXML workbook without changing its original bytes."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import xml.etree.ElementTree as ET
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/validation/gse287751"
path = OUT / "GSM9286835_scRNAseq_perCell_assignment_to_hto_cluster_umap.xlsx"
with ZipFile(path) as archive:
    strings = ["".join(node.itertext()) for node in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
    def sheet(name):
        records = []
        for row in ET.fromstring(archive.read(name)).findall(".//{*}row"):
            cells = {}
            for cell in row.findall("{*}c"):
                value = cell.find("{*}v")
                if value is not None:
                    cells[cell.attrib["r"].rstrip("0123456789")] = strings[int(value.text)] if cell.attrib.get("t") == "s" else value.text
            if cells:
                records.append(cells)
        return records
    notes = sheet("xl/worksheets/sheet1.xml")
    records = sheet("xl/worksheets/sheet2.xml")
header, body = records[0], records[1:]
metadata = pd.DataFrame([{header[key]: value for key, value in row.items()} for row in body])
assert metadata["index"].is_unique and metadata.hto.nunique() == 8
metadata.to_csv(OUT / "author_cell_metadata.tsv", sep="\t", index=False)
metadata["barcode"] = metadata["index"].str.replace(r"-1$", "", regex=True)
raw = pd.read_csv(ROOT.parent / "modules/module10/tables/M10_r2_HTO_assignment_raw.tsv", sep="\t")
matched = metadata.merge(raw, on="barcode", how="left", validate="one_to_one")
cross = matched.groupby(["hto", "hto_tag"], dropna=False).size().reset_index(name="n_cells")
cross.to_csv(OUT / "author_condition_to_raw_hto.tsv", sep="\t", index=False)
counts = metadata.groupby("hto").size().reset_index(name="author_retained_cells")
counts.to_csv(OUT / "author_condition_counts.tsv", sep="\t", index=False)
(OUT / "metadata_provenance.json").write_text(json.dumps({
    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    "format": "Strict OOXML; parsed using original XML namespaces",
    "author_readme": notes, "n_cells": len(metadata),
    "n_conditions": metadata.hto.nunique(),
    "warning": "Eight conditions are not eight biological replicates"}, indent=2), encoding="utf-8")
print(counts.to_string(index=False))
print(cross.to_string(index=False))
