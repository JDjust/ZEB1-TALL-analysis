#!/usr/bin/env python3
"""Extract sample-level QC proxies from public GDC STAR gene-count TSVs.

Run on Hulu via Slurm; parses one file at a time and writes aggregate QC only.
Counts are unstranded; mapping/alignment quality beyond these summary rows is
not available here. No sequence reads are accessed.
"""
from pathlib import Path
import argparse
import csv
import json
import hashlib


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def scan(path):
    specials = {}
    assigned = mt = protein = 0
    detected = detected_protein = 0
    with path.open("r", encoding="utf-8") as fh:
        first = fh.readline().strip()
        if not first.startswith("# gene-model:"):
            raise ValueError(f"missing gene model header: {path}")
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"]
            val = int(row["unstranded"])
            if gene_id.startswith("N_"):
                specials[gene_id] = val
                continue
            assigned += val
            detected += val > 0
            if row["gene_type"] == "protein_coding":
                protein += val
                detected_protein += val > 0
            if row["gene_name"].startswith("MT-"):
                mt += val
    required = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}
    if not required.issubset(specials) or assigned == 0:
        raise ValueError(f"incomplete STAR summary: {path}")
    total = assigned + sum(specials.values())
    return {"gene_model": first.split(":", 1)[1].strip(),
            "assigned_unstranded": assigned,
            "detected_genes_count_gt0": detected,
            "detected_protein_coding_count_gt0": detected_protein,
            "mitochondrial_count_fraction": mt / assigned,
            "protein_coding_count_fraction": protein / assigned,
            "assigned_fraction_all_reported": assigned / total,
            **specials}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--raw-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    with args.metadata.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 265 or len({r["usi"] for r in rows}) != 265:
        raise ValueError("unexpected TARGET metadata sample count/IDs")
    result = []
    for index, sample in enumerate(rows, 1):
        file_name = sample["file_name"]
        if Path(file_name).name != file_name:
            raise ValueError("metadata file_name is not a basename")
        source = args.raw_dir / file_name
        if not source.is_file():
            raise FileNotFoundError(source)
        result.append({"usi": sample["usi"], "file_name": file_name, **scan(source)})
        if index % 25 == 0:
            print(f"processed {index}/{len(rows)}", flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    out = args.output / "target_star_raw_sample_qc.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(result[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(result)
    manifest = {"n_samples": len(result), "metadata_sha256": sha256(args.metadata),
                "raw_directory": str(args.raw_dir), "output_sha256": sha256(out),
                "count_field": "unstranded", "status": "complete"}
    (args.output / "target_star_raw_qc_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
