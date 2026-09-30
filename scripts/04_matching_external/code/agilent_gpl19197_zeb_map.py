#!/usr/bin/env python3
"""Sequence-level ZEB1/ZEB2/TLX1/LMO2 mapping for GPL19197.

GSE62141 (TLX1 siRNA, ALL-SIL) and GSE62143 (TLX1 overexpression, CD34+
thymocytes) were left unmapped because the limma tables are Agilent probe
IDs. GPL19197 stores the 60-mer Sequence next to Gene Symbol. A probe is
sequence-supported when that 60-mer, or its reverse complement, occurs
exactly in a current RefSeq of the named gene. Symbol-only rows are reported
separately and are not used for the bidirectional call.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "data" / "validation" / "zeb_family_state"
GPL = VAL / "GPL19197_data.txt"
FASTA = VAL / "refseq_zeb_tlx_lmo.fa"
KD = ROOT.parent / "modules" / "module9" / "tables" / "M9_gse62141_KD.tsv"
OE = ROOT.parent / "modules" / "module9" / "tables" / "M9_gse62143_TLX1_OE.tsv"
GENES = {"ZEB1", "ZEB2", "TLX1", "LMO2"}
ACCESSION_GENE = {
    "NM_001128128": "ZEB1",
    "NM_030751": "ZEB1",
    "NM_014795": "ZEB2",
    "NM_001171653": "ZEB2",
    "NM_005521": "TLX1",
    "NM_005574": "LMO2",
}


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def parse_fasta(path: Path) -> dict[str, str]:
    out, acc, buf = {}, None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if acc:
                out[acc] = "".join(buf).upper().replace("U", "T")
            acc = line[1:].split()[0].split(".")[0]
            buf = []
        else:
            buf.append(line.strip())
    if acc:
        out[acc] = "".join(buf).upper().replace("U", "T")
    return out


def load_platform(path: Path) -> pd.DataFrame:
    rows = []
    header = None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        started = False
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                started = True
                continue
            if not started:
                continue
            if line.startswith("!platform_table_end"):
                break
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    df = pd.DataFrame(rows).rename(columns={"ID": "probe", "Gene Symbol": "symbol", "Sequence": "sequence"})
    df["sequence"] = df["sequence"].str.upper().str.replace("U", "T", regex=False)
    return df


def hamming_best(oligo: str, target: str) -> int:
    n = len(oligo)
    if len(target) < n:
        return n
    best = n
    for seq in (oligo, rc(oligo)):
        for i in range(len(target) - n + 1):
            mm = sum(a != b for a, b in zip(seq, target[i:i + n]))
            if mm < best:
                best = mm
            if best == 0:
                return 0
    return best


def main() -> None:
    targets = parse_fasta(FASTA)
    missing = set(ACCESSION_GENE) - set(targets)
    if missing:
        raise SystemExit(f"missing RefSeqs: {missing}")
    platform = load_platform(GPL)
    platform = platform[platform["sequence"].str.fullmatch(r"[ACGTN]+", na=False)].copy()

    hits = []
    by_gene = {gene: [(acc, targets[acc]) for acc, g in ACCESSION_GENE.items() if g == gene] for gene in GENES}
    for rec in platform.itertuples(index=False):
        seq = rec.sequence
        matched = []
        for gene, seqs in by_gene.items():
            for acc, target in seqs:
                if seq in target or rc(seq) in target:
                    matched.append(f"{gene}:{acc}")
        symbol = str(rec.symbol).split("|")[0].strip()
        symbol_gene = symbol if symbol in GENES else ""
        if matched or symbol_gene:
            hits.append({
                "probe": rec.probe,
                "symbol": symbol,
                "sequence": seq,
                "sequence_hits": ";".join(matched),
                "sequence_gene": ";".join(sorted({m.split(":")[0] for m in matched})),
                "symbol_in_query": symbol_gene,
            })
    mapped = pd.DataFrame(hits)
    # Near-match audit only for symbol hits that failed the exact sequence test.
    near = []
    for rec in mapped.itertuples(index=False):
        if rec.sequence_gene or not rec.symbol_in_query:
            continue
        gene = rec.symbol_in_query
        best = min(hamming_best(rec.sequence, target) for _, target in by_gene[gene])
        near.append({"probe": rec.probe, "symbol": gene, "best_mismatches_to_symbol_refseq": best})
    near_df = pd.DataFrame(near)
    mapped.to_csv(VAL / "gpl19197_candidate_probes.tsv", sep="\t", index=False)
    near_df.to_csv(VAL / "gpl19197_symbol_without_exact_match.tsv", sep="\t", index=False)

    supported = mapped[mapped["sequence_gene"].ne("")].copy()
    kd = pd.read_csv(KD, sep="\t").rename(columns={"gene": "probe"})
    oe = pd.read_csv(OE, sep="\t").rename(columns={"gene": "probe"})
    effects = supported.merge(kd[["probe", "log2FoldChange", "pvalue", "padj"]], on="probe", how="left")
    effects = effects.merge(
        oe[["probe", "log2FoldChange", "pvalue", "padj"]],
        on="probe", how="left", suffixes=("_KD", "_OE"),
    )
    effects.to_csv(VAL / "gpl19197_sequence_supported_effects.tsv", sep="\t", index=False)

    summary = []
    for gene in ("TLX1", "ZEB1", "ZEB2", "LMO2"):
        part = effects[effects["sequence_gene"].str.contains(gene)]
        summary.append({
            "gene": gene,
            "n_exact_probes": int(len(part)),
            "KD_median_log2FC": float(part["log2FoldChange_KD"].median()) if len(part) else None,
            "KD_min_log2FC": float(part["log2FoldChange_KD"].min()) if len(part) else None,
            "KD_max_log2FC": float(part["log2FoldChange_KD"].max()) if len(part) else None,
            "KD_n_positive": int((part["log2FoldChange_KD"] > 0).sum()) if len(part) else 0,
            "OE_median_log2FC": float(part["log2FoldChange_OE"].median()) if len(part) else None,
            "OE_min_log2FC": float(part["log2FoldChange_OE"].min()) if len(part) else None,
            "OE_max_log2FC": float(part["log2FoldChange_OE"].max()) if len(part) else None,
            "OE_n_negative": int((part["log2FoldChange_OE"] < 0).sum()) if len(part) else 0,
            "n_symbol_conflict": int(((mapped.symbol.eq(gene)) & (~mapped.sequence_gene.str.contains(gene))).sum()),
        })
    pd.DataFrame(summary).to_csv(VAL / "gpl19197_bidirectional_summary.tsv", sep="\t", index=False)
    (VAL / "gpl19197_map_methods.json").write_text(json.dumps({
        "platform": "GPL19197 Agilent-041648 CMGG Human V1.1 60k",
        "series": ["GSE62141 TLX1 siRNA ALL-SIL", "GSE62143 TLX1 OE CD34+ thymocytes"],
        "gpl_sha256": hashlib.sha256(GPL.read_bytes()).hexdigest(),
        "fasta_sha256": hashlib.sha256(FASTA.read_bytes()).hexdigest(),
        "rule": "Exact 60-mer or reverse complement contained in the listed RefSeq. Mismatch scan is limited to symbol-annotated probes that failed this test.",
        "not_genome_wide": "Uniqueness is relative to these six RefSeqs only.",
        "limma_sources": [str(KD), str(OE)],
    }, indent=2), encoding="utf-8")
    print(pd.DataFrame(summary).to_string(index=False))
    print(effects[["probe", "symbol", "sequence_gene", "log2FoldChange_KD", "padj_KD", "log2FoldChange_OE", "padj_OE"]].to_string(index=False))
    if len(near_df):
        print(near_df.to_string(index=False))


if __name__ == "__main__":
    main()
