#!/usr/bin/env python3
"""Freeze GSE29959 NOTCH signatures, then score the 8 DepMap T-ALL models.

Contrasts, per cell line, from the GEO design (one array each):
  parental DMSO, parental GSI (Compound E),
  ICN+GSI, c-Myc+GSI.
A drop is log2(GSI) - log2(DMSO) <= -0.5.
Rescue means the GSI sample is restored at least halfway to DMSO.
A gene is called only when 4 of 5 lines agree.
This threshold is not chosen from the DepMap AUC.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

MATRIX = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE29959/GSE29959_series_matrix.txt.gz")
ANNOT = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE13159/GPL570.annot.gz")
DEPMAP = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
OUT = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/_zeb1_forensic")

DROP = -0.5
MIN_LINES = 4
LINES = ["ALLSIL", "DND41", "HPBALL", "KOPTK1", "TALL-1"]

GROUP = {
    "ACH-000981": ("DND-41", "dependent"),
    "ACH-000197": ("TALL-1", "dependent"),
    "ACH-000995": ("JURKAT", "dependent"),
    "ACH-000953": ("SUP-T1", "dependent"),
    "ACH-000937": ("PF-382", "nondependent"),
    "ACH-000519": ("PEER", "nondependent"),
    "ACH-000101": ("KE-37", "nondependent"),
    "ACH-001737": ("CCRF-HSB-2-DM", "nondependent"),
}
CONTROLS = ["HES1", "HES4", "HEY1", "DTX1", "NRARP", "MYC", "PTCRA", "NOTCH3", "IL7R", "SHQ1"]


def gene_symbol(col: str) -> str:
    return str(col).split(" (")[0].strip()


def auc(values: pd.Series, dep_ids, non_ids) -> float:
    d = values.reindex(dep_ids).to_numpy(float)
    n = values.reindex(non_ids).to_numpy(float)
    if np.isnan(d).any() or np.isnan(n).any():
        return np.nan
    wins = ties = 0
    for a in d:
        for b in n:
            if a > b:
                wins += 1
            elif a == b:
                ties += 1
    return (wins + 0.5 * ties) / (len(d) * len(n))


def load_expr() -> pd.DataFrame:
    # GEO matrix: comment lines, then a table. Values are submitted processed intensities.
    df = pd.read_csv(MATRIX, sep="\t", comment="!", quotechar='"')
    df = df.rename(columns={df.columns[0]: "probe"})
    df = df[df["probe"].ne("ID_REF")]
    df = df.set_index("probe")
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def sample_map(columns) -> dict:
    # Titles are in GEO order: 6 samples x 5 lines, matching GSM order in the matrix.
    # Confirmed from !Sample_title before this script is applied to DepMap.
    roles = ["DMSO", "GSI", "ICN_DMSO", "ICN_GSI", "MYC_DMSO", "MYC_GSI"]
    if len(columns) != 30:
        raise SystemExit(f"expected 30 samples, got {len(columns)}")
    mapping = {}
    for i, col in enumerate(columns):
        line = LINES[i // 6]
        role = roles[i % 6]
        mapping[col] = (line, role)
    return mapping


def load_symbols() -> pd.Series:
    # GPL570.annot.gz: comment header, then ID / Gene symbol.
    sym = {}
    with ANNOT.open("rb") as raw:
        import gzip
        with gzip.open(raw, "rt", errors="replace") as fh:
            header = None
            for line in fh:
                if line.startswith("#") or line.startswith("!") or line.strip() in {"Annotation", "^Annotation"}:
                    continue
                if not line.startswith("ID\t") and header is None:
                    continue
                if header is None:
                    header = line.rstrip("\n").split("\t")
                    id_i = header.index("ID")
                    try:
                        g_i = header.index("Gene symbol")
                    except ValueError:
                        g_i = header.index("Gene Symbol")
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(id_i, g_i):
                    continue
                gene = parts[g_i].split("///")[0].strip()
                if gene and gene not in {"---", "NA"}:
                    sym[parts[id_i]] = gene
    return pd.Series(sym, name="symbol")


def line_effects(expr: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    mapping = sample_map(expr.columns)
    # average probes per gene
    ann = symbols.reindex(expr.index).dropna()
    use = expr.loc[ann.index].copy()
    use["symbol"] = ann.values
    gene = use.groupby("symbol").mean(numeric_only=True)
    rows = []
    for line in LINES:
        cols = {role: col for col, (ln, role) in mapping.items() if ln == line}
        d = gene[cols["DMSO"]]
        g = gene[cols["GSI"]]
        i = gene[cols["ICN_GSI"]]
        m = gene[cols["MYC_GSI"]]
        # Submitted values are unlogged and include negatives. Floor at 1, then log2.
        # The 0.5 cutoff is then a 1.4-fold change, fixed before DepMap scoring.
        d = np.log2(np.clip(d, 1, None))
        g = np.log2(np.clip(g, 1, None))
        i = np.log2(np.clip(i, 1, None))
        m = np.log2(np.clip(m, 1, None))
        drop = g - d
        need = -drop
        icn_back = i - g
        myc_back = m - g
        part = pd.DataFrame({
            f"{line}_drop": drop,
            f"{line}_gsi_down": drop <= DROP,
            f"{line}_icn_rescue": (drop <= DROP) & (icn_back >= 0.5 * need),
            f"{line}_myc_rescue": (drop <= DROP) & (myc_back >= 0.5 * need),
            f"{line}_gsi_changed": drop.abs() >= 0.5,
            f"{line}_icn_reverses": (
                ((drop <= DROP) & (icn_back >= 0.5 * need))
                | ((drop >= 0.5) & ((g - i) >= 0.5 * drop))
            ),
        })
        rows.append(part)
    return pd.concat(rows, axis=1)


def call_sets(effects: pd.DataFrame) -> dict[str, list[str]]:
    down = effects[[c for c in effects.columns if c.endswith("_gsi_down")]].sum(axis=1)
    icn = effects[[c for c in effects.columns if c.endswith("_icn_rescue")]].sum(axis=1)
    myc = effects[[c for c in effects.columns if c.endswith("_myc_rescue")]].sum(axis=1)
    changed = effects[[c for c in effects.columns if c.endswith("_gsi_changed")]].sum(axis=1)
    reverses = effects[[c for c in effects.columns if c.endswith("_icn_reverses")]].sum(axis=1)
    direct = effects.index[(down >= MIN_LINES) & (icn >= MIN_LINES) & (myc <= 1)]
    via_myc = effects.index[(down >= MIN_LINES) & (icn >= MIN_LINES) & (myc >= MIN_LINES)]
    nonspecific = effects.index[(changed >= MIN_LINES) & (reverses <= 1)]
    # keep the three definitions disjoint for the gate
    nonspecific = nonspecific.difference(direct).difference(via_myc)
    return {
        "NOTCH_direct": sorted(direct),
        "NOTCH_to_MYC": sorted(via_myc),
        "GSI_nonspecific": sorted(nonspecific),
    }


def score_depmap(genesets: dict[str, list[str]]) -> pd.DataFrame:
    df = pd.read_csv(DEPMAP, low_memory=False)
    df = df[df["ModelID"].astype(str).isin(GROUP)].copy()
    df = df.rename(columns={c: gene_symbol(c) for c in df.columns if c != "ModelID"})
    df["ModelID"] = df["ModelID"].astype(str)
    mat = df.groupby("ModelID").mean(numeric_only=True).reindex(list(GROUP))
    dep = [i for i, (_, g) in GROUP.items() if g == "dependent"]
    non = [i for i, (_, g) in GROUP.items() if g == "nondependent"]
    rows = []
    for name, genes in genesets.items():
        present = [g for g in genes if g in mat.columns]
        if len(present) < 3:
            rows.append({"signature": name, "n_genes_defined": len(genes), "n_genes_in_depmap": len(present), "auc_dep_higher": np.nan})
            continue
        sub = mat[present]
        z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
        score = z.mean(axis=1)
        rows.append({
            "signature": name,
            "n_genes_defined": len(genes),
            "n_genes_in_depmap": len(present),
            "auc_dep_higher": auc(score, dep, non),
            **{GROUP[i][0]: float(score.loc[i]) for i in GROUP},
        })
    return pd.DataFrame(rows)


def main() -> None:
    expr = load_expr()
    print("matrix", expr.shape, "median", float(np.nanmedian(expr.to_numpy())))
    symbols = load_symbols()
    effects = line_effects(expr, symbols)
    print("control gene drops (GSI-DMSO, log2); negative means GSI down")
    for gene in CONTROLS:
        if gene not in effects.index:
            print(gene, "ABSENT")
            continue
        drops = [float(effects.loc[gene, f"{ln}_drop"]) for ln in LINES]
        flags = []
        for ln in LINES:
            flags.append(
                f"{ln[0]}:{int(effects.loc[gene, f'{ln}_gsi_down'])}/"
                f"{int(effects.loc[gene, f'{ln}_icn_rescue'])}/"
                f"{int(effects.loc[gene, f'{ln}_myc_rescue'])}"
            )
        print(gene, " ".join(f"{v:+.2f}" for v in drops), " | ", " ".join(flags))
    genesets = call_sets(effects)
    for name, genes in genesets.items():
        print(f"{name} n={len(genes)} head={genes[:12]}")
        hit = [g for g in CONTROLS if g in genes]
        print("  controls inside:", hit)
    sets_path = OUT / "gse29959_signatures.tsv"
    pd.DataFrame(
        [(name, g) for name, genes in genesets.items() for g in genes],
        columns=["signature", "gene"],
    ).to_csv(sets_path, sep="\t", index=False)
    direct_wo_overlap = [g for g in genesets["NOTCH_direct"] if g not in {"HES4", "DTX1"}]
    genesets["NOTCH_direct_without_HES4_DTX1"] = direct_wo_overlap
    scored = score_depmap(genesets)
    scored.to_csv(OUT / "gse29959_depmap_auc.tsv", sep="\t", index=False)
    print(scored.round(3).to_string(index=False))
    # Per-gene AUC for the frozen direct set, so one familiar gene cannot hide.
    df = pd.read_csv(DEPMAP, low_memory=False)
    df = df[df["ModelID"].astype(str).isin(GROUP)].copy()
    df = df.rename(columns={c: gene_symbol(c) for c in df.columns if c != "ModelID"})
    mat = df.groupby(df["ModelID"].astype(str)).mean(numeric_only=True)
    dep = [i for i, (_, g) in GROUP.items() if g == "dependent"]
    non = [i for i, (_, g) in GROUP.items() if g == "nondependent"]
    print("per-gene AUC in NOTCH_direct")
    for g in genesets["NOTCH_direct"]:
        if g not in mat.columns:
            print(g, "missing")
            continue
        print(f"{g}\t{auc(mat[g], dep, non):.3f}")


if __name__ == "__main__":
    main()
