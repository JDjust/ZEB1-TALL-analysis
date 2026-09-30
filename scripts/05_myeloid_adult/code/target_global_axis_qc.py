#!/usr/bin/env python3
"""Check whether TARGET ZEB1 tracks broad sample expression/QC axes.

Uses the existing processed log2(TPM+1) matrix, not raw counts. Any association
is a diagnostic flag rather than proof of its technical or biological cause.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.stats import rankdata
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/validation"


def main():
    meta = pd.read_csv(DATA / "target_samples.tsv", sep="\t")
    expr = pd.read_csv(DATA / "target_log2tpm_filtered.tsv.gz", sep="\t")
    if expr.columns[0] != "gene" or expr.gene.duplicated().any():
        raise ValueError("gene schema or duplicate gene IDs")
    samples = expr.columns[1:].tolist()
    if len(samples) != 265 or meta.usi.duplicated().any() or set(samples) != set(meta.usi):
        raise ValueError("TARGET patient/sample alignment failed")
    meta = meta.set_index("usi").loc[samples]
    x = expr.drop(columns="gene").to_numpy(dtype=np.float32)
    if not np.isfinite(x).all() or (x < 0).any():
        raise ValueError("unexpected expression scale or nonfinite values")
    zrow = expr.index[expr.gene.eq("ZEB1")]
    if len(zrow) != 1 or not np.allclose(x[zrow[0]], meta.ZEB1_log, atol=.02):
        raise ValueError("ZEB1 source row disagrees with metadata")
    qc = pd.DataFrame(index=samples)
    qc.index.name = "usi"
    qc["ZEB1_log"] = meta.ZEB1_log.to_numpy()
    qc["subtype"] = meta.subtype.to_numpy()
    qc["age_years"] = meta.age_years.to_numpy()
    qc["sample_type"] = meta.sample_type.to_numpy()
    qc["n_genes_detected_TPM_gt0"] = (x > 0).sum(axis=0)
    qc["median_log2tpm"] = np.median(x, axis=0)
    qc["mean_log2tpm"] = x.mean(axis=0)
    qc["frac_genes_log2tpm_gt1"] = (x > 1).mean(axis=0)
    raw_path = DATA / "target_star_raw_sample_qc.tsv"
    if raw_path.is_file():
        raw_qc = pd.read_csv(raw_path, sep="\t").set_index("usi")
        if len(raw_qc) != 265 or set(raw_qc.index) != set(qc.index):
            raise ValueError("STAR raw QC sample alignment failed")
        raw_metrics = ["assigned_unstranded", "detected_genes_count_gt0",
                       "detected_protein_coding_count_gt0", "mitochondrial_count_fraction",
                       "protein_coding_count_fraction", "assigned_fraction_all_reported",
                       "N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]
        qc = qc.join(raw_qc[raw_metrics])
    else:
        raw_metrics = []
    # Gene-centered PCA; ZEB1 is excluded to avoid a trivially circular PC.
    keep = expr.gene.ne("ZEB1").to_numpy()
    matrix = x[keep, :].T
    pca = PCA(n_components=3, svd_solver="randomized", random_state=20260924)
    pcs = pca.fit_transform(matrix)
    for i in range(3):
        qc[f"expression_PC{i+1}"] = pcs[:, i]
    qc.to_csv(DATA / "target_global_axis_sample_qc.tsv", sep="\t")
    metrics = ["n_genes_detected_TPM_gt0", "median_log2tpm", "mean_log2tpm",
               "frac_genes_log2tpm_gt1", "expression_PC1", "expression_PC2",
               "expression_PC3"] + raw_metrics
    rows = []
    for metric in metrics:
        r, p = spearmanr(qc.ZEB1_log, qc[metric])
        rows.append({"metric": metric, "n": len(qc), "rho": r, "p_exploratory": p,
                     "unit": "patient/sample", "input": "STAR raw counts" if metric in raw_metrics
                     else "processed log2(TPM+1)"})
    out = pd.DataFrame(rows)
    out.to_csv(DATA / "target_global_axis_associations.tsv", sep="\t", index=False)
    # Match the 242 known-subtype sample design used in module 4. A rank-based
    # residual association is diagnostic; conditioning on expression PC1 may
    # remove biology as well as technical variation.
    known = qc[qc.subtype.notna() & qc.subtype.ne("Unknown") & qc.age_years.notna()].copy()
    if len(known) != 242:
        raise ValueError("known-subtype module 4 sample count changed")
    cov = pd.get_dummies(known["subtype"], drop_first=True, dtype=float)
    cov["age_years"] = known.age_years.to_numpy(dtype=float)
    cov["peripheral_blood"] = known.sample_type.str.contains("Peripheral").astype(float).to_numpy()
    design = np.column_stack([np.ones(len(known)), cov.to_numpy(dtype=float)])
    def residual_rank(v):
        y = rankdata(np.asarray(v, dtype=float))
        return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    z_resid = residual_rank(known.ZEB1_log)
    partial_rows = []
    for metric in metrics:
        r, p = spearmanr(z_resid, residual_rank(known[metric]))
        partial_rows.append({"metric": metric, "n": len(known), "partial_rho_rank_residual": r,
                             "p_exploratory": p,
                             "covariates": "molecular subtype, age_years, blood vs marrow"})
    partial = pd.DataFrame(partial_rows)
    partial.to_csv(DATA / "target_global_axis_partial.tsv", sep="\t", index=False)
    summary = {"n_samples": len(qc), "n_genes": len(expr),
               "pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
               "caution": "STAR count summaries are available, but alignment QC, batch and tumor purity remain unmeasured"}
    (DATA / "target_global_axis_run.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(out.to_string(index=False))
    print(partial.to_string(index=False))
    print(summary)


if __name__ == "__main__":
    main()
