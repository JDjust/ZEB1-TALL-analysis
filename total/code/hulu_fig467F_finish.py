#!/usr/bin/env python3
"""Finish Fig4/6/7 F analyses: TARGET Cox, GDSC MTX, Park Human T cells."""
from __future__ import annotations

import json
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
OUT = ROOT / "total_fig467F"
RAW = DATA / "external_fig467F"
PYDEPS = OUT / "pydeps"
sys.path.insert(0, str(PYDEPS))
(OUT / "tables").mkdir(parents=True, exist_ok=True)

UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36"
WGET = [
    "wget", "-c", "--tries=20", "--timeout=60", "--waitretry=8",
    "--retry-connrefused", "--no-verbose", f"--user-agent={UA}",
]


def log(m):
    print(m, flush=True)


def write_tsv(df, name):
    path = OUT / "tables" / name
    df.to_csv(path, sep="\t", index=False)
    log(f"wrote {path.name} nrow={len(df)}")
    return path


def gene_token(col):
    return str(col).split(" ")[0].split("(")[0]


def pick_col(cols, needles):
    cols = list(cols)
    low = [str(c).replace("\n", " ").strip() for c in cols]
    for nd in needles:
        for c, lc in zip(cols, low):
            if lc.lower() == nd.lower():
                return c
    for nd in needles:
        for c, lc in zip(cols, low):
            if nd.lower() in lc.lower():
                return c
    return None


def spearman_ci(x, y, n_boot=1500, seed=1):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 4:
        return dict(n=n, rho=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        r, _ = stats.spearmanr(x[i], y[i])
        if np.isfinite(r):
            boots.append(r)
    lo, hi = np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    return dict(n=n, rho=float(rho), p=float(p), ci_lo=float(lo), ci_hi=float(hi))


def zscore(x):
    x = np.asarray(x, float)
    s = np.nanstd(x)
    if not np.isfinite(s) or s == 0:
        return np.full_like(x, np.nan, dtype=float)
    return (x - np.nanmean(x)) / s


def coxph_fit(time, event, X, maxiter=40):
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    m = np.isfinite(t) & np.isfinite(e) & np.all(np.isfinite(X), axis=1) & (t > 0)
    t, e, X = t[m], e[m], X[m]
    n, p = X.shape
    n_event = int(e.sum())
    if n < 20 or n_event < 8:
        return dict(n=int(n), n_event=n_event, beta=np.full(p, np.nan), se=np.full(p, np.nan))
    order = np.argsort(t, kind="mergesort")
    t, e, X = t[order], e[order], X[order]
    beta = np.zeros(p)
    hess = np.eye(p)
    for _ in range(maxiter):
        eta = np.clip(X @ beta, -20, 20)
        risk = np.exp(eta)
        R = np.cumsum(risk[::-1])[::-1]
        SX = np.cumsum((X * risk[:, None])[::-1], axis=0)[::-1]
        grad = np.zeros(p)
        H = np.zeros((p, p))
        i = 0
        while i < n:
            j = i
            while j < n and t[j] == t[i]:
                j += 1
            d = e[i:j].sum()
            if d > 0 and R[i] > 0:
                mean = SX[i] / R[i]
                wX = X[i:] * risk[i:, None]
                second = (wX.T @ X[i:]) / R[i]
                var = second - np.outer(mean, mean)
                grad += X[i:j][e[i:j] == 1].sum(axis=0) - d * mean
                H += d * var
            i = j
        try:
            step = np.linalg.solve(H + 1e-8 * np.eye(p), grad)
        except np.linalg.LinAlgError:
            break
        beta = beta + step
        hess = H
        if np.max(np.abs(step)) < 1e-8:
            break
    try:
        se = np.sqrt(np.clip(np.diag(np.linalg.inv(hess + 1e-8 * np.eye(p))), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)
    return dict(n=int(n), n_event=n_event, beta=beta, se=se)


def cox_row(endpoint, model, time, event, names, cols, covariates):
    fit = coxph_fit(time, event, cols)
    rows = []
    for i, name in enumerate(names):
        b, se = fit["beta"][i], fit["se"][i]
        hr = np.exp(b) if np.isfinite(b) else np.nan
        z = b / se if (np.isfinite(b) and np.isfinite(se) and se > 0) else np.nan
        p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
        rows.append({
            "endpoint": endpoint, "model": model, "term": name,
            "n": fit["n"], "n_event": fit["n_event"],
            "hr": float(hr) if np.isfinite(hr) else np.nan,
            "ci_lo": float(np.exp(b - 1.96 * se)) if np.isfinite(hr) and np.isfinite(se) else np.nan,
            "ci_hi": float(np.exp(b + 1.96 * se)) if np.isfinite(hr) and np.isfinite(se) else np.nan,
            "p": p, "covariates": covariates,
        })
    return rows


def key6(s):
    return pd.Series(s).astype(str).str.replace("-", "", regex=False).str.extract(r"([A-Z0-9]{5,6})$")[0]


def analyze_cox():
    clin = pd.read_csv(ROOT / "module1/processed/target_clinical_combined.tsv", sep="\t")
    samp = pd.read_csv(ROOT / "module4/processed/target_samples.tsv", sep="\t")
    log("clin " + str(clin.shape) + " samp " + str(samp.shape))
    samp = samp.copy()
    clin = clin.copy()
    samp["_k"] = key6(samp["usi"] if "usi" in samp.columns else samp["case_submitter_id"])
    clin["_k"] = key6(clin["TARGET USI"])
    m = samp.merge(clin, on="_k", how="inner", suffixes=("", "_clin"))
    m = m.drop_duplicates("_k")
    log(f"merged {len(m)} unique USI")

    def os_event(s):
        x = s.astype(str).str.lower()
        return np.where(x.str.contains("dead|deceas|died"), 1,
                        np.where(x.str.contains("alive"), 0, np.nan))

    def efs_event(s):
        x = s.astype(str).str.lower()
        return np.where(x.str.contains("censor|none|no event"), 0,
                        np.where(x.str.contains("relapse|death|dead|fail|progress"), 1, np.nan))

    m["os_time"] = pd.to_numeric(m["Overall Survival Time in Days"], errors="coerce")
    m["os_event"] = os_event(m["Vital Status"])
    m["efs_time"] = pd.to_numeric(m["Event Free Survival Time in Days"], errors="coerce")
    m["efs_event"] = efs_event(m["First Event"])
    m["ZEB1"] = pd.to_numeric(m["ZEB1"], errors="coerce")
    agec = "Age at Diagnosis in Days" if "Age at Diagnosis in Days" in m.columns else "age_at_diagnosis"
    m["age_days"] = pd.to_numeric(m[agec], errors="coerce")
    m["subtype"] = m["subtype"].astype(str) if "subtype" in m.columns else "NA"
    slim = m[["_k", "ZEB1", "os_time", "os_event", "efs_time", "efs_event", "age_days", "subtype"]].copy()
    write_tsv(slim, "F7_target_surv_merged.tsv")
    z = zscore(slim["ZEB1"])
    agez = zscore(slim["age_days"])
    dummies = pd.get_dummies(slim["subtype"].replace({"nan": "NA"}), drop_first=True).astype(float)
    rows = []
    for ep, tt, ee in (("OS", "os_time", "os_event"), ("EFS", "efs_time", "efs_event")):
        rows += cox_row(ep, "ZEB1_unadj", slim[tt], slim[ee], ["ZEB1_z"], z, "ZEB1 (per SD)")
        rows += cox_row(ep, "ZEB1_age", slim[tt], slim[ee], ["ZEB1_z", "age_z"],
                        np.column_stack([z, agez]), "ZEB1 (per SD) + age")
        if dummies.shape[1] >= 1:
            rows += cox_row(ep, "ZEB1_subtype", slim[tt], slim[ee],
                            ["ZEB1_z"] + [f"sub_{c}" for c in dummies.columns],
                            np.column_stack([z, dummies.to_numpy()]),
                            "ZEB1 (per SD) + subtype")
    out = pd.DataFrame(rows)
    write_tsv(out, "F7_target_cox.tsv")
    log(out[out.term.eq("ZEB1_z")].to_string(index=False))
    return out


def load_depmap():
    model = pd.read_csv(DATA / "DepMap_26Q1/Model.csv")
    model["ModelID"] = model["ModelID"].astype(str)
    model["lineage"] = np.where(
        model["OncotreePrimaryDisease"].eq("T-Lymphoblastic Leukemia/Lymphoma"), "T-ALL",
        np.where(model["OncotreePrimaryDisease"].eq("B-Lymphoblastic Leukemia/Lymphoma"), "B-ALL",
                 np.where(model["OncotreeLineage"].eq("Lymphoid"), "other_lymphoid",
                          np.where(model["OncotreeLineage"].eq("Myeloid"), "myeloid", "other"))),
    )
    expr_p = DATA / "DepMap_26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    header = pd.read_csv(expr_p, nrows=0)
    cols = list(header.columns)
    zeb_col = next((c for c in cols if gene_token(c) == "ZEB1"), None)
    idcol = cols[0]
    expr = pd.read_csv(expr_p, usecols=[idcol, zeb_col] if zeb_col else [idcol])
    expr = expr.rename(columns={idcol: "ModelID", zeb_col: "ZEB1"} if zeb_col else {idcol: "ModelID"})
    expr["ModelID"] = expr["ModelID"].astype(str)
    m = model.merge(expr, on="ModelID", how="left")
    m["_key"] = m["StrippedCellLineName"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    log("T-ALL n=" + str(int(m.lineage.eq("T-ALL").sum())))
    return m


def analyze_mtx(dep):
    points, stats_rows = [], []

    def add(src, metric, df, line_col, val_col):
        d = df.copy()
        d["_key"] = d[line_col].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
        merged = d.merge(dep[["ModelID", "CellLineName", "lineage", "ZEB1", "_key"]], on="_key", how="left")
        for _, r in merged.iterrows():
            points.append({
                "source": src, "metric": metric, "drug": "Methotrexate",
                "cell_line": r.get("CellLineName") if pd.notna(r.get("CellLineName")) else r[line_col],
                "model_id": r.get("ModelID"), "lineage": r.get("lineage"),
                "ZEB1": r.get("ZEB1"), "value": r[val_col],
            })
        mp = pd.DataFrame([p for p in points if p["source"] == src])
        for lab, sub in (("T-ALL", mp[mp.lineage.eq("T-ALL")]),
                         ("ALL_TB", mp[mp.lineage.isin(["T-ALL", "B-ALL"])]),
                         ("all_matched", mp.dropna(subset=["ZEB1"]))):
            sp = spearman_ci(sub["ZEB1"], sub["value"])
            stats_rows.append({"source": src, "metric": metric, "subset": lab, **sp})

    for fname, src in (
        ("GDSC2_fitted_dose_response_24Jul22.csv", "GDSC2"),
        ("GDSC1_fitted_dose_response_24Jul22.csv", "GDSC1"),
    ):
        p = RAW / "GDSC" / fname
        if not p.exists():
            continue
        df = pd.read_csv(p, low_memory=False)
        dcol = pick_col(df.columns, ["DRUG_NAME"])
        mtx = df[df[dcol].astype(str).str.contains("methotrexate", case=False, na=False)].copy()
        log(f"{src} MTX rows={len(mtx)}")
        if mtx.empty:
            continue
        line = pick_col(mtx.columns, ["CELL_LINE_NAME"])
        val = pick_col(mtx.columns, ["LN_IC50"])
        mtx = mtx.groupby(line, as_index=False)[val].mean()
        add(src, "LN_IC50", mtx, line, val)
        write_tsv(mtx, f"F7_{src}_mtx_raw.tsv")
    pts, st = pd.DataFrame(points), pd.DataFrame(stats_rows)
    if len(pts):
        write_tsv(pts, "F7_mtx_points.tsv")
    if len(st):
        write_tsv(st, "F7_mtx_stats.tsv")
        log(st.to_string(index=False))
    return pts, st


def wget(url, dest):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 10**7:
        log(f"exists {dest} {dest.stat().st_size}")
        return True
    cmd = WGET + ["-O", str(dest), url]
    log("WGET " + url)
    rc = subprocess.run(cmd).returncode
    log(f"  rc={rc} size={dest.stat().st_size if dest.exists() else 0}")
    return dest.exists() and dest.stat().st_size > 10**7


def analyze_park():
    dest = RAW / "Park_HTA" / "Human_T_cells.h5ad"
    url = "https://datasets.cellxgene.cziscience.com/f9f498e1-11c1-4056-a896-40690435010f.h5ad"
    wget(url, dest)
    if not dest.exists() or dest.stat().st_size < 10**7:
        raise RuntimeError("Human T cells h5ad missing")
    import anndata as ad
    log(f"read {dest} {dest.stat().st_size}")
    adata = ad.read_h5ad(dest, backed="r")
    log(f"n_obs={adata.n_obs} n_vars={adata.n_vars} obs={list(adata.obs.columns)}")
    genes = ["ZEB1", "ZEB2", "LMO2", "TCF7", "GATA3", "BCL11B", "CD3E", "CD1A", "DNTT", "CD34"]
    var_names = np.array(adata.var_names.astype(str))
    symbols = var_names
    for c in ("feature_name", "gene_symbols", "symbol", "Gene"):
        if c in adata.var.columns:
            symbols = np.array(adata.var[c].astype(str))
            break
    want_idx, want_gene = [], []
    for g in genes:
        hit = np.where(var_names == g)[0]
        if len(hit) == 0:
            hit = np.where(symbols == g)[0]
        if len(hit):
            want_idx.append(int(hit[0]))
            want_gene.append(g)
    log("genes " + str(want_gene))
    sub = adata[:, want_idx].to_memory()
    X = sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    expr = pd.DataFrame(np.asarray(X), columns=want_gene, index=sub.obs_names.astype(str))
    obs = sub.obs
    ctype = pick_col(obs.columns, [
        "cell_type", "Anno_level_fig1", "annotation", "celltype", "author_cell_type",
        "cell_type_ontology_term_id",
    ])
    # prefer finer T labels if present
    for c in obs.columns:
        if str(c).lower() in ("anno_level_fig1", "anno_level_fig2", "cell_type", "author_cell_type"):
            ctype = c
            break
    donor = pick_col(obs.columns, ["donor_id", "donor", "sample", "development_stage", "individual"])
    log(f"ctype={ctype} donor={donor}")
    expr["celltype"] = obs[ctype].astype(str).values if ctype else "unknown"
    expr["donor"] = obs[donor].astype(str).values if donor else "unknown"
    agg = expr.groupby(["donor", "celltype"], as_index=False)[want_gene].mean()
    ncell = expr.groupby(["donor", "celltype"]).size().rename("n_cells").reset_index()
    agg = agg.merge(ncell, on=["donor", "celltype"])
    write_tsv(agg, "F4_park_donor_celltype.tsv")
    ct = expr.groupby("celltype", as_index=False).agg(
        n_cells=("donor", "size"), n_donors=("donor", "nunique"),
        **{g: (g, "mean") for g in want_gene},
    )
    write_tsv(ct, "F4_park_celltype_means.tsv")
    log("celltypes\n" + ct.sort_values("n_cells", ascending=False).head(20).to_string(index=False))
    return agg


def main():
    summary = {"started": time.strftime("%Y-%m-%d %H:%M:%S")}
    try:
        analyze_cox()
        summary["cox"] = "ok"
    except Exception as e:
        summary["cox"] = str(e)
        traceback.print_exc()
    try:
        dep = load_depmap()
        analyze_mtx(dep)
        summary["mtx"] = "ok"
    except Exception as e:
        summary["mtx"] = str(e)
        traceback.print_exc()
    try:
        analyze_park()
        summary["park"] = "ok"
    except Exception as e:
        summary["park"] = str(e)
        traceback.print_exc()
    # decoupler tables from first job
    tdir = OUT / "tables"
    summary["tables"] = sorted(p.name for p in tdir.glob("*"))
    summary["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    (tdir / "RUN_SUMMARY_finish.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log("SUMMARY " + json.dumps(summary))
    log("FINISH_DONE")


if __name__ == "__main__":
    main()
