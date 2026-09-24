#!/usr/bin/env python3
"""Rematch GDSC MTX via COSMIC/Sanger IDs; recode TARGET EFS; run decoupler+CollecTRI."""
from __future__ import annotations

import subprocess
import sys
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


def log(m):
    print(m, flush=True)


def write_tsv(df, name):
    path = OUT / "tables" / name
    df.to_csv(path, sep="\t", index=False)
    log(f"wrote {path.name} nrow={len(df)}")
    return path


def gene_token(col):
    return str(col).split(" ")[0].split("(")[0]


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
    if n < 20 or n_event < 5:
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
                vmat = second - np.outer(mean, mean)
                grad += X[i:j][e[i:j] == 1].sum(axis=0) - d * mean
                H += d * vmat
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


def fdr_bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    mask = np.isfinite(p)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranked = np.empty(n, dtype=float)
    prev = 1.0
    for i, k in enumerate(order[::-1]):
        rank = n - i
        val = min(prev, pv[k] * n / rank)
        ranked[k] = val
        prev = val
    out[mask] = ranked
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
    log(f"expr idcol={idcol!r} zeb={zeb_col!r}")
    expr = pd.read_csv(expr_p, usecols=[idcol] + ([zeb_col] if zeb_col else []))
    expr = expr.rename(columns={idcol: "ModelID"})
    if zeb_col:
        expr = expr.rename(columns={zeb_col: "ZEB1"})
    expr["ModelID"] = expr["ModelID"].astype(str)
    m = model.merge(expr, on="ModelID", how="left")
    log("T-ALL with ZEB1 " + str(int(m.loc[m.lineage.eq("T-ALL"), "ZEB1"].notna().sum()))
        + " / " + str(int(m.lineage.eq("T-ALL").sum())))
    return m


def analyze_mtx(dep):
    points, stats_rows = [], []
    dep = dep.copy()
    dep["COSMICID_n"] = pd.to_numeric(dep["COSMICID"], errors="coerce")
    dep["SangerModelID"] = dep["SangerModelID"].astype(str)
    dep["_key"] = dep["StrippedCellLineName"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)

    for fname, src in (
        ("GDSC2_fitted_dose_response_24Jul22.csv", "GDSC2"),
        ("GDSC1_fitted_dose_response_24Jul22.csv", "GDSC1"),
    ):
        p = RAW / "GDSC" / fname
        df = pd.read_csv(p, low_memory=False)
        mtx = df[df["DRUG_NAME"].astype(str).str.contains("methotrexate", case=False, na=False)].copy()
        mtx["LN_IC50"] = pd.to_numeric(mtx["LN_IC50"], errors="coerce")
        mtx["COSMIC_ID"] = pd.to_numeric(mtx["COSMIC_ID"], errors="coerce")
        g = mtx.groupby("CELL_LINE_NAME", as_index=False).agg(
            LN_IC50=("LN_IC50", "mean"),
            COSMIC_ID=("COSMIC_ID", "first"),
            SANGER_MODEL_ID=("SANGER_MODEL_ID", "first") if "SANGER_MODEL_ID" in mtx.columns else ("CELL_LINE_NAME", "first"),
        )
        g["SANGER_MODEL_ID"] = g["SANGER_MODEL_ID"].astype(str)
        g["_key"] = g["CELL_LINE_NAME"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
        merged = g.merge(dep[["ModelID", "CellLineName", "lineage", "ZEB1", "COSMICID_n"]],
                         left_on="COSMIC_ID", right_on="COSMICID_n", how="left")
        miss = merged["ZEB1"].isna()
        extra = g.merge(dep[["ModelID", "CellLineName", "lineage", "ZEB1", "SangerModelID"]],
                        left_on="SANGER_MODEL_ID", right_on="SangerModelID", how="left")
        for c in ("ModelID", "CellLineName", "lineage", "ZEB1"):
            merged.loc[miss, c] = extra.loc[miss, c].values
        miss = merged["ZEB1"].isna()
        extra2 = g.merge(dep[["ModelID", "CellLineName", "lineage", "ZEB1", "_key"]], on="_key", how="left")
        for c in ("ModelID", "CellLineName", "lineage", "ZEB1"):
            merged.loc[miss, c] = extra2.loc[miss, c].values
        log(f"{src} MTX lines={len(merged)} ZEB1={int(merged.ZEB1.notna().sum())} "
            f"T-ALL={int(merged.lineage.eq('T-ALL').sum())}")
        for _, r in merged.iterrows():
            points.append({
                "source": src, "metric": "LN_IC50", "drug": "Methotrexate",
                "cell_line": r["CELL_LINE_NAME"], "model_id": r.get("ModelID"),
                "lineage": r.get("lineage"), "ZEB1": r.get("ZEB1"), "value": r["LN_IC50"],
            })
        mp = pd.DataFrame([x for x in points if x["source"] == src])
        for lab, sub in (("T-ALL", mp[mp.lineage.eq("T-ALL")]),
                         ("ALL_TB", mp[mp.lineage.isin(["T-ALL", "B-ALL"])]),
                         ("all_matched", mp.dropna(subset=["ZEB1"]))):
            sp = spearman_ci(sub["ZEB1"], sub["value"])
            stats_rows.append({"source": src, "metric": "LN_IC50", "subset": lab, **sp})
    pts, st = pd.DataFrame(points), pd.DataFrame(stats_rows)
    write_tsv(pts, "F7_mtx_points.tsv")
    write_tsv(st, "F7_mtx_stats.tsv")
    log(st.to_string(index=False))
    return pts, st


def recode_efs():
    slim = pd.read_csv(OUT / "tables" / "F7_target_surv_merged.tsv", sep="\t")
    clin = pd.read_csv(ROOT / "module1/processed/target_clinical_combined.tsv", sep="\t", keep_default_na=False)
    clin["_k"] = clin["TARGET USI"].astype(str).str.replace("-", "", regex=False).str.extract(r"([A-Z0-9]{5,6})$")[0]
    slim = slim.merge(clin[["_k", "First Event", "Vital Status"]], on="_k", how="left", validate="one_to_one")
    fe = slim["First Event"].fillna("")
    event_map = {"None": 0, "Censored": 0, "Relapse": 1, "Death": 1,
                 "Progression": 1, "Second Malignant Neoplasm": 1}
    unknown = set(fe) - set(event_map) - {""}
    if unknown:
        raise ValueError(f"Unrecognized First Event categories: {unknown}")
    efs = fe.map(event_map).to_numpy(dtype=float)
    slim["efs_event"] = efs
    log("First Event " + str(clin["First Event"].value_counts(dropna=False).head(12).to_dict()))
    log("efs_event recoded " + str(pd.Series(efs).value_counts(dropna=False).to_dict()))
    write_tsv(slim.drop(columns=[c for c in ("First Event", "Vital Status") if c in slim.columns], errors="ignore"),
              "F7_target_surv_merged.tsv")
    z = zscore(slim["ZEB1"])
    agez = zscore(slim["age_days"])
    dummies = pd.get_dummies(slim["subtype"].astype(str).replace({"nan": "NA"}), drop_first=True).astype(float)
    rows = []
    for ep, tt, ee in (("OS", "os_time", "os_event"), ("EFS", "efs_time", "efs_event")):
        rows += cox_row(ep, "ZEB1_unadj", slim[tt], slim[ee], ["ZEB1_z"], z, "ZEB1 (per SD)")
        rows += cox_row(ep, "ZEB1_age", slim[tt], slim[ee], ["ZEB1_z", "age_z"],
                        np.column_stack([z, agez]), "ZEB1 (per SD) + age")
        if dummies.shape[1]:
            rows += cox_row(ep, "ZEB1_subtype", slim[tt], slim[ee],
                            ["ZEB1_z"] + [f"sub_{c}" for c in dummies.columns],
                            np.column_stack([z, dummies.to_numpy()]),
                            "ZEB1 (per SD) + subtype")
    out = pd.DataFrame(rows)
    write_tsv(out, "F7_target_cox.tsv")
    log(out[out.term.eq("ZEB1_z")].to_string(index=False))


def analyze_decoupler():
    net_p = RAW / "CollecTRI_omnipath.tsv"
    if not net_p.exists() or net_p.stat().st_size < 1000:
        urls = [
            "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=sources,references",
            "https://omnipathdb.org/interactions/?datasets=collectri&genesymbols=1",
        ]
        for u in urls:
            log("WGET " + u)
            rc = subprocess.run([
                "wget", "-c", "--tries=10", "--timeout=60", "--no-verbose",
                "-O", str(net_p), u,
            ]).returncode
            if rc == 0 and net_p.exists() and net_p.stat().st_size > 1000:
                break
    if not net_p.exists() or net_p.stat().st_size < 1000:
        raise RuntimeError("CollecTRI download failed")
    net = pd.read_csv(net_p, sep="\t")
    log("collectri cols " + str(list(net.columns)[:20]) + " nrow=" + str(len(net)))
    rename = {}
    for c in net.columns:
        cl = c.lower()
        if cl in ("source_genesymbol", "source", "tf"):
            rename[c] = "source"
        elif cl in ("target_genesymbol", "target", "target_genesymbol"):
            rename[c] = "target"
    net = net.rename(columns=rename)
    if "weight" not in net.columns:
        stim = net["is_stimulation"].astype(str).isin(["1", "True", "true"]) if "is_stimulation" in net.columns else False
        inh = net["is_inhibition"].astype(str).isin(["1", "True", "true"]) if "is_inhibition" in net.columns else False
        w = np.ones(len(net))
        if "is_inhibition" in net.columns:
            w = np.where(inh, -1.0, w)
        net["weight"] = w
    net = net[["source", "target", "weight"]].dropna()
    write_tsv(net.head(20), "F6_collectri_preview.tsv")

    import decoupler as dc
    counts = pd.read_csv(ROOT / "module9/processed/gse110635_symbol_counts.tsv.gz", sep="\t", index_col=0)
    cd = pd.read_csv(ROOT / "module9/processed/gse110635_symbol_coldata.tsv", sep="\t")
    lib = counts.sum(axis=0).replace(0, np.nan)
    mat = np.log1p(counts.divide(lib, axis=1) * 1e6).T
    mat.index = mat.index.astype(str)
    cd["_id"] = cd["gsm"].astype(str)
    cd = cd.set_index("_id")
    common = [i for i in mat.index if i in cd.index]
    mat = mat.loc[common]
    cd = cd.loc[common]
    log(f"mat {mat.shape} groups {cd['group'].value_counts().to_dict()}")
    acts = None
    try:
        out = dc.run_ulm(mat=mat, net=net, source="source", target="target", weight="weight", min_n=5)
        acts = out[0] if isinstance(out, tuple) else out
    except Exception as e:
        log("run_ulm fail " + str(e))
        try:
            out = dc.mt.ulm(mat, net)
            acts = out[0] if isinstance(out, tuple) else out
        except Exception as e2:
            log("mt.ulm fail " + str(e2))
            traceback.print_exc()
    if acts is None:
        raise RuntimeError("ULM failed")
    if acts.shape[0] != mat.shape[0] and acts.shape[1] == mat.shape[0]:
        acts = acts.T
    grp = cd["group"].astype(str)
    ctrl = grp.str.lower().str.contains("control|scrambled|ntc|neg")
    kd = grp.str.lower().str.contains("sitlx")
    if int(kd.sum()) == 0:
        kd = ~ctrl
    ctrl_ids = [i for i in cd.index[ctrl] if i in acts.index]
    kd_ids = [i for i in cd.index[kd] if i in acts.index]
    log(f"ctrl={len(ctrl_ids)} kd={len(kd_ids)}")
    rows = []
    for tf in acts.columns:
        a = pd.to_numeric(acts.reindex(ctrl_ids)[tf], errors="coerce").dropna()
        b = pd.to_numeric(acts.reindex(kd_ids)[tf], errors="coerce").dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        t, p = stats.ttest_ind(b, a, equal_var=False)
        rows.append({
            "tf": tf, "mean_ctrl": float(a.mean()), "mean_kd": float(b.mean()),
            "delta_kd_minus_ctrl": float(b.mean() - a.mean()),
            "t": float(t), "p": float(p), "n_ctrl": int(len(a)), "n_kd": int(len(b)),
        })
    out = pd.DataFrame(rows)
    out["fdr"] = fdr_bh(out["p"])
    out = out.sort_values("p")
    write_tsv(out, "F6_decoupler_ulm_all.tsv")
    focus = ["TLX1", "ZEB1", "GATA3", "TCF7", "BCL11B", "MYC", "MEF2C",
             "RUNX1", "TAL1", "LMO2", "LYL1", "LEF1", "NOTCH1", "MYB", "HOXA9"]
    write_tsv(out[out.tf.isin(focus)], "F6_decoupler_ulm_focus.tsv")
    log(out[out.tf.isin(focus)].to_string(index=False))


def main():
    recode_efs()
    dep = load_depmap()
    analyze_mtx(dep)
    try:
        analyze_decoupler()
    except Exception:
        traceback.print_exc()
    log("REMATCH_DONE")


if __name__ == "__main__":
    main()
