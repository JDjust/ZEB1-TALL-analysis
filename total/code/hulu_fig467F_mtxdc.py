#!/usr/bin/env python3
"""MTX via Module 8 DepMap line table + decoupler CollecTRI on GSE110635."""
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
sys.path.insert(0, str(OUT / "pydeps"))
(OUT / "tables").mkdir(parents=True, exist_ok=True)


def log(m):
    print(m, flush=True)


def write_tsv(df, name):
    path = OUT / "tables" / name
    df.to_csv(path, sep="\t", index=False)
    log(f"wrote {path.name} nrow={len(df)}")
    return path


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


def load_lines():
    p = ROOT / "module8/processed/line_table_26q1.tsv"
    if not p.exists():
        p = ROOT / "module8/tables/M8_8.1_8.5_line_table.tsv"
    m = pd.read_csv(p, sep="\t")
    m["ModelID"] = m["ModelID"].astype(str)
    if "lineage" not in m.columns:
        m["lineage"] = np.where(m.get("is_tall", False), "T-ALL", "other")
    m["COSMICID_n"] = pd.to_numeric(m.get("COSMICID"), errors="coerce")
    m["SangerModelID"] = m["SangerModelID"].astype(str) if "SangerModelID" in m.columns else ""
    m["_key"] = m["StrippedCellLineName"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    log("lines " + str(m.shape) + " T-ALL ZEB1 "
        + str(int(m.loc[m.lineage.eq("T-ALL"), "ZEB1"].notna().sum())))
    return m


def analyze_mtx(dep):
    points, stats_rows = [], []
    for fname, src in (
        ("GDSC2_fitted_dose_response_24Jul22.csv", "GDSC2"),
        ("GDSC1_fitted_dose_response_24Jul22.csv", "GDSC1"),
    ):
        df = pd.read_csv(RAW / "GDSC" / fname, low_memory=False)
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
        dep_c = dep.dropna(subset=["COSMICID_n"]).drop_duplicates("COSMICID_n")
        dep_s = dep[dep["SangerModelID"].ne("") & dep["SangerModelID"].ne("nan")].drop_duplicates("SangerModelID")
        dep_k = dep.drop_duplicates("_key")
        a = g.merge(dep_c[["ModelID", "CellLineName", "lineage", "ZEB1", "COSMICID_n"]],
                    left_on="COSMIC_ID", right_on="COSMICID_n", how="left")
        b = g.merge(dep_s[["ModelID", "CellLineName", "lineage", "ZEB1", "SangerModelID"]],
                    left_on="SANGER_MODEL_ID", right_on="SangerModelID", how="left")
        c = g.merge(dep_k[["ModelID", "CellLineName", "lineage", "ZEB1", "_key"]], on="_key", how="left")
        merged = a.copy()
        for other in (b, c):
            for col in ("ModelID", "CellLineName", "lineage", "ZEB1"):
                merged[col] = merged[col].fillna(other[col])
        log(f"{src} MTX={len(merged)} ZEB1={int(merged.ZEB1.notna().sum())} "
            f"T-ALL={int((merged.lineage=='T-ALL').sum())}")
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
            stats_rows.append({"source": src, "metric": "LN_IC50", "subset": lab,
                               **spearman_ci(sub["ZEB1"], sub["value"])})
    pts, st = pd.DataFrame(points), pd.DataFrame(stats_rows)
    write_tsv(pts, "F7_mtx_points.tsv")
    write_tsv(st, "F7_mtx_stats.tsv")
    log(st.to_string(index=False))


def analyze_decoupler():
    net_p = RAW / "CollecTRI_omnipath.tsv"
    if not net_p.exists() or net_p.stat().st_size < 1000:
        for u in (
            "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1",
            "https://omnipathdb.org/interactions/?datasets=collectri&genesymbols=1",
        ):
            log("WGET " + u)
            subprocess.run(["wget", "-c", "--tries=12", "--timeout=60", "--no-verbose", "-O", str(net_p), u])
            if net_p.exists() and net_p.stat().st_size > 1000:
                break
    if not net_p.exists() or net_p.stat().st_size < 1000:
        raise RuntimeError("no CollecTRI")
    net = pd.read_csv(net_p, sep="\t")
    log("net " + str(net.shape) + " " + str(list(net.columns)[:15]))
    rename = {}
    for col in net.columns:
        cl = col.lower()
        if cl in ("source_genesymbol", "source"):
            rename[col] = "source"
        elif cl in ("target_genesymbol", "target"):
            rename[col] = "target"
    net = net.rename(columns=rename)
    if "weight" not in net.columns:
        inh = net["is_inhibition"].astype(str).isin(["1", "True", "true"]) if "is_inhibition" in net.columns else False
        net["weight"] = np.where(inh, -1.0, 1.0)
    net = net[["source", "target", "weight"]].dropna().drop_duplicates()
    import decoupler as dc
    counts = pd.read_csv(ROOT / "module9/processed/gse110635_symbol_counts.tsv.gz", sep="\t", index_col=0)
    cd = pd.read_csv(ROOT / "module9/processed/gse110635_symbol_coldata.tsv", sep="\t")
    lib = counts.sum(axis=0).replace(0, np.nan)
    mat = np.log1p(counts.divide(lib, axis=1) * 1e6).T
    mat.index = mat.index.astype(str)
    cd = cd.set_index(cd["gsm"].astype(str))
    common = [i for i in mat.index if i in cd.index]
    mat, cd = mat.loc[common], cd.loc[common]
    log("mat " + str(mat.shape) + " " + str(cd["group"].value_counts().to_dict()))
    acts = None
    try:
        out = dc.run_ulm(mat=mat, net=net, source="source", target="target", weight="weight", min_n=5)
        acts = out[0] if isinstance(out, tuple) else out
    except Exception as e:
        log("run_ulm " + str(e))
        out = dc.mt.ulm(mat, net)
        acts = out[0] if isinstance(out, tuple) else out
    if acts.shape[0] != mat.shape[0] and acts.shape[1] == mat.shape[0]:
        acts = acts.T
    grp = cd["group"].astype(str)
    ctrl = grp.str.lower().str.contains("control|scrambled|ntc|neg")
    kd = grp.str.lower().str.contains("sitlx")
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
    write_tsv(out[out.tf.isin(focus)].copy(), "F6_decoupler_ulm_focus.tsv")
    log(out[out.tf.isin(focus)].to_string(index=False))


def main():
    analyze_mtx(load_lines())
    try:
        analyze_decoupler()
    except Exception:
        traceback.print_exc()
    log("MTX_DC_DONE")


if __name__ == "__main__":
    main()
