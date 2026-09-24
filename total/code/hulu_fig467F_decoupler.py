#!/usr/bin/env python3
"""decoupler 2 ULM on GSE110635 with gene-symbol CollecTRI columns."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
RAW = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F")
OUT = ROOT / "total_fig467F"
sys.path.insert(0, str(OUT / "pydeps"))


def log(m):
    print(m, flush=True)


def write_tsv(df, name):
    path = OUT / "tables" / name
    df.to_csv(path, sep="\t", index=False)
    log(f"wrote {path.name} nrow={len(df)}")


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


def main():
    import decoupler as dc
    raw = pd.read_csv(RAW / "CollecTRI_omnipath.tsv", sep="\t")
    inh = raw["is_inhibition"].astype(str).isin(["1", "True", "true"]) if "is_inhibition" in raw.columns else False
    net = pd.DataFrame({
        "source": raw["source_genesymbol"].astype(str),
        "target": raw["target_genesymbol"].astype(str),
        "weight": np.where(inh, -1.0, 1.0),
    })
    net = net[(net.source != "nan") & (net.target != "nan")]
    net = net.groupby(["source", "target"], as_index=False)["weight"].mean()
    log("net " + str(net.shape) + " TFs=" + str(net.source.nunique()))
    counts = pd.read_csv(ROOT / "module9/processed/gse110635_symbol_counts.tsv.gz", sep="\t", index_col=0)
    cd = pd.read_csv(ROOT / "module9/processed/gse110635_symbol_coldata.tsv", sep="\t")
    lib = counts.sum(axis=0).replace(0, np.nan)
    mat = np.log1p(counts.divide(lib, axis=1) * 1e6).T
    mat.index = mat.index.astype(str)
    cd = cd.set_index(cd["gsm"].astype(str))
    common = [i for i in mat.index if i in cd.index]
    mat, cd = mat.loc[common], cd.loc[common]
    log("mat " + str(mat.shape) + " " + str(cd["group"].value_counts().to_dict()))
    log("decoupler " + str(getattr(dc, "__version__", "?")))
    log("mt attrs " + str([x for x in dir(dc.mt) if not x.startswith("_")]))
    res = dc.mt.ulm(mat, net, tmin=5)
    log("ulm type " + str(type(res)))
    if isinstance(res, tuple):
        acts = res[0]
        log("tuple0 " + str(getattr(acts, "shape", type(acts))))
    elif hasattr(res, "obsm"):
        acts = res.obsm["score_ulm"] if "score_ulm" in getattr(res, "obsm", {}) else None
        log("anndata uns/obsm " + str(list(getattr(res, "obsm", {}).keys())))
    else:
        acts = res
    # decoupler 2 returns AnnData with obsm['ulm_estimate'] or similar
    if acts is None or not isinstance(acts, pd.DataFrame):
        if hasattr(res, "obsm"):
            keys = list(res.obsm.keys())
            log("obsm keys " + str(keys))
            acts = pd.DataFrame(res.obsm[keys[0]], index=res.obs_names, columns=getattr(res, "var_names", None))
            # actually scores are often in uns or a DataFrame
        if hasattr(res, "X") and hasattr(res, "obs_names"):
            # maybe result is samples x TFs AnnData
            X = res.X
            if hasattr(X, "toarray"):
                X = X.toarray()
            acts = pd.DataFrame(np.asarray(X), index=np.array(res.obs_names).astype(str),
                                columns=np.array(res.var_names).astype(str))
    if not isinstance(acts, pd.DataFrame):
        raise RuntimeError("could not parse ULM output " + str(type(res)))
    log("acts " + str(acts.shape) + " idx " + str(list(acts.index)[:5]))
    if acts.shape[0] != mat.shape[0] and acts.shape[1] == mat.shape[0]:
        acts = acts.T
    grp = cd["group"].astype(str)
    ctrl_ids = [i for i in cd.index[grp.str.lower().str.contains("control|scrambled|ntc|neg")] if i in acts.index]
    kd_ids = [i for i in cd.index[grp.str.lower().str.contains("sitlx")] if i in acts.index]
    log(f"ctrl={ctrl_ids} kd={kd_ids}")
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
    sub = out[out.tf.isin(focus)].copy()
    write_tsv(sub, "F6_decoupler_ulm_focus.tsv")
    log(sub.to_string(index=False) if len(sub) else "no focus TFs")
    log("DC_OK")


if __name__ == "__main__":
    main()
