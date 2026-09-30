#!/usr/bin/env python3
"""Go/No-Go pilot: ZEB1 vs a pre-specified ferroptosis/lipid-redox program.

Five questions, tables only. No figures.
A cohort replication, B subtype adjustment, C T-ALL cell lines,
D CRISPR dependency and ferroptosis drugs, E malignant pseudobulk.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
OUT = ROOT / "pilot_ferroptosis"
OUT.mkdir(parents=True, exist_ok=True)

SUBSTRATE = ["ACSL4", "FADS2", "ELOVL5", "LPCAT3"]
DEFENSE = ["GPX4", "SLC7A11", "AIFM2", "GCLC", "GCLM", "GSS"]
IRON = ["TFRC", "NCOA4"]
STORE = ["FTH1", "FTL"]
OTHER = ["DHODH", "NFE2L2"]
GENES = ["ZEB1"] + SUBSTRATE + DEFENSE + IRON + STORE + OTHER
ARMS = {"substrate": SUBSTRATE, "defense": DEFENSE, "iron": IRON, "store": STORE}
ENSG = {
    "ZEB1": "ENSG00000148516", "ACSL4": "ENSG00000068366", "FADS2": "ENSG00000134824",
    "ELOVL5": "ENSG00000012660", "LPCAT3": "ENSG00000111684", "GPX4": "ENSG00000167468",
    "SLC7A11": "ENSG00000151012", "AIFM2": "ENSG00000042286", "GCLC": "ENSG00000001084",
    "GCLM": "ENSG00000023909", "GSS": "ENSG00000100983", "TFRC": "ENSG00000072274",
    "NCOA4": "ENSG00000138293", "FTH1": "ENSG00000167996", "FTL": "ENSG00000087086",
    "DHODH": "ENSG00000102967", "NFE2L2": "ENSG00000116044",
}
REV = {v: k for k, v in ENSG.items()}


def log(msg):
    print(msg, flush=True)


def symbol_of(token):
    t = str(token).strip().strip('"')
    base = t.split(".")[0]
    if base in ENSG or base in GENES:
        return base if base in GENES else REV.get(base, base)
    if base in REV:
        return REV[base]
    head = t.split(" (")[0].split("|")[0]
    if head in GENES:
        return head
    return None


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 5 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
        return n, np.nan, np.nan
    rho, p = spearmanr(x[ok], y[ok])
    return n, float(rho), float(p)


def arm_score(df, genes):
    cols = [g for g in genes if g in df.columns]
    if len(cols) < 2:
        return None, cols
    z = df[cols].apply(lambda s: (s - s.mean()) / s.std(ddof=0), axis=0)
    return z.mean(axis=1), cols


def load_gene_by_sample(path, genes):
    """First column is the gene id. Keep only requested genes. Samples are columns."""
    want = set(genes) | set(ENSG)
    opener = gzip.open if str(path).endswith(".gz") else open
    found = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = [h.strip().strip('"') for h in header[1:]]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sym = symbol_of(parts[0])
            if sym is None or sym not in genes or sym in found:
                continue
            vals = np.array([np.nan if v in {"", "NA", "NaN"} else float(v) for v in parts[1:]], float)
            if len(vals) != len(samples):
                continue
            found[sym] = vals
            if len(found) == len(genes):
                break
    if not found:
        raise RuntimeError(f"no panel genes in {path}")
    mat = pd.DataFrame(found, index=samples)
    log(f"loaded {path.name} samples={mat.shape[0]} genes={sorted(found)}")
    return mat


def load_depmap_wide(path, model_ids, genes):
    ids = set(map(str, model_ids))
    header = list(pd.read_csv(path, nrows=0).columns)
    if "ModelID" in header:
        key = "ModelID"
    else:
        key = header[0]
    keep = [key]
    rename = {} if key == "ModelID" else {key: "ModelID"}
    if "IsDefaultEntryForModel" in header:
        keep.append("IsDefaultEntryForModel")
    seen = set()
    for col in header:
        if col in keep:
            continue
        sym = symbol_of(col)
        if sym in genes and sym not in seen:
            keep.append(col)
            rename[col] = sym
            seen.add(sym)
    chunks = []
    for ch in pd.read_csv(path, usecols=keep, chunksize=500):
        ch = ch.rename(columns=rename)
        if "IsDefaultEntryForModel" in ch.columns:
            flag = ch["IsDefaultEntryForModel"]
            ch = ch[flag.isin([True, "True", "true", 1, "1", "Yes", "yes"])]
            ch = ch.drop(columns=["IsDefaultEntryForModel"])
        sub = ch[ch["ModelID"].astype(str).isin(ids)]
        if len(sub):
            chunks.append(sub)
    if not chunks:
        return pd.DataFrame(columns=["ModelID"] + [g for g in genes if g in seen])
    out = pd.concat(chunks, ignore_index=True).drop_duplicates("ModelID")
    log(f"depmap {path.name} rows={len(out)} genes={[c for c in out.columns if c != 'ModelID']}")
    return out


def residualize(y, groups):
    y = np.asarray(y, float)
    g = pd.Series(groups).astype(str)
    dummies = pd.get_dummies(g, drop_first=True)
    X = np.column_stack([np.ones(len(y)), dummies.to_numpy(float)])
    ok = np.isfinite(y)
    beta = np.linalg.lstsq(X[ok], y[ok], rcond=None)[0]
    fitted = X @ beta
    out = y - fitted
    out[~ok] = np.nan
    return out


def gene_rows(cohort, df, adjust=None, adjust_name=""):
    rows = []
    z = df["ZEB1"]
    if adjust is not None:
        z = pd.Series(residualize(z, adjust), index=df.index)
    for g in GENES:
        if g == "ZEB1" or g not in df.columns:
            continue
        y = df[g]
        if adjust is not None:
            y = pd.Series(residualize(y, adjust), index=df.index)
        n, rho, p = spearman(z, y)
        rows.append(dict(cohort=cohort, test=adjust_name or "unadjusted", gene=g, n=n, rho=rho, p=p))
    for arm, members in ARMS.items():
        if adjust is None:
            score, used = arm_score(df, members)
            if score is None:
                continue
            n, rho, p = spearman(df["ZEB1"], score)
        else:
            parts = []
            for g in members:
                if g not in df.columns:
                    continue
                parts.append(pd.Series(residualize(df[g], adjust), index=df.index))
            if len(parts) < 2:
                continue
            score = pd.concat(parts, axis=1)
            score = score.apply(lambda s: (s - np.nanmean(s)) / (np.nanstd(s) if np.nanstd(s) else np.nan), axis=0)
            score = score.mean(axis=1)
            used = [g for g in members if g in df.columns]
            n, rho, p = spearman(z, score)
        rows.append(dict(cohort=cohort, test=adjust_name or "unadjusted", gene="ARM_" + arm,
                         n=n, rho=rho, p=p, genes=",".join(used)))
    return rows


def call_arm(rows, gene):
    sub = [r for r in rows if r["gene"] == gene and r["test"] == "unadjusted" and np.isfinite(r["rho"])]
    pos = [r for r in sub if r["rho"] > 0]
    neg = [r for r in sub if r["rho"] < 0]
    major = pos if len(pos) >= len(neg) else neg
    same = len(major) >= 2 and len(sub) >= 2
    any_sig = any(r["p"] < 0.05 for r in major)
    return same and any_sig, len(sub), len(major), "pos" if major is pos else "neg"


def main():
    rows = []
    notes = []

    target = load_gene_by_sample(
        ROOT / "module4/processed/target_log2tpm_filtered.tsv.gz", GENES)
    lab = pd.read_csv(ROOT / "module3/processed/target_tall_subtype.tsv", sep="\t")
    log("TARGET subtype columns " + ",".join(lab.columns))
    idcol = "usi" if "usi" in lab.columns else lab.columns[0]
    lab[idcol] = lab[idcol].astype(str)
    target = target.copy()
    target.index = target.index.astype(str)
    subcol = "subtype" if "subtype" in lab.columns else None
    m = target.join(lab.set_index(idcol)[[c for c in [subcol, "MolecularSubtype", "molecular_subtype"] if c and c in lab.columns]], how="inner")
    # prefer a column that actually separates TLX/TAL/HOXA
    use = None
    for c in m.columns:
        if c in GENES:
            continue
        vals = set(m[c].dropna().astype(str))
        if {"TLX", "TAL", "HOXA"} & vals or any(v.startswith("TLX") or v.startswith("TAL") for v in vals):
            use = c
            break
    if use is None and subcol:
        use = subcol
    keep = m["ZEB1"].notna()
    if use:
        bad = m[use].isna() | m[use].astype(str).isin(["", "Unknown", "nan"])
        keep = keep & ~bad
    tdf = m.loc[keep, [c for c in GENES if c in m.columns]].copy()
    if tdf.empty:
        tdf = target.dropna(subset=["ZEB1"])
        use = None
        notes.append(f"B_TARGET label join empty; expression n={len(tdf)} unadjusted only")
    rows += gene_rows("TARGET", tdf)
    if use and len(tdf):
        rows += gene_rows("TARGET", tdf, m.loc[keep, use], f"residual_{use}")
        notes.append(f"B_TARGET adjustment={use} levels={m.loc[keep, use].astype(str).value_counts().to_dict()}")
    else:
        notes.append("B_TARGET no subtype column")

    ph = load_gene_by_sample(
        ROOT / "module4/processed/pharmacotype_log2fpkm_filtered.tsv.gz", GENES)
    plab = pd.read_csv(ROOT / "module3/processed/pharmacotype_tall_subtype.tsv", sep="\t")
    log("PHARMA columns " + ",".join(map(str, plab.columns[:12])))
    sid = "sample" if "sample" in plab.columns else plab.columns[0]
    plab[sid] = plab[sid].astype(str)
    ph.index = ph.index.astype(str)
    pm = ph.join(plab.set_index(sid)[[c for c in ["etp", "subtype", "Immunophenotype"] if c in plab.columns]], how="inner")
    pdf = pm[[c for c in GENES if c in pm.columns]].copy()
    if pdf.empty:
        pdf = ph.dropna(subset=["ZEB1"])
        notes.append(f"B_StJude label join empty; expression n={len(pdf)} unadjusted only")
    rows += gene_rows("StJude", pdf)
    if "etp" in pm.columns and len(pm):
        rows += gene_rows("StJude", pdf, pm["etp"], "residual_etp")
        notes.append(f"B_StJude adjustment=etp counts={pm['etp'].astype(str).value_counts().to_dict()} "
                     "TLX/TAL/HOXA labels are not in this table")
    else:
        notes.append("B_StJude no etp column")

    nopho_path = DATA / "GSE272023/GSE272023_MatrixRNAseqVSTcountsFilted.txt.gz"
    if nopho_path.exists():
        nop = load_gene_by_sample(nopho_path, GENES)
        # drop non-tumor if a cimp map marks them; keep all finite ZEB1 samples
        ndf = nop.dropna(subset=["ZEB1"])
        rows += gene_rows("NOPHO", ndf)
        notes.append("B_NOPHO not adjusted: GSE272023 local labels are CIMP, not TLX/TAL/HOXA/ETP")
    else:
        notes.append(f"A_NOPHO missing {nopho_path}")

    # C and D: DepMap T-ALL
    model_csv = DATA / "DepMap_26Q1/Model.csv"
    expr_csv = DATA / "DepMap_26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    if not expr_csv.exists():
        alt = DATA / "DepMap_26Q1/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
        expr_csv = alt if alt.exists() else expr_csv
    ge_csv = DATA / "DepMap_26Q1/CRISPRGeneEffect.csv"
    if model_csv.exists() and expr_csv.exists():
        models = pd.read_csv(model_csv, low_memory=False)
        disease = models["OncotreePrimaryDisease"].astype(str)
        tall = models.loc[disease.eq("T-Lymphoblastic Leukemia/Lymphoma"), ["ModelID", "CellLineName"]].copy()
        expr = load_depmap_wide(expr_csv, tall["ModelID"], GENES)
        edf = tall.merge(expr, on="ModelID", how="inner")
        edf = edf.dropna(subset=["ZEB1"])
        for c in GENES:
            if c in edf.columns and c != "ZEB1":
                edf[c] = pd.to_numeric(edf[c], errors="coerce")
        edf["ZEB1"] = pd.to_numeric(edf["ZEB1"], errors="coerce")
        rows += gene_rows("DepMap_TALL", edf)
        notes.append(f"C_TALL_lines n={edf['ZEB1'].notna().sum()} of {len(tall)} Oncotree T-ALL models")
        # random-gene null for substrate and defense mean |rho| is not available without the full matrix.
        # Structure call uses sign concordance inside the pre-specified arms.
        if ge_csv.exists():
            dep_genes = ["GPX4", "SLC7A11", "AIFM2", "GCLC", "GCLM", "GSS", "ACSL4", "FADS2", "ELOVL5", "LPCAT3"]
            ge = load_depmap_wide(ge_csv, edf["ModelID"], dep_genes)
            gdf = edf[["ModelID", "ZEB1"]].merge(ge, on="ModelID", how="inner")
            for g in dep_genes:
                if g not in gdf.columns:
                    continue
                n, rho, p = spearman(gdf["ZEB1"], pd.to_numeric(gdf[g], errors="coerce"))
                rows.append(dict(cohort="DepMap_CRISPR", test="GeneEffect_vs_ZEB1_RNA", gene=g, n=n, rho=rho, p=p))
            notes.append(f"D_CRISPR intersect n={gdf['ZEB1'].notna().sum()} "
                         "GeneEffect more negative means more dependent; "
                         "ZEB1-high GPX4 dependency predicts rho<0")
        else:
            notes.append(f"D_CRISPR missing {ge_csv}")
    else:
        notes.append(f"C missing model={model_csv.exists()} expr={expr_csv.exists()}")

    # D drugs: look only in known project folders for ferroptosis compound names
    drug_hits = []
    needles = ("rsl3", "ml210", "ml162", "erastin", "ferrostatin", "ike", "fin56")
    search_roots = [
        DATA / "DepMap_26Q1",
        DATA / "GDSC",
        ROOT / "module7",
        ROOT / "module8",
    ]
    for base in search_roots:
        if not base.exists():
            continue
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [d for d in dirnames if d not in {"logs", "figures", ".git"}]
            if dirpath.count(os.sep) - str(base).count(os.sep) > 2:
                dirnames.clear()
                continue
            for fn in filenames:
                low = fn.lower()
                if not any(k in low for k in ("drug", "prism", "gdsc", "ctrp", "compound", "auc", "ic50")):
                    continue
                drug_hits.append(str(Path(dirpath) / fn))
    notes.append("D_drug_files " + ("; ".join(drug_hits[:30]) if drug_hits else "none named drug/prism/gdsc/ctrp"))
    # scan small text headers for compound names
    found_compounds = []
    for fp in drug_hits[:40]:
        p = Path(fp)
        if p.stat().st_size > 80_000_000:
            continue
        try:
            opener = gzip.open if fp.endswith(".gz") else open
            with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
                blob = fh.read(200000).lower()
            hit = [k for k in needles if k in blob]
            if hit:
                found_compounds.append(f"{p.name}:{','.join(hit)}")
        except Exception as exc:
            notes.append(f"D_drug_read_fail {p.name} {exc}")
    notes.append("D_ferroptosis_compounds " + ("; ".join(found_compounds) if found_compounds else "not found in scanned headers"))

    # E malignant pseudobulk
    h5 = ROOT / "module5/processed/GSE227122_annotated.h5ad"
    if h5.exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(h5, backed="r")
            var = pd.Series(a.var_names.astype(str))
            idx = {}
            for g in GENES:
                hit = np.where(var.map(symbol_of) == g)[0]
                if len(hit):
                    idx[g] = int(hit[0])
            obs = a.obs
            mal = np.ones(a.n_obs, dtype=bool)
            if "malignant" in obs.columns:
                mal = obs["malignant"].astype(str).str.lower().isin(["malignant", "true", "1", "yes"]).to_numpy()
            if "timepoint" in obs.columns:
                mal = mal & obs["timepoint"].astype(str).eq("Dx").to_numpy()
            patient = obs["patient"].astype(str).to_numpy() if "patient" in obs.columns else np.array(["all"] * a.n_obs)
            take = np.where(mal)[0]
            X = a.X[take][:, list(idx.values())]
            if hasattr(X, "toarray"):
                X = X.toarray()
            X = np.asarray(X, float)
            pb = pd.DataFrame(X, columns=list(idx))
            pb["patient"] = patient[take]
            pb = pb.groupby("patient", observed=True).mean(numeric_only=True)
            rows += gene_rows("GSE227122_malignant_Dx", pb)
            notes.append(f"E_GSE227122 malignant Dx patients={pb.shape[0]} genes={list(idx)}")
            a.file.close()
        except Exception as exc:
            notes.append(f"E_GSE227122 failed {type(exc).__name__}: {exc}")
    else:
        notes.append(f"E missing {h5}")

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "pilot_associations.tsv", sep="\t", index=False)
    (OUT / "pilot_notes.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")

    # decisions
    dec = []
    sub_ok, ncoh, nmajor, direction = call_arm(rows, "ARM_substrate")
    def_ok, ncoh_d, nmajor_d, direction_d = call_arm(rows, "ARM_defense")

    def rho_of(cohort, gene, test):
        hit = [r for r in rows if r["cohort"] == cohort and r["gene"] == gene and r["test"] == test]
        return hit[0] if hit else None

    # structure: opposite signs, not both positive
    cohorts = []
    for cohort in ["TARGET", "StJude", "NOPHO"]:
        s = rho_of(cohort, "ARM_substrate", "unadjusted")
        d = rho_of(cohort, "ARM_defense", "unadjusted")
        if s and d and np.isfinite(s["rho"]) and np.isfinite(d["rho"]):
            cohorts.append((cohort, s["rho"], d["rho"], np.sign(s["rho"]) != np.sign(d["rho"])))
    n_opposite = sum(1 for c in cohorts if c[3])
    both_pos = sum(1 for c in cohorts if c[1] > 0 and c[2] > 0)
    a_call = "GO" if (sub_ok or def_ok) and n_opposite >= 2 else "NO-GO"
    if both_pos >= 2 and n_opposite < 2:
        a_reason = "arms replicate but both move with ZEB1; not a substrate-up/defense-down program"
        a_call = "NO-GO"
    elif n_opposite >= 2 and (sub_ok or def_ok):
        a_reason = f"opposite arm signs in {n_opposite} cohorts; substrate {direction} x{nmajor}/{ncoh}; defense {direction_d} x{nmajor_d}/{ncoh_d}"
    else:
        a_reason = f"no two-cohort same-direction arm with p<0.05 and opposite structure; opposite={n_opposite}"
    dec.append(dict(pilot="A", call=a_call, reason=a_reason))

    # B: if unadjusted TARGET arm is the interesting one, does residual keep the sign?
    b_bits = []
    b_fail = False
    b_ok = False
    for gene in ["ARM_substrate", "ARM_defense"]:
        raw = rho_of("TARGET", gene, "unadjusted")
        tests = [r for r in rows if r["cohort"] == "TARGET" and r["gene"] == gene and r["test"].startswith("residual_")]
        if not raw or not tests or not np.isfinite(raw["rho"]):
            continue
        adj = tests[0]
        kept = np.isfinite(adj["rho"]) and np.sign(adj["rho"]) == np.sign(raw["rho"]) and abs(adj["rho"]) >= 0.5 * abs(raw["rho"])
        b_bits.append(f"{gene} raw {raw['rho']:.3f} residual {adj['rho']:.3f} p {adj['p']:.3g} kept={kept}")
        if abs(raw["rho"]) >= 0.15 and raw["p"] < 0.05 and not kept:
            b_fail = True
        if kept and adj["p"] < 0.05:
            b_ok = True
    st = rho_of("StJude", "ARM_substrate", "residual_etp")
    if st:
        raws = rho_of("StJude", "ARM_substrate", "unadjusted")
        b_bits.append(f"StJude substrate raw {raws['rho']:.3f} etp-residual {st['rho']:.3f}")
    if b_fail:
        b_call = "NO-GO"
    elif b_ok:
        b_call = "GO"
    else:
        b_call = "NO-GO"
    dec.append(dict(pilot="B", call=b_call, reason="; ".join(b_bits) if b_bits else "adjustment not computed"))

    cell = [r for r in rows if r["cohort"] == "DepMap_TALL" and r["gene"] in ("ARM_substrate", "ARM_defense", "ACSL4", "FADS2", "ELOVL5", "GPX4", "SLC7A11", "AIFM2")]
    if not cell:
        c_call, c_reason = "INSUFFICIENT", "no DepMap T-ALL expression rows"
    else:
        n = cell[0]["n"]
        signs = {r["gene"]: r["rho"] for r in cell}
        sub_genes = [signs.get(g) for g in SUBSTRATE if g in signs]
        def_genes = [signs.get(g) for g in ["GPX4", "SLC7A11", "AIFM2"] if g in signs]
        sub_same = len(sub_genes) >= 3 and (sum(x > 0 for x in sub_genes) >= 3 or sum(x < 0 for x in sub_genes) >= 3)
        mixed = False
        if sub_genes and def_genes and np.isfinite(sub_genes).all():
            mixed = (np.nanmean(sub_genes) * np.nanmean(def_genes)) < 0
        if n < 8:
            c_call = "INSUFFICIENT"
        elif sub_same and mixed:
            c_call = "GO"
        else:
            c_call = "NO-GO"
        c_reason = f"n={n} substrate_concordant={sub_same} opposite_defense={mixed} " + \
                   " ".join(f"{k}={v:.2f}" for k, v in signs.items() if np.isfinite(v))
    dec.append(dict(pilot="C", call=c_call, reason=c_reason))

    crispr = [r for r in rows if r["cohort"] == "DepMap_CRISPR"]
    if not crispr:
        d_call, d_reason = "INSUFFICIENT", "no CRISPR overlap"
    else:
        n = min(r["n"] for r in crispr)
        gpx = rho_of("DepMap_CRISPR", "GPX4", "GeneEffect_vs_ZEB1_RNA")
        drug = "ferroptosis compounds present" if found_compounds else "no RSL3/ML210/erastin in scanned drug files"
        if n < 8:
            d_call = "INSUFFICIENT"
        elif gpx and np.isfinite(gpx["rho"]) and gpx["rho"] < 0 and gpx["p"] < 0.05:
            d_call = "GO"
        else:
            d_call = "NO-GO"
        d_reason = f"n={n}; {drug}; " + " ".join(f"{r['gene']} rho={r['rho']:.2f} p={r['p']:.3g}" for r in crispr if r["gene"] in ("GPX4", "SLC7A11", "AIFM2", "GCLC"))
    dec.append(dict(pilot="D", call=d_call, reason=d_reason))

    sc = [r for r in rows if r["cohort"] == "GSE227122_malignant_Dx" and str(r["gene"]).startswith("ARM_")]
    if not sc:
        e_call, e_reason = "INSUFFICIENT", "malignant pseudobulk not computed"
    else:
        n = sc[0]["n"]
        txt = " ".join(f"{r['gene']} rho={r['rho']:.2f} p={r['p']:.3g} n={r['n']}" for r in sc)
        # compare sign to TARGET unadjusted
        agree = 0
        tested = 0
        for r in sc:
            raw = rho_of("TARGET", r["gene"], "unadjusted")
            if raw and np.isfinite(r["rho"]) and np.isfinite(raw["rho"]):
                tested += 1
                agree += int(np.sign(r["rho"]) == np.sign(raw["rho"]))
        if n < 8:
            e_call = "INSUFFICIENT"
        elif tested and agree == tested and any(r["p"] < 0.05 for r in sc):
            e_call = "GO"
        else:
            e_call = "NO-GO"
        e_reason = txt + f"; sign_agree_with_TARGET {agree}/{tested}"
    dec.append(dict(pilot="E", call=e_call, reason=e_reason))

    calls = {d["pilot"]: d["call"] for d in dec}
    if calls["A"] == "NO-GO" or calls["B"] == "NO-GO":
        overall = "NO-GO"
    elif calls["A"] == "GO" and calls["B"] == "GO":
        overall = "GO"
    else:
        overall = "HOLD"
    dec.append(dict(pilot="OVERALL", call=overall,
                    reason="stop if A or B is NO-GO; C/D/E insufficient does not rescue a failed patient program"))
    pd.DataFrame(dec).to_csv(OUT / "GO_NOGO.tsv", sep="\t", index=False)
    log("WROTE " + str(OUT / "GO_NOGO.tsv"))
    log(pd.DataFrame(dec).to_string(index=False))


if __name__ == "__main__":
    main()
