"""A5 pre-Gate: frozen ZEB programs on GSE234608 bulk. No new DEG/Hallmark."""
from pathlib import Path
import gzip
import io
import re
import zlib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026")
CNT = ROOT / "data/deepen_2026/a5_gse234608/GSE234608_raw_counts.csv.gz"
META = ROOT / "data/deepen_2026/a5_gse234608/GSE234608_series_matrix.txt.gz"
FROZEN = ROOT / "data/deepen_2026/gate2_program/frozen_ZEB_side50.tsv"
OUT = ROOT / "data/deepen_2026/a5_gse234608"
OUT.mkdir(parents=True, exist_ok=True)


def parse_quoted(line):
    return re.findall(r'"([^"]*)"', line)


def _read_counts_text(path: Path) -> str:
    data = path.read_bytes()
    try:
        text = gzip.decompress(data).decode("utf-8", "replace")
        complete = True
    except (OSError, EOFError):
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        text = d.decompress(data).decode("utf-8", "replace")
        complete = False
    if not text.endswith("\n"):
        text = text.rsplit("\n", 1)[0] + "\n"
    print("counts_gzip_complete", complete, "n_lines", text.count("\n"), flush=True)
    return text


def main():
    with gzip.open(META, "rt", encoding="utf-8", errors="replace") as f:
        titles, descs = None, None
        for line in f:
            if line.startswith("!Sample_title"):
                titles = parse_quoted(line)
            elif line.startswith("!Sample_description"):
                descs = parse_quoted(line)
    lab = pd.DataFrame({"sample_id": descs, "geo_title": titles})
    lab["class"] = np.where(
        lab.geo_title.str.startswith("ETP"), "ETP",
        np.where(lab.geo_title.str.startswith("MPAL"), "MPAL", "T-ALL"),
    )
    lab["group"] = np.where(lab["class"].isin(["ETP", "MPAL"]), "ETP_or_MPAL", "conventional_TALL")
    lab.to_csv(OUT / "a5_sample_labels.tsv", sep="\t", index=False)

    raw = pd.read_csv(io.StringIO(_read_counts_text(CNT)), header=None)
    # row0: empty, empty, sample ids
    samples = [str(x).strip().strip('"') for x in raw.iloc[0, 2:].tolist()]
    genes = raw.iloc[1:, 0].astype(str).str.strip().str.strip('"')
    ens = raw.iloc[1:, 1].astype(str).str.strip().str.strip('"')
    mat = raw.iloc[1:, 2:].apply(pd.to_numeric, errors="coerce")
    mat.columns = samples
    mat.index = genes.values
    cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
    logcpm = np.log2(cpm + 1)

    frozen = pd.read_csv(FROZEN, sep="\t")
    avail = {}
    for side in ["ZEB1-side50", "ZEB2-side50"]:
        syms = frozen.loc[frozen.side == side, "symbol"].astype(str)
        hit = [s for s in syms if s in logcpm.index]
        avail[side] = hit
        print(side, "recovered", len(hit), "/", len(syms), flush=True)

    def zmean(cols, gene_list):
        sub = logcpm.loc[gene_list, cols]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
        return z.mean(axis=0)

    def gene_z(sym):
        if sym not in logcpm.index:
            return pd.Series(np.nan, index=logcpm.columns)
        x = logcpm.loc[sym]
        s = x.std(ddof=0)
        return (x - x.mean()) / s if s and np.isfinite(s) and s > 0 else x * np.nan

    score = pd.DataFrame({
        "sample_id": logcpm.columns,
        "ZEB1_side50": zmean(logcpm.columns, avail["ZEB1-side50"]).values,
        "ZEB2_side50": zmean(logcpm.columns, avail["ZEB2-side50"]).values,
        "z_ZEB1": gene_z("ZEB1").values,
        "z_ZEB2": gene_z("ZEB2").values,
        "z_CD1A": gene_z("CD1A").values,
        "z_CD34": gene_z("CD34").values,
        "z_LYL1": gene_z("LYL1").values,
    })
    score["balance"] = score["z_ZEB1"] - score["z_ZEB2"]
    score["dev"] = score["z_CD1A"] - (score["z_CD34"] + score["z_LYL1"]) / 2
    score = score.merge(lab, on="sample_id", how="left")
    if score["class"].isna().any():
        raise ValueError("unlabeled samples: " + ",".join(score.loc[score["class"].isna(), "sample_id"]))
    score.to_csv(OUT / "a5_patient_frozen_scores.tsv", sep="\t", index=False)

    a = score.loc[score.group == "ETP_or_MPAL"]
    b = score.loc[score.group == "conventional_TALL"]
    rows = []
    for col, expect in [
        ("ZEB2_side50", "ETP_MPAL_higher"),
        ("ZEB1_side50", "ETP_MPAL_lower"),
        ("balance", "ETP_MPAL_lower"),
    ]:
        ua, ub = a[col].to_numpy(float), b[col].to_numpy(float)
        p = float(mannwhitneyu(ua, ub, alternative="two-sided").pvalue)
        delta = float(np.nanmedian(ua) - np.nanmedian(ub))
        ok = (delta > 0) if "higher" in expect else (delta < 0)
        rows.append(dict(
            metric=col, n_ETP_MPAL=len(a), n_TALL=len(b),
            median_ETP_MPAL=float(np.nanmedian(ua)), median_TALL=float(np.nanmedian(ub)),
            delta=delta, p=p, expected=expect, direction_ok=ok,
        ))
    tests = pd.DataFrame(rows)
    z2_ok, z1_ok = tests.loc[tests.metric == "ZEB2_side50", "direction_ok"].iloc[0], tests.loc[tests.metric == "ZEB1_side50", "direction_ok"].iloc[0]
    z2_p, z1_p = tests.loc[tests.metric == "ZEB2_side50", "p"].iloc[0], tests.loc[tests.metric == "ZEB1_side50", "p"].iloc[0]
    if z2_ok and z1_ok and min(z2_p, z1_p) < 0.05:
        grade = "Pass"
    elif z2_ok and z1_ok:
        grade = "Partial"
    else:
        grade = "Fail"
    tests["pregate"] = grade
    tests.to_csv(OUT / "a5_pregate_tests.tsv", sep="\t", index=False)
    summary = pd.DataFrame([
        {"item": "n_samples", "value": len(score)},
        {"item": "n_TALL", "value": int((score["class"] == "T-ALL").sum())},
        {"item": "n_ETP", "value": int((score["class"] == "ETP").sum())},
        {"item": "n_MPAL", "value": int((score["class"] == "MPAL").sum())},
        {"item": "ZEB1_side50_recovered", "value": len(avail["ZEB1-side50"])},
        {"item": "ZEB2_side50_recovered", "value": len(avail["ZEB2-side50"])},
        {"item": "primary_contrast", "value": "ETP_or_MPAL vs conventional_TALL"},
        {"item": "iTALL_secondary", "value": "not_in_GEO_titles; not invented"},
        {"item": "pregate", "value": grade},
        {"item": "counts_matrix", "value": "truncated_gzip_recovered; drop last incomplete line"},
        {"item": "cite_seq_next", "value": "download GSE234610 only if Pass or Partial"},
    ])
    summary.to_csv(OUT / "a5_pregate_summary.tsv", sep="\t", index=False)
    print(lab["class"].value_counts().to_string())
    print(tests.to_string(index=False))
    print("grade", grade)


if __name__ == "__main__":
    main()
