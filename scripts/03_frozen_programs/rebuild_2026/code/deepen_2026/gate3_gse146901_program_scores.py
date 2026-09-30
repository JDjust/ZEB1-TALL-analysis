"""Gate 3B: frozen program scores in GSE146901 ETP vs non-ETP FPKM."""
from pathlib import Path
import gzip
import io
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]
FR = ROOT / "data" / "deepen_2026" / "gate2_program" / "frozen_ZEB_side50.tsv"
XLS = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\data\validation\gse146901\raw\GSE146901_TALL_count_FPKM.xls.gz")
OUT = ROOT / "data" / "deepen_2026" / "gate3_validation"
OUT.mkdir(parents=True, exist_ok=True)

ETP = {"076", "077", "093", "097", "098", "107", "108", "115"}
NON = {"102", "103", "116", "117", "118", "121", "122", "123", "124", "132"}


def main():
    fr = pd.read_csv(FR, sep="\t")
    raw = gzip.open(XLS, "rb").read()
    df = pd.read_excel(io.BytesIO(raw), sheet_name="FPKM", engine="xlrd")
    df.columns = [str(c) for c in df.columns]
    symcol = None
    for key in ("gene_name", "gene symbol", "symbol", "genename", "gene"):
        hits = [c for c in df.columns if c.lower() == key]
        if hits:
            symcol = hits[0]
            break
    if symcol is None:
        symcol = df.columns[0]
    df["symbol"] = df[symcol].astype(str).str.upper()
    sample_cols = []
    for c in df.columns:
        key = c.replace("_RNAseq", "").replace("_HiC", "")
        if key in ETP or key in NON:
            sample_cols.append(c)
    if len(sample_cols) != 18:
        raise ValueError(f"expected 18 leukemia FPKM columns, got {sample_cols}")

    def score(genes):
        hit = df[df.symbol.isin({g.upper() for g in genes})].copy()
        if hit.empty:
            return pd.Series(np.nan, index=sample_cols), 0
        m = hit[sample_cols].apply(pd.to_numeric, errors="coerce")
        z = (m.sub(m.mean(axis=1), axis=0)).div(m.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
        return z.mean(axis=0, skipna=True), int((~z.isna().all(axis=1)).sum())

    s1, n1 = score(fr.loc[fr.side.eq("ZEB1-side50"), "symbol"])
    s2, n2 = score(fr.loc[fr.side.eq("ZEB2-side50"), "symbol"])
    rec = []
    for c in sample_cols:
        key = c.replace("_RNAseq", "").replace("_HiC", "")
        rec.append(dict(sample=key, group="ETP" if key in ETP else "non-ETP",
                        ZEB1_side50=float(s1[c]), ZEB2_side50=float(s2[c])))
    sc = pd.DataFrame(rec)
    sc.to_csv(OUT / "gse146901_program_scores.tsv", sep="\t", index=False)

    def test(col, etp_higher):
        a = sc.loc[sc.group.eq("ETP"), col].to_numpy()
        b = sc.loc[sc.group.eq("non-ETP"), col].to_numpy()
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        delta = float(np.median(a) - np.median(b))
        return dict(contrast=f"ETP_minus_nonETP_{col}", n_ETP=8, n_nonETP=10,
                    n_genes_scored=n1 if "ZEB1" in col else n2,
                    median_difference=delta, mann_whitney_p=float(p),
                    expected="ETP higher" if etp_higher else "ETP lower",
                    direction_matches=bool(delta > 0) if etp_higher else bool(delta < 0))

    tests = pd.DataFrame([test("ZEB2_side50", True), test("ZEB1_side50", False)])
    tests.to_csv(OUT / "gse146901_program_tests.tsv", sep="\t", index=False)
    print(sc.to_string(index=False))
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
