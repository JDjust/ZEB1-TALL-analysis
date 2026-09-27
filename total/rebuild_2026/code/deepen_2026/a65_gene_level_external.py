"""A6.5: gene-level external polarity. No reselection. No drop of reverse genes."""
from pathlib import Path
import gzip
import io
import zlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026")
FR = ROOT / "data/deepen_2026/gate2_program/frozen_ZEB_side50.tsv"
OUT = ROOT / "data/deepen_2026/a65_external"
OUT.mkdir(parents=True, exist_ok=True)
frozen = pd.read_csv(FR, sep="\t")
frozen["symbol_u"] = frozen["symbol"].astype(str).str.upper()
beta = frozen.set_index("symbol_u")["beta_R"].to_dict()
side = frozen.set_index("symbol_u")["side"].to_dict()


def summarize(cohort, ge):
    ge = ge.dropna(subset=["delta_external", "beta_R"]).copy()
    ge["neg_beta"] = -ge["beta_R"]
    ge["expected"] = np.sign(ge["delta_external"]) == np.sign(ge["neg_beta"])
    # zero delta is not expected-direction
    ge.loc[ge["delta_external"] == 0, "expected"] = False
    rho, p = spearmanr(ge["delta_external"], ge["neg_beta"])
    rows = []
    for s, sub in [("all", ge)] + [(k, ge[ge.side == k]) for k in ("ZEB1-side50", "ZEB2-side50")]:
        if s != "all":
            rows.append({
                "cohort": cohort, "set": s, "n": len(sub),
                "expected_frac": float(sub.expected.mean()) if len(sub) else np.nan,
                "rho_delta_vs_neg_beta": np.nan, "rho_p": np.nan,
            })
        else:
            rows.append({
                "cohort": cohort, "set": "all_recovered",
                "n": len(sub),
                "expected_frac": float(sub.expected.mean()),
                "rho_delta_vs_neg_beta": float(rho),
                "rho_p": float(p),
            })
    return ge, pd.DataFrame(rows)


# --- GSE146901 FPKM ---
XLS = Path(r"D:\_bioinformation\ZEB1\total\data\validation\gse146901\raw\GSE146901_TALL_count_FPKM.xls.gz")
ETP = {"076", "077", "093", "097", "098", "107", "108", "115"}
NON = {"102", "103", "116", "117", "118", "121", "122", "123", "124", "132"}
raw = gzip.open(XLS, "rb").read()
df = pd.read_excel(io.BytesIO(raw), sheet_name="FPKM", engine="xlrd")
df.columns = [str(c) for c in df.columns]
symcol = next((c for c in df.columns if c.lower() in {"gene_name", "gene symbol", "symbol", "genename", "gene"}), df.columns[0])
df["symbol_u"] = df[symcol].astype(str).str.upper()
etp_c = [c for c in df.columns if c.replace("_RNAseq", "").replace("_HiC", "") in ETP]
non_c = [c for c in df.columns if c.replace("_RNAseq", "").replace("_HiC", "") in NON]
m = df.set_index("symbol_u")[etp_c + non_c].apply(pd.to_numeric, errors="coerce")
m = m.groupby(level=0).mean()
rows = []
for g in frozen["symbol_u"]:
    if g not in m.index:
        continue
    d = float(m.loc[g, etp_c].mean()) - float(m.loc[g, non_c].mean())
    rows.append({"cohort": "GSE146901", "symbol": g, "side": side[g], "beta_R": beta[g], "delta_external": d})
g1, s1 = summarize("GSE146901", pd.DataFrame(rows))
g1.to_csv(OUT / "a65_GSE146901_genes.tsv", sep="\t", index=False)

# --- GSE234608 recovered counts ---
CNT = ROOT / "data/deepen_2026/a5_gse234608/GSE234608_raw_counts.csv.gz"
LAB = ROOT / "data/deepen_2026/a5_gse234608/a5_sample_labels.tsv"
data = CNT.read_bytes()
try:
    text = gzip.decompress(data).decode("utf-8", "replace")
except (OSError, EOFError):
    text = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data).decode("utf-8", "replace")
if not text.endswith("\n"):
    text = text.rsplit("\n", 1)[0] + "\n"
raw2 = pd.read_csv(io.StringIO(text), header=None)
samples = [str(x).strip().strip('"') for x in raw2.iloc[0, 2:].tolist()]
genes = raw2.iloc[1:, 0].astype(str).str.strip().str.strip('"').str.upper()
mat = raw2.iloc[1:, 2:].apply(pd.to_numeric, errors="coerce")
mat.columns = samples
mat.index = genes.values
mat = mat.groupby(level=0).sum()
cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
logcpm = np.log2(cpm + 1)
lab = pd.read_csv(LAB, sep="\t")
a = lab.loc[lab.group == "ETP_or_MPAL", "sample_id"].astype(str)
b = lab.loc[lab.group == "conventional_TALL", "sample_id"].astype(str)
rows = []
for g in frozen["symbol_u"]:
    if g not in logcpm.index:
        continue
    d = float(logcpm.loc[g, a].mean()) - float(logcpm.loc[g, b].mean())
    rows.append({"cohort": "GSE234608", "symbol": g, "side": side[g], "beta_R": beta[g], "delta_external": d})
g2, s2 = summarize("GSE234608", pd.DataFrame(rows))
g2.to_csv(OUT / "a65_GSE234608_genes.tsv", sep="\t", index=False)

summ = pd.concat([s1, s2], ignore_index=True)
summ.to_csv(OUT / "a65_summary.tsv", sep="\t", index=False)
print(summ.to_string(index=False), flush=True)
both = (float(s1.loc[s1.set == "all_recovered", "expected_frac"].iloc[0]) >= 0.6
        and float(s2.loc[s2.set == "all_recovered", "expected_frac"].iloc[0]) >= 0.6
        and float(s1.loc[s1.set == "all_recovered", "rho_delta_vs_neg_beta"].iloc[0]) > 0
        and float(s2.loc[s2.set == "all_recovered", "rho_delta_vs_neg_beta"].iloc[0]) > 0)
lock = pd.DataFrame([
    {"item": "question", "value": "is external polarity gene-level or a few-gene score artifact?"},
    {"item": "expected", "value": "delta_ETP_minus_conv ~ -beta_R"},
    {"item": "reselection", "value": "none; reverse genes kept"},
    {"item": "gene_level_concordance", "value": "Pass" if both else "Partial" if (
        float(s1.loc[s1.set == "all_recovered", "rho_delta_vs_neg_beta"].iloc[0]) > 0
        and float(s2.loc[s2.set == "all_recovered", "rho_delta_vs_neg_beta"].iloc[0]) > 0
    ) else "Fail"},
    {"item": "upgrade_allowed", "value": "program-score and constituent-gene polarity" if both else "keep program-score wording"},
])
lock.to_csv(OUT / "a65_lock.tsv", sep="\t", index=False)
print(lock.to_string(index=False), flush=True)
