"""Write a clean AIEOP S9 label table. S9 is the fusion-annotated subset only."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
XLSX = ROOT / "data" / "deepen_2026" / "clintall" / "SupplementaryTables1-10.xlsx"
OUT = ROOT / "data" / "deepen_2026" / "clintall" / "aieop_s9_labels.tsv"


def normalize_subtype(x):
    if pd.isna(x):
        return ""
    compact = "".join(ch for ch in str(x).upper() if ch.isalnum())
    if compact.startswith("TAL1") and "DP" in compact:
        return "TAL1 DP-like"
    if compact.startswith("TAL1"):
        return "TAL1 αβ-like"
    if "LMO2" in compact and ("G" in compact or "Γ" in compact or "γ" in s or "δ" in s):
        return "LMO2 γδ-like"
    if "STAG2" in compact:
        return "STAG2&LMO2"
    return str(x).strip()


def main():
    raw = pd.read_excel(XLSX, sheet_name="SupplementaryTableS9", header=None)
    # Row 2 is the real header in 0-based indexing after title rows.
    header = [str(v).strip() if pd.notna(v) else f"col_{i}" for i, v in enumerate(raw.iloc[2])]
    d = raw.iloc[3:].copy()
    d.columns = header
    d = d[d["Sample"].notna()].copy()
    d["sample_id"] = d["Sample"].astype(str).str.strip()
    d["predicted_subtype"] = d["Predicted Subtype"].map(normalize_subtype)
    d["inferred_subtype"] = d["Infered Subtype"].map(normalize_subtype)
    d["fusion"] = d["FUSION"].astype(str).str.strip()
    d["tasc"] = d["TASC"].map(normalize_subtype)
    keep = d[["sample_id", "predicted_subtype", "inferred_subtype", "fusion", "tasc"]].copy()
    keep["label_source"] = "clinTALL_S9_fusion_annotated_subset"
    keep["n_published_prediction_rows"] = len(keep)
    keep.to_csv(OUT, sep="\t", index=False)
    print(keep.to_string(index=False))
    print("wrote", OUT, "n=", len(keep))
    print(keep.predicted_subtype.value_counts().to_string())


if __name__ == "__main__":
    main()
