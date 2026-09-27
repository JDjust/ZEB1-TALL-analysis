"""Gate 1 step 1: published clinTALL calls vs fusion-anchored subtypes.

S9 is the fusion-annotated AIEOP subset only. This is classifier QC,
not residual validation and not a 120-patient label set.
Split inferred labels on 'or' before normalizing, so STIL::TAL1 can
anchor both TAL1 DP-like and TAL1 αβ-like.
"""
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
XLSX = ROOT / "data" / "deepen_2026" / "clintall" / "SupplementaryTables1-10.xlsx"
OUT = ROOT / "data" / "deepen_2026" / "gate1_classifier"
OUT.mkdir(parents=True, exist_ok=True)


def normalize_one(x):
    if pd.isna(x) or not str(x).strip():
        return ""
    compact = "".join(ch for ch in str(x).upper() if ch.isalnum())
    if compact.startswith("TAL1") and "DP" in compact:
        return "TAL1 DP-like"
    if compact.startswith("TAL1"):
        return "TAL1 αβ-like"
    if compact.startswith("LMO2"):
        return "LMO2 γδ-like"
    if "STAG2" in compact:
        return "STAG2&LMO2"
    return str(x).strip()


def parts(text):
    if pd.isna(text) or not str(text).strip():
        return []
    return [normalize_one(p) for p in re.split(r"\s+or\s+", str(text), flags=re.I) if p.strip()]


def main():
    raw = pd.read_excel(XLSX, sheet_name="SupplementaryTableS9", header=None)
    header = [str(v).strip() if pd.notna(v) else f"col_{i}" for i, v in enumerate(raw.iloc[2])]
    d = raw.iloc[3:].copy()
    d.columns = header
    d = d[d["Sample"].notna()].copy()
    d["sample_id"] = d["Sample"].astype(str).str.strip()
    d["predicted_subtype"] = d["Predicted Subtype"].map(normalize_one)
    d["inferred_raw"] = d["Infered Subtype"]
    d["fusion"] = d["FUSION"].astype(str).str.strip()
    d["tasc"] = d["TASC"].map(normalize_one)
    rows = []
    for r in d.itertuples(index=False):
        allowed = [p for p in parts(r.inferred_raw) if p]
        predicted = str(r.predicted_subtype).strip()
        exact = len(allowed) == 1 and predicted == allowed[0]
        compatible = predicted in allowed
        singleton = len(allowed) == 1
        if not allowed:
            call = "unanchored"
        elif exact:
            call = "exact"
        elif compatible:
            call = "compatible_ambiguous_anchor"
        else:
            call = "discordant"
        rows.append(dict(
            sample_id=r.sample_id, fusion=r.fusion,
            predicted_subtype=predicted, inferred_subtype=str(r.inferred_raw),
            tasc=r.tasc, n_allowed=len(allowed),
            singleton_anchor=singleton, exact=exact, compatible=compatible,
            call=call
        ))
    out = pd.DataFrame(rows)
    n = len(out)
    summary = pd.DataFrame([
        dict(item="n_fusion_annotated", value=n),
        dict(item="n_exact", value=int(out.exact.sum())),
        dict(item="n_compatible_including_exact", value=int(out.compatible.sum())),
        dict(item="n_discordant", value=int((out.call == "discordant").sum())),
        dict(item="n_ambiguous_anchors", value=int((~out.singleton_anchor).sum())),
        dict(item="concordance_compatible", value=float(out.compatible.mean())),
        dict(item="concordance_exact_among_singletons",
             value=float(out.loc[out.singleton_anchor, "exact"].mean()) if out.singleton_anchor.any() else None),
        dict(item="paper_reported_accuracy", value="27/33 = 0.818"),
        dict(item="label_status", value="classifier-assigned; fusion-anchored QC only"),
        dict(item="full_120_ready", value="FALSE"),
    ])
    out.to_csv(OUT / "s9_classifier_vs_fusion_anchor.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "s9_classifier_qc.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    print(summary.to_string(index=False))
    print("discordant:")
    print(out.loc[out.call == "discordant", ["sample_id", "fusion", "predicted_subtype", "inferred_subtype"]].to_string(index=False))


if __name__ == "__main__":
    main()
