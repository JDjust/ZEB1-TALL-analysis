#!/usr/bin/env python3
"""Count paired patient MRD transitions for an exploratory display.

The 0.01 percent cutpoint follows Lee et al., PMID 36604538; this does not define a clinical
intervention threshold or infer missing outcomes.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "modules/module7/tables/M7_r2_core_scores.tsv"
DEST = ROOT / "data/validation/mrd_transitions.tsv"
EARLY = "Day 15 MRD (%)"
LATE = "Day 42 or 46 MRD (%)"
CUTOFF = .01


def main():
    d = pd.read_csv(SOURCE, sep="\t")
    if d["sample"].duplicated().any():
        raise ValueError("sample identifiers are not unique")
    paired = d.dropna(subset=[EARLY, LATE]).copy()
    if (paired[[EARLY, LATE]] < 0).any().any():
        raise ValueError("MRD percentages cannot be negative")
    paired["day15"] = paired[EARLY].ge(CUTOFF).map({True: "positive", False: "below_threshold"})
    paired["later"] = paired[LATE].ge(CUTOFF).map({True: "positive", False: "below_threshold"})
    out = paired.groupby(["day15", "later"]).size().rename("n_patients").reset_index()
    all_pairs = pd.MultiIndex.from_product(
        [["below_threshold", "positive"], ["below_threshold", "positive"]],
        names=["day15", "later"])
    out = out.set_index(["day15", "later"]).reindex(all_pairs, fill_value=0).reset_index()
    out["threshold_percent_provisional"] = CUTOFF
    out["source_patients"] = len(d)
    out["paired_patients"] = len(paired)
    out["unit"] = "patient"
    if len(paired) != 92 or out.n_patients.sum() != 92:
        raise ValueError("source cohort changed; reassess MRD transition design")
    DEST.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(DEST, sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
