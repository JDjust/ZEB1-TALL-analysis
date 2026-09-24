"""Build a source-linked registry for central numerical claims across all main figures.

Rows preserve source values and uncertainty; this is not a certification that
all upstream models, sample independence, or prose claims have been audited.
"""
from pathlib import Path
import hashlib
import json
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]


def main():
    rows = []
    def ingest(name, figure, unit, effect, context, limits, n=None, lo=None, hi=None,
               p=None, q=None, family="Not established in this registry", select=None):
        path = ROOT / "data" / name
        data = pd.read_csv(path, sep="\t")
        if select:
            data = data.loc[select(data)]
        if data.empty or effect not in data:
            raise ValueError(f"Missing result for {figure}: {name}")
        required = [column for column in [n, lo, hi, p, q] if column]
        missing = set(required) - set(data.columns)
        if missing:
            raise ValueError(f"Missing registry fields {missing}: {name}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        for number, (_, row) in enumerate(data.iterrows(), 1):
            cohort = row.get("cohort", row.get("dataset", row.get("source", context)))
            label = next((row[key] for key in ["contrast", "feature", "drug_short", "CellLineName", "gene", "model", "subtype", "Subtype"] if key in row), "")
            n_value = row.get(n, "") if n else ""
            if not n and "n_a" in row and "n_b" in row:
                n_value = f"{int(row.n_a)} versus {int(row.n_b)}"
            payload = {k: None if pd.isna(v) else v for k, v in row.items()}
            rows.append({"claim_id": f"{figure}_{number:03d}", "claim_context": context,
                         "result_label": label, "figure": figure, "dataset": cohort,
                         "independent_unit": unit, "n": n_value,
                         "effect_measure": effect, "effect": row[effect],
                         "ci95_lo": row.get(lo, "") if lo else "",
                         "ci95_hi": row.get(hi, "") if hi else "",
                         "p": row.get(p, "") if p else "",
                         "fdr": row.get(q, "") if q else "",
                         "testing_family": row.get("fdr_family", family),
                         "source_table": f"data/{name}", "source_sha256": digest,
                         "source_row_json": json.dumps(payload, ensure_ascii=False),
                         "limitations": limits,
                         "verification_scope": "Current exported row; upstream design audit separate"})
    ingest("SourceData_Fig1C_cohort_effects.tsv", "F1C", "patient/sample record", "hedges_g",
           "T-ALL versus B-ALL ZEB1 difference", "Lineage comparison does not isolate leukemia-acquired change",
           lo="ci95_lo", hi="ci95_hi", p="wilcox_p", q="fdr")
    ingest("SourceData_Fig2C_maturity_sensitivity.tsv", "F2C", "expression library", "coef_disease",
           "Disease coefficient under different maturity adjustments", "Six normal stage libraries are not six independent donors; platform comparability remains limited",
           n="n_total_libraries", lo="ci95_lo", hi="ci95_hi", p="p")
    ingest("SourceData_Fig3BC_subtype_effects.tsv", "F3BC", "patient/sample record", "hedges_g",
           "Gene expression by TARGET molecular subtype", "Subtype-defined contrasts are observational",
           lo="ci95_lo", hi="ci95_hi", p="wilcox_p", q="fdr")
    ingest("SourceData_Fig3E_meta.tsv", "F3E", "cohort", "hedges_g",
           "TLX versus rest across cohorts", "Four-study heterogeneity estimates are imprecise; non-equivalent immature groups are not pooled",
           n="k", lo="ci95_lo", hi="ci95_hi", select=lambda d: d.contrast.eq("TLX vs rest"))
    ingest("SourceData_Fig3G_subtype_summary.tsv", "F3G", "patient", "Median_ZEB1",
           "TARGET original molecular-subtype ZEB1 medians", "Same 265 patients as grouped analyses; Q1/Q3 in source rows are interquartile ranges, not confidence intervals; user-completed summaries reused",
           n="n", family="Descriptive medians and IQR; no new tests for this panel")
    ingest("SourceData_Fig4B_crosscohort.tsv", "F4B", "patient", "rho",
           "Key-gene association with ZEB1 in three bulk cohorts", "Unadjusted association; FDR families differ; GSE272023 was previously examined and is not newly held out",
           n="n_patients", p="p", q="fdr")
    ingest("validation/lmo2_within_subtype/within_subtype.tsv", "ST1", "patient", "rho",
           "Within-original-subtype LMO2-ZEB1 association", "Exploratory; same patients as cohort analyses; expression-derived subtype labels; no causal inference",
           n="n", lo="ci_low", hi="ci_high", p="p_perm", q="q_bh",
           family="BH across 14 eligible known-label cohort-subtype strata; 9999 permutations",
           select=lambda d: d.inference)
    ingest("validation/lmo2_within_subtype/cohort_interactions.tsv", "ST1interaction", "patient", "wald_chi2",
           "Within-cohort subtype heterogeneity in standardized rank slopes", "HC3 asymptotic omnibus test; small-stratum precision limited; not independent validation",
           n="n", p="p_hc3", q="q_bh", family="BH across three evaluable cohorts",
           select=lambda d: d.p_hc3.notna())
    ingest("SourceData_Fig5C_malignant_assoc.tsv", "F5C", "patient", "rho",
           "GSE248287 author-malignant patient pseudobulk", "Small cohort; sorted/re-pooled cell mixture; null tests do not establish absent mechanism",
           n="n_patients", p="p", q="fdr_bh", family="BH across 12 tested features")
    ingest("validation/three_gene_patient/three_gene_associations.tsv", "S15", "patient", "rho",
           "Fixed GATA3/TCF7/BCL11B score and constituent genes", "Malignant definitions differ; within-cohort scaling; no cohort pooling",
           n="n_patients", lo="bootstrap_ci95_lo", hi="bootstrap_ci95_hi", p="permutation_p",
           q="fdr_bh_4_features", family="BH across score plus 3 genes within cohort")
    ingest("SourceData_Fig6E_perturbation_forest.tsv", "F6E", "see source n_note; library versus biological replicate audit required", "log2FC",
           "ZEB1/Zeb1 perturbation effects by factor and developmental stage", "Contrasts within one experiment are dependent; manipulations and developmental contexts differ",
           n="n_note", lo="ci_lo", hi="ci_hi", q="fdr")
    ingest("validation/gse287751/author_condition_contrasts.tsv", "S12C", "experimental condition pair", "delta_KO_minus_EV",
           "GSE287751 author-labelled KO versus control by timepoint and priming",
           "Descriptive log2(1+CPM) contrast; four contexts are not independent biological replicates")
    ingest("SourceData_Fig7A_drug_spearman.tsv", "F7A", "patient", "rho",
           "ZEB1 versus ex-vivo drug response", "MTX absent; drug-specific missingness; complete tested family retained",
           n="n", p="p", q="fdr", family="BH within ZEB1 across 17 analyzable drugs",
           select=lambda d: d.gene.eq("ZEB1_log"))
    ingest("SourceData_Fig7B_MRD.tsv", "F7B", "patient", "rho",
           "ZEB1 versus MRD at each timepoint", "Nominal association; protocol-dependent timing and 0.005 below-limit substitution documented in clinical_endpoint_audit; not external validation",
           n="n", p="p", family="Two displayed MRD associations; nominal P")
    ingest("SourceData_Fig7C_ZEB1_dependency.tsv", "F7C", "cell line", "ZEB1_GE",
           "DepMap 26Q1 ZEB1 gene effect by T-ALL model", "Individual cell-line effect is not patient clinical targetability",
           select=lambda d: d.ZEB1_GE.notna())
    ingest("SourceData_Fig7F_mtx_stats.tsv", "F7F", "cell line", "rho",
           "ZEB1 versus MTX within each GDSC T-ALL screen", "Eight shared ModelIDs; screens are not independent biological validation; permutation sensitivity reported separately",
           n="n", lo="ci_lo", hi="ci_hi", p="p", select=lambda d: d.subset.eq("T-ALL"))
    ingest("validation/tlx_meta_sensitivity.tsv", "F3E_sensitivity", "cohort", "estimate",
           "TLX meta-analysis: DL and REML/modified Hartung-Knapp",
           "Four cohorts; prediction intervals exploratory; sensitivity methods are not independent evidence",
           n="k", lo="ci95_lo", hi="ci95_hi", p="p", family="Nominal meta-analysis sensitivity")
    ingest("validation/matched_gene_sets/competitive_summary.tsv", "S14EF", "patient", "observed_rho",
           "Three-gene score versus expression-matched background sets",
           "Post hoc competitive test on selected genes; upper-tail P is not a patient-null association test",
           n="n_patients", p="competitive_upper_tail_p", family="Four exploratory cohort/matching comparisons; nominal competitive probabilities")
    ingest("validation/scrna_patient_qc/cell_threshold_sensitivity.tsv", "S15E", "patient", "rho",
           "Fixed-score sensitivity across minimum malignant-cell counts",
           "Overlapping retained patient sets; no threshold chosen; discovery labels heuristic",
           n="n_patients", lo="bootstrap_ci95_lo", hi="bootstrap_ci95_hi", p="permutation_p",
           family="Ten exploratory threshold comparisons, including repeated patient sets; nominal P")
    ingest("validation/harmony_label_sensitivity/score_association_sensitivity.tsv", "S6B_sensitivity", "patient", "rho",
           "Discovery score under extended Harmony and alternate clustering resolutions",
           "Four computational settings in the same ten patients; original heuristic labels are not independent ground truth",
           n="n_patients", lo="bootstrap_ci95_lo", hi="bootstrap_ci95_hi", p="permutation_p",
           family="Exploratory computational sensitivities; nominal P, no setting selected")
    model_path = ROOT / "data/SourceData_Fig7G_MRD_persistence.tsv"
    model = pd.read_csv(model_path, sep="\t")
    outcome = model.late_positive.astype(int)
    assert set(outcome.unique()) == {0, 1} and model["sample"].is_unique
    for feature in ["loo_prob_early_mrd", "loo_prob_early_mrd_plus_zeb1"]:
        ranks = rankdata(model[feature])
        positive = outcome.eq(1).to_numpy()
        npos, nneg = positive.sum(), (~positive).sum()
        auc = (ranks[positive].sum() - npos * (npos + 1) / 2) / (npos * nneg)
        rows.append({"claim_id": f"F7G_{feature}", "claim_context": "Internal leave-one-out MRD persistence discrimination",
                     "result_label": feature, "figure": "F7G", "dataset": "StJude_Pharmacotype",
                     "independent_unit": "patient", "n": len(model), "effect_measure": "LOO AUC", "effect": auc,
                     "source_table": "data/SourceData_Fig7G_MRD_persistence.tsv",
                     "source_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
                     "limitations": "66 early-positive patients; 18 events; source-defined 0.01% threshold; protocol-specific timing; no independent validation",
                     "verification_scope": "AUC recomputed from exported patient predictions"})
    pd.DataFrame(rows).to_csv(ROOT / "claims.tsv", sep="\t", index=False)
    print(f"Central result registry: {len(rows)} rows, all 7 main figures represented")


if __name__ == "__main__":
    main()
