"""Apply source-traceable reviewer-round-three prose changes to the English MD."""
from pathlib import Path

root = Path(__file__).resolve().parents[3]
path = root / "total/rebuild_2026/manuscript/ZEB1_ZEB2_TALL_manuscript_polished.md"
s = path.read_text(encoding="utf-8")

def one(old: str, new: str):
    global s
    count = s.count(old)
    if count != 1:
        raise AssertionError(f"Expected one occurrence, got {count}: {old[:85]!r}")
    s = s.replace(old, new, 1)

one(
    "**Artificial Intelligence:** OpenAI ChatGPT, OpenAI Codex and Cursor (Anysphere) assisted with English drafting and revision, R plotting-code and figure-layout editing, and manuscript formatting. The affected material includes the text and data-visualization figures; no generative tool supplied primary observations. [AUTHOR CONFIRMATION: tool versions/access dates, commercial status, exact use in analysis, prompts/functions, verification of numerical results against source tables, and author responsibility for the final content.]",
    "**Artificial Intelligence:** OpenAI ChatGPT and Codex (OpenAI) and Cursor (Anysphere) assisted with language revision, analysis-code drafting and debugging, data-figure layout, and document assembly during manuscript preparation. The numerical analyses were executed on the cited data with documented scripts; these tools did not supply primary observations. The authors retain responsibility for the accuracy of the analysis, figures, references and final text."
)

start = s.index("## Abstract\n\n") + len("## Abstract\n\n")
end = s.index("\n\n---\n\n## Introduction", start)
s = s[:start] + (
    "ZEB-family transcription factors distinguish immature and differentiated T-lineage states, but whether relative ZEB1–ZEB2 expression differs across T-cell acute lymphoblastic leukemia (T-ALL) subtypes beyond developmental expression remains unclear. We integrated four normal-thymus resources with diagnostic RNA profiles from 1,309 patients in 17 molecular subtypes. A CD34/LYL1/CD1A score excluding ZEB genes served as a developmental-expression proxy. Primary residuals were estimated using a cohort-internal natural spline and 3,000 stratified patient bootstraps. Normal thymus showed an early-to-cortical rise in ZEB balance with stage-specific reversals. The developmental proxy explained 25.6% of balance variation; subtype added 13.9 percentage points to the combined model (nested P = 7.8 × 10⁻⁴⁸). BCL11B had the lowest median residual (−1.76; 95% confidence interval, −2.47 to −0.98) and TLX3 the highest (+0.64). Leave-one-subtype-out fitting retained these poles. Under marker omissions and an expanded expression score, both remained among the two most extreme subtypes in their respective directions. Developmentally matched BCL11B versus ETP-like cases differed in residual balance (mean paired difference −1.92; bootstrap 95% confidence interval, −2.98 to −0.67). All 12 lesion-defined BCL11B-rearranged cases were ZEB2-dominant, including seven without a ZEB2 rearrangement partner. Independent bulk data reproduced ETP/non-ETP polarity, and patient-paired pseudobulks localized the balance to ZBTB16-defined malignant states. No independent cohort reproduced the complete 17-subtype residual ranking. These findings support subtype-associated reconfiguration of a developmentally patterned expression axis while leaving developmental identity, mechanism and clinical utility unresolved."
) + s[end:]

one("Critically, neither ZEB gene contributed to the developmental coordinate, preventing circular inference.",
    "Neither ZEB gene contributed to the developmental coordinate, avoiding direct reuse of the outcome genes. The coordinate remains an expression-based proxy, since leukemia drivers may alter its component markers.")
one("Here we show that molecular subtypes systematically reconfigure a developmentally inherited ZEB1–ZEB2 axis, identify BCL11B-rearranged leukemia as a convergent extreme case independent of rearrangement partner, and demonstrate that independent bulk, single-cell, and chromatin measurements localize consistent expression polarity.",
    "We therefore test whether subtype-associated balance differences persist under alternative proxy definitions and when each subtype is held out of baseline fitting. Independent bulk and single-cell resources address narrower expression contrasts; chromatin data are reported as exploratory context.")

one("The resulting normal reference spline, fitted to 20 stage units from GSE195812 and GSE206710, captured 72.8% of stage-unit balance variation (Fig. 3A).",
    "The descriptive normal-reference spline used 20 measured units: eight pooled FACS stage libraries from GSE195812 and 12 sample–stage aggregates nested within three GSE206710 donors. Its in-sample R² was 72.8%; these units are not 20 independent donors (Fig. 3A; Supplementary Table S2).")
one("Immunophenotype alone explained 23.2% of balance variation, developmental position 25.6%, molecular subtype 34.3%, and the combined developmental-plus-subtype model 39.5% (Fig. 2F). Adding subtype after the developmental spline increased explained variance by 13.9 percentage points (nested-model F test, P = 7.8 × 10⁻⁴⁸). Thus, developmental position and molecular identity overlap substantially but neither substitutes fully for the other.",
    "The descriptive raw R² values were 23.2% for author-assigned immunophenotype, 25.6% for the developmental-expression proxy, 34.3% for molecular subtype and 39.5% for their developmental-plus-subtype combination (Fig. 2F). Because these models use different numbers of parameters, Figure 7C compares adjusted R² (23.0%, 25.4%, 33.5%, 38.6%) and repeated five-fold cross-validated R² (22.7%, 25.0%, 32.4%, 37.2%), respectively. Adding subtype after the spline increased raw explained variance by 13.9 percentage points (partial R² = 18.7%; nested-model F test, P = 7.8 × 10⁻⁴⁸). These are overlapping associations, not causal variance components.")
one("Each patient's primary residual—observed minus cohort-expected balance—quantified the deviation from developmental expectation (Fig. 3C).",
    "Each primary residual was observed balance minus the cohort-fitted expectation at the same three-gene expression score (Fig. 3C). It is conditional on a developmental-expression proxy, not a measured distance from a normal thymocyte or a directly observed maturation stage.")
one("The complete 17-subtype estimates appear in Supplementary Table S3.",
    "The complete 17-subtype estimates appear in Supplementary Table S3. The residual median and joint-model subtype coefficient answer different questions: the former summarizes patient residuals from a model without subtype, whereas the latter compares a subtype intercept with the equal-weighted mean of 17 subtype intercepts conditional on the common spline. They can differ in magnitude or sign (for example, TLX1, HOXA9 TCR and NKX2-1) without being inconsistent.")
one("Given BCL11B's extreme developmental residual, we examined the 12-case GSE162280 series of BCL11B-rearranged immature leukemias [8,16].",
    "Given BCL11B's lowest observed median residual, we examined the 12-case GSE162280 series of BCL11B-rearranged immature leukemias [8,16]. The rank itself does not establish a difference from every neighboring subtype: fixed bootstrap comparisons with SPI1 and LMO2 γδ-like span zero (Fig. 7E; Supplementary Table S8D).")
one("### Independent bulk expression validates ETP/non-ETP ZEB polarity",
    "### Independent bulk expression reproduces the narrower ETP/non-ETP ZEB contrast")
one("In the Lim single-cell resource, analyzed at the biological-sample level [6,29,30,31,32], ZBTB16-positive malignant-cell pseudobulks",
    "In the Lim single-cell resource, analyzed at the biological-sample level [6,29,30,31,32], ZBTB16 detection provided an operational within-patient malignant-cell state label with sufficient cells for paired pseudobulks; it was not assumed to measure developmental stage or causal regulation. ZBTB16-positive malignant-cell pseudobulks")

start = s.index("### HiChIP provides orthogonal chromatin support with defined assay boundaries")
end = s.index("### Clinical associations are largely explained", start)
s = s[:start] + (
    "### Residual polarity persists across held-out subtype fits and alternative expression proxies\n\n"
    "Refitting the developmental spline after excluding each tested subtype preserved BCL11B as the lowest median residual (−2.22) and TLX3 as the highest (+0.78); no BCL11B or TLX3 patient lay beyond the training cohort's coordinate range (Fig. 7A; Supplementary Table S8A). Across three leave-one-marker-out scores and an expanded score selected from normal-source marker directions, BCL11B ranked first or second most negative and TLX3 first or second most positive (Fig. 7B; Supplementary Table S8B). Exact ranks and ETP-like residuals varied by score, so the coordinate should not be equated with a driver-independent developmental stage.\n\n"
    "Eighteen BCL11B cases could be matched one-to-one without replacement to ETP-like patients by developmental coordinate (median distance 0.016; maximum 0.090). The mean paired difference in primary residual was −1.92 (3,000-pair bootstrap 95% CI, −2.98 to −0.67; Fig. 7D). The unpaired BCL11B-minus-ETP-like difference in median residual was −1.53 (95% CI, −2.32 to −0.67), whereas comparisons with SPI1 and LMO2 γδ-like spanned zero (Fig. 7E; Supplementary Table S8D). These within-cohort sensitivity analyses strengthen the BCL11B versus developmentally similar ETP-like contrast but do not provide independent external replication or causal matching.\n\n"
    "HiChIP is retained as exploratory chromatin context (Supplementary Figs. S5 and S9). Three ETP and four non-ETP libraries showed complete rank separation for a ZEB2/ZEB1 contact ratio, but the exact two-sided P value was 0.057. Pooled loops and seven scATAC peak sets did not independently reproduce that polarity; chromatin mechanism is therefore not claimed.\n\n"
) + s[end:]
one("fully adjusted sensitivities also retained associations for MRD ≥0.01%, EFS and OS, whereas induction failure did not",
    "fully adjusted sensitivities also retained associations for MRD ≥0.01%, EFS and OS, whereas induction failure did not; rare-subtype Cox coefficient warnings make the survival estimates exploratory")

start = s.index("**Second, developmental position and molecular subtype")
end = s.index("**Third, BCL11B-rearranged leukemia", start)
s = s[:start] + (
    "**Second, a developmental-expression proxy and molecular subtype supply complementary descriptive information.** A combined model explains 39.5% of balance variation in-sample, and its advantage persists with parameter-adjusted and held-out prediction measures. The three-gene score is supported by normal-stage directionality, yet it remains vulnerable to driver-related marker regulation. Leave-one-marker-out and expanded-score analyses preserve the broad BCL11B/TLX3 polarity while changing exact ranks and the ETP-like residual. Cross-fitting excludes the tested subtype from its own expected-balance curve, but it does not convert this cohort-internal analysis into an independent validation.\n\n"
) + s[end:]
one("**Third, BCL11B-rearranged leukemia presents a convergent extreme.**",
    "**Third, BCL11B-associated leukemia shows an especially strong contrast with developmentally matched ETP-like disease.**")
one("Several important boundaries circumscribe these conclusions. The developmental residual framework requires within-cohort standardization; residual magnitudes cannot be treated as absolute distances from normal cells, particularly for BCL11B cases that fall outside the observed normal-stage coordinate range. The independent datasets each strengthen different aspects of the story—GSE146901 validates ETP/non-ETP polarity, the Lim data localize the axis within patient-matched malignant states, and HiChIP supports chromatin-level organization—but none individually recapitulates the complete framework. The HiChIP P = 0.057 represents the discrete minimum for complete rank separation with groups of three and four, and cross-assay support (pooled loops, scATAC) is negative. Clinical associations attenuate with subtype adjustment, while some fully adjusted endpoint associations persist; no calibrated or external prediction analysis establishes clinical utility.",
    "Several boundaries circumscribe these conclusions. A residual conditioned on three expression markers is neither a normal-cell distance nor direct proof of a maturation state; 14 of 18 primary-cohort BCL11B patients fall beyond the measured normal-reference coordinate range. Cohort-internal cross-fitting and alternative scores test sensitivity, whereas none of the external datasets reproduces the complete 17-subtype residual framework. GSE146901 addresses ETP polarity, Lim addresses paired malignant-cell states and the mixed-phenotype BCL11B series addresses within-lesion partner convergence without an ordinary T-ALL control. HiChIP complete rank separation has exact P = 0.057 and is not independently reproduced by pooled loops or scATAC. Clinical associations attenuate with subtype adjustment; sparse Cox coefficients, absent external calibration and no prediction analysis preclude a clinical-utility claim.")

one("The coordinate markers showed the expected early-to-cortical direction in at least three normal resources. Within the 1,309-patient diagnostic cohort, the primary developmental residual was",
    "The coordinate markers showed the expected early-to-cortical direction in at least three normal resources; this is an expression proxy, not a direct stage measure. Within the diagnostic cohort, the primary residual was")
one("In 3,000 subtype-stratified patient bootstraps (seed 20260926), both models were refitted; percentile intervals quantified uncertainty.",
    "In 3,000 subtype-stratified patient bootstraps (seed 20260926), both models were refitted; percentile intervals quantified uncertainty. Post-review cross-fitting, alternative scores, repeated patient-level cross-validation and focused BCL11B contrasts are detailed in the Supplementary Methods.")
one("Figures were produced in R 4.6.0 with ggplot2 and patchwork.",
    "Figures were produced in R 4.6.0 with ggplot2 and patchwork.")

one("| GSE206710 | normal thymus | three-donor scRNA stage aggregate | 3 stage summaries | DN-to-DP shape; reference units |",
    "| GSE206710 | normal thymus | sample-stage aggregate nested within donor | 12 sample-stage aggregates; 3 donors | DN-to-DP shape; reference units |")

legend7_start = s.index("**Figure 7.")
legend7_end = s.index("\n\n", legend7_start)
s = s[:legend7_start] + (
    "**Figure 7. Developmental-expression proxy robustness and evidence boundaries.** "
    "(**A**) For each of 17 subtypes, median primary residual (grey) and median residual from a natural spline trained on the other 16 subtypes (colored); all 1,309 patients appear only in their own held-out evaluation group. "
    "(**B**) BCL11B, ETP-like and TLX3 residuals under the three-gene score, three leave-one-marker-out definitions and an expanded marker score. Ranks are descriptive across 17 subtypes. "
    "(**C**) Adjusted and repeated five-fold cross-validated R² for author immunophenotype, developmental proxy, molecular subtype and combined models. All 1,309 patients enter full fits; the singleton author-IP Other class is retained in every CV training set, yielding 1,308 held-out predictions per repeat. "
    "(**D**) One-to-one developmental-coordinate matching without replacement for 18 BCL11B and 18 ETP-like cases, with primary residuals joined by pair. Printed interval is a 3,000-pair bootstrap for the mean paired difference. "
    "(**E**) Three fixed BCL11B-minus-comparator differences in median primary residual with 3,000 within-group patient-bootstrap 95% intervals; SPI1 and LMO2 γδ-like intervals include zero. "
    "(**F**) Scope of the primary and supporting observations. Expanded-score markers were selected for early/cortical direction in two normal sources, not for association with ZEB balance. These analyses are cohort-internal sensitivities, not independent replication or causal estimates."
) + s[legend7_end:]
s = s.replace("Plotting code and panel layout were AI-assisted; displayed values derive from the cited datasets and figure source tables and require final author verification.", "")
s = s.replace("Plotting code and panel layout were AI-assisted; displayed values derive from frozen patient-level source tables and require final author verification.", "")
s = s.replace("  ", "  ")
path.write_text(s, encoding="utf-8")
print(path)
