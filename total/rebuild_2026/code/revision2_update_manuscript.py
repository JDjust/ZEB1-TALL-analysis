"""Apply the focused statistical Revision 2 to reviewed manuscript text."""
from pathlib import Path
import re

base=Path(__file__).resolve().parents[1]
man=base/'manuscript'
main=man/'ZEB1_ZEB2_TALL_manuscript.md'
supp=man/'Supplementary_Information.md'

def between(text,start,end,replacement):
    a=text.index(start)
    b=text.index(end,a+len(start))
    return text[:a]+replacement.rstrip()+'\n\n'+text[b:]

m=main.read_text(encoding='utf-8')
abstract=("## Abstract\n\n"
"**Background:** ZEB-family expression distinguishes immature and more differentiated T-lineage leukemia states, but its relationship to normal thymopoiesis and modern molecular subtypes has remained unclear. "
"**Methods:** We compared four normal-human-thymus datasets with diagnostic RNA profiles from 1,309 patients assigned to 17 T-ALL molecular subtypes. A prespecified developmental coordinate used CD34, LYL1 and CD1A, excluding ZEB1 and ZEB2. Primary residuals and subtype effects were estimated within the diagnostic cohort; 3,000 stratified patient bootstraps refitted the developmental model. Independent bulk RNA, malignant-cell pseudobulk, BCL11B-rearranged case RNA and chromatin-contact series provided complementary observations. "
"**Results:** Normal thymopoiesis showed recurrent early-low and cortical/double-positive-high ZEB1–ZEB2 balance with later stage-specific reversals. In T-ALL, developmental position explained 25.6% of balance variation, and subtype added 13.9 percentage points after development (global P = 7.8 × 10⁻⁴⁸). Within-cohort residual medians were −1.76 for BCL11B (bootstrap 95% CI, −2.47 to −0.98), −1.68 for SPI1, −1.50 for LMO2 γδ-like and −1.21 for TME-enriched, compared with −0.23 for ETP-like and +0.64 for TLX3. The joint-model BCL11B adjusted effect was −2.48 (bootstrap 95% CI, −3.08 to −1.91). MLLT10 had a negative residual median but no separately significant joint-model coefficient after 17-subtype correction. All 12 BCL11B-rearranged cases in a lesion-defined series were ZEB2-dominant, including seven without a ZEB2 rearrangement partner. Independent ETP/non-ETP RNA polarity and patient-paired malignant-cell state localization supported the expression axis. HiChIP showed exploratory sample-level contact organization; independent pooled loops and scATAC did not replicate that polarity. "
"**Conclusions:** T-ALL molecular subtypes organize a developmentally patterned ZEB axis in distinct ways. Cohort-internal developmental adjustment distinguishes subtype-associated expression configuration from inherited developmental position without requiring cross-cohort residual scales to be biologically interchangeable.\n")
m=between(m,'## Abstract','## Introduction',abstract)
m=m.replace('We therefore treated the relative ZEB1–ZEB2 expression balance as a descriptive axis and anchored its expected value to normal human thymopoiesis.',
            'We therefore treated relative ZEB1–ZEB2 expression as a descriptive axis, used normal thymopoiesis to establish its developmental shape, and estimated subtype deviations within one diagnostic RNA cohort.')
m=m.replace('Here, we ask whether molecular subtypes differ after normal developmental position is taken into account,',
            'Here, we ask whether molecular subtypes differ after a within-cohort developmental coordinate is taken into account,')
results=("### Developmental residuals separate inherited position from subtype-associated reconfiguration\n\n"
"Normal thymus established a reproducible trajectory shape, but its stage libraries and diagnostic leukemia RNA were standardized separately. We therefore made the 1,309-patient Pölönen cohort the inferential reference: a three-degree-of-freedom developmental spline was fitted to patient balance, and each patient's primary residual was observed minus cohort-expected balance (Fig. 3A–C). The model reproduced the frozen developmental R² of 0.256. In 3,000 patient bootstraps stratified by molecular subtype, the spline was refitted before each subtype median was recomputed (Fig. 3D; Supplementary Table S4).\n\n"
"BCL11B had the most negative within-cohort median residual (−1.76; bootstrap 95% CI, −2.47 to −0.98), followed by SPI1 (−1.68; −2.61 to −1.18), LMO2 γδ-like (−1.50; −2.17 to −0.78) and TME-enriched (−1.21; −2.03 to −0.55). ETP-like was closer to its developmental expectation (−0.23; −0.40 to −0.03), despite its low raw balance, while TLX3 was positive (+0.64; +0.40 to +0.82). MLLT10 also had a negative residual median (−0.75; −1.28 to −0.46), but this descriptive contrast did not translate into a separately significant subtype coefficient in the joint model. These estimates describe expression configuration, not regulatory causality.\n\n"
"To test subtype contribution jointly rather than relying only on one-sample residual summaries, we fitted balance on the same spline plus all 17 subtypes with sum-to-zero contrasts. Subtype increased R² by 0.139 (global nested-model P = 7.8 × 10⁻⁴⁸; Fig. 2F). Relative to the equal-weighted mean subtype intercept, BCL11B had an adjusted effect of −2.48 (bootstrap 95% CI, −3.08 to −1.91; BH FDR = 1.1 × 10⁻¹³), whereas TLX3 was +1.25 (+1.04 to +1.46; BH FDR = 2.3 × 10⁻²⁶; Fig. 3E). The median residual and the joint-model coefficient are distinct estimands, and all 17 estimates and multiplicity-adjusted tests are provided in Supplementary Table S4.\n\n"
"For biological context only, the earlier normal-stage projection retained a similar subtype ordering (Fig. 3F; Supplementary Fig. S2E–H). Its residual magnitude is not an absolute normal-versus-leukemia offset: the cohorts used separate z-score origins and different platforms. In particular, 14 of 18 BCL11B coordinates were outside the observed 20-stage normal range, so the former normal-reference median of −2.64 is reported only as an extrapolated sensitivity value. The primary BCL11B inference is the cohort-internal deviation above.\n")
m=between(m,'### Developmental residuals separate inherited position from subtype-associated reconfiguration',
          '### BCL11B-rearranged leukemias converge',results)
m=m.replace('BCL11B was the clearest negative-residual molecular subtype in the primary T-ALL cohort',
            'BCL11B had the strongest negative within-cohort developmental residual in the primary T-ALL cohort')
clinical=("### Clinical associations largely follow subtype composition\n\n"
"Diagnostic balance showed unadjusted associations with induction failure and MRD that attenuated after molecular-subtype adjustment. The induction-failure odds ratio moved from 0.67 to 0.92, and the MRD ≥0.1% estimate from 0.70 to 0.92. Full six-endpoint estimates, confidence intervals and endpoint-level multiplicity correction are in Supplementary Table S1. These exploratory associations were not developed as an independent prognostic model.\n")
m=between(m,'### Clinical associations largely follow subtype composition','## Discussion',clinical)
discussion=("Immature T-ALL and ETP-associated ZEB2 expression are not new concepts. ETP leukemia was defined by an early developmental program and later connected to specific genomic abnormalities [12,13,15,58]. The advance here is the combination of a ZEB-independent normal developmental scaffold with subtype comparisons estimated entirely within one diagnostic cohort. ETP-like cases have low raw balance but a relatively modest cohort-internal developmental residual. BCL11B, SPI1, LMO2 γδ-like and TME-enriched groups show larger negative deviations; TLX3 departs in the opposite direction. MLLT10 is negative by the residual-median description, although its separate joint-model coefficient is uncertain after subtype correction. The normal-stage projection supports the ordering as a sensitivity only. Its separately standardized scales and the 14 of 18 BCL11B cases outside the measured normal domain preclude reading −2.64 as a biological offset from normal cells.\n\n")
a=m.index('Immature T-ALL and ETP-associated ZEB2 expression are not new concepts.')
b=m.index('BCL11B-rearranged leukemia makes the distinction concrete.',a)
m=m[:a]+discussion+m[b:]
methods=("### ZEB balance and developmental residual\n\n"
"In the primary cohort, balance was z(ZEB1) − z(ZEB2), with each gene standardized across patients after TMM-normalized log2 counts per million. The developmental coordinate was z(CD1A) − [z(CD34) + z(LYL1)]/2. These three genes showed the expected early-to-cortical direction in at least three of four normal resources; ZEB1, ZEB2 and LMO2 were excluded from the coordinate. TMM normalization and digital expression methods follow established work [47–49].\n\n"
"For primary inference, we fitted balance ~ ns(developmental coordinate, df = 3) in the 1,309-patient Pölönen cohort and calculated each patient's observed-minus-fitted residual. We then fitted balance ~ ns(developmental coordinate, df = 3) + molecular subtype, with sum-to-zero subtype contrasts. A subtype coefficient is therefore relative to the equal-weighted mean of the 17 subtype-specific intercepts at the same coordinate. The nested-model F test assessed the global subtype increment; individual adjusted effects used two-sided HC3 robust t tests followed by Benjamini–Hochberg correction across 17 subtypes. Residual medians and model coefficients are reported as different estimands.\n\n"
"For uncertainty, 3,000 patient bootstraps resampled with replacement within each subtype, preserving its observed size. The developmental-only and combined models were refitted in every replicate; percentile 2.5th and 97.5th quantiles formed the intervals for subtype residual medians and adjusted coefficients. The random seed was 20260926. Rare subtypes consequently retain wide intervals. The refitted developmental model reproduced the frozen patient residuals to numerical tolerance and the original model R² values.\n\n"
"Separately, the normal developmental curve used 20 stage units from GSE195812 and GSE206710 (R² = 0.7282), after within-normal-source standardization. Its projection into leukemia is retained only as a biological-shape and ordering sensitivity. Normal and leukemia z scores have different reference populations; normal-projection residual magnitudes are not absolute cross-cohort biological offsets. Normal-domain extrapolation is marked explicitly, including 14 of 18 BCL11B patients beyond the observed stage-coordinate range.\n")
m=between(m,'### ZEB balance and developmental residual',
          '### Independent expression and single-cell comparisons',methods)
m=m.replace('Supplementary Table S1 gives the unadjusted and molecular-subtype-adjusted estimates, confidence intervals, sample sizes, event counts and P values.',
            'Supplementary Table S1 gives the unadjusted and molecular-subtype-adjusted estimates, confidence intervals, sample sizes, event counts, P values and BH correction across six endpoints within each adjustment level.')
m=m.replace('All reported estimates are traced to the completed analysis tables rather than recalculated from visualizations.',
            'All reported estimates are traced to the completed analysis tables or the prespecified Revision 2 patient bootstrap, never inferred from plot geometry.')

fig3=("**Figure 3. Cohort-internal developmental adjustment resolves subtype-specific ZEB-family configurations.**\n\n"
"**A,** Normal-thymus stage units and their spline show the developmental trajectory shape after within-source standardization; the curve is a biological scaffold, not the primary numerical reference for leukemia residuals. **B,** All 1,309 diagnostic Pölönen patients against the spline fitted within that same RNA cohort; BCL11B, ETP-like and TLX3 are highlighted. **C,** Patient-level observed-minus-cohort-expected residuals across 17 subtypes. Diamonds mark medians and densities are drawn only where n ≥10. **D,** Subtype residual medians with percentile 95% intervals from 3,000 stratified patient bootstraps that refit the spline. **E,** Development-adjusted subtype coefficients from the joint spline-plus-subtype model, relative to the equal-weighted mean subtype intercept; intervals are bootstrap percentiles and filled markers indicate HC3 tests meeting BH FDR <0.05 across 17 subtypes. Coefficients and residual medians are different estimands. **F,** The former normal-stage projection compared with primary within-cohort medians for selected subtypes. These reference scales have separate z-score origins. Fourteen of 18 BCL11B developmental coordinates fall outside the observed normal-stage range, so the normal-projection magnitude is a sensitivity result only.\n")
fig4=("**Figure 4. BCL11B-rearranged leukemias converge on a ZEB2-dominant configuration.**\n\n"
"**A,** BCL11B has the most negative median within-cohort developmental residual among 17 Pölönen T-ALL subtypes (18 diagnostic patients; median −1.76). **B,** Paired ZEB1 and ZEB2 TMM log2 CPM values in those 18 patients. **C,** Case-level log2(ZEB1 CPM/ZEB2 CPM) in the 12-case GSE162280 BCL11B-rearranged series; negative values denote ZEB2 dominance. **D,** Rearrangement-class and reported disease-phenotype composition, including AML and MPAL. **E,** Case-group medians and observed ranges: all five ZEB2-fusion and all seven non-ZEB2-partner cases have higher ZEB2 than ZEB1. **F,** Recorded lesion flags and balance in the primary-cohort BCL11B subtype; flags can overlap. GSE162280 is a lesion-defined, mixed-phenotype series without ordinary T-ALL controls and supports convergence, not causation or disease-level replication.\n")

def legends(text):
    text=between(text,'**Figure 3.','**Figure 4.',fig3)
    text=between(text,'**Figure 4.','**Figure 5.',fig4)
    return text
m=legends(m)
main.write_text(m,encoding='utf-8')
(man/'Figure_3_legend.md').write_text(fig3,encoding='utf-8')
(man/'Figure_4_legend.md').write_text(fig4,encoding='utf-8')
legend_file=man/'Figure_Legends.md'
legend_file.write_text(legends(legend_file.read_text(encoding='utf-8')),encoding='utf-8')

s=supp.read_text(encoding='utf-8')
s=s.replace('The normal spline was not refit for the submitted plots.',
            'The normal spline was retained as a biological-shape and sensitivity reference; it is not the primary inferential baseline.')
s=between(s,'### Patient-level subtype and residual source data',
          '### Lim single-cell pseudobulk',
"### Patient-level subtype and residual source data\n\n"
"Figure S2A–D retains all constituent-gene and balance distributions. The primary residual is observed minus the developmental spline fitted within the 1,309-patient cohort. Figure 3D and Table S4 report 3,000 within-subtype patient-bootstrap percentile intervals after refitting that spline in every replicate. Figure 3E and Table S4 additionally report sum-contrast joint-model subtype effects with HC3 P values and 17-test BH FDR. The normal-stage projection is an ordering sensitivity only: normal and leukemia were separately standardized, and 14 of 18 BCL11B patients fall beyond the observed normal-stage coordinate range. Figure S2E–H makes that limitation and both definitions visible without equating their numerical units. Small subtypes retain their individual points and appropriately wide intervals.\n")
s=s.replace('Table S1 contains the frozen six-endpoint model estimates and denominators.',
            'Table S1 contains the frozen six-endpoint model estimates, denominators and BH correction across the six endpoints within each adjustment level.')
s=s.replace('The extreme negative BCL11B position is preserved by the within-cohort sensitivity, but normal-reference extrapolation remains a limit.',
            'The cohort-internal BCL11B position is primary; normal-reference projection is a separately standardized sensitivity limited by extrapolation.')
s=s.replace('The table reports both residual definitions, normal-reference range counts, and frozen P/FDR values. Its statistical labels should be read with the effect magnitude and reference-domain caveat described above.',
            'The table reports primary within-cohort residual medians, bootstrap 95% intervals, joint-model adjusted effects, HC3/BH statistics, and the normal-stage projection as a separate sensitivity. The patient-bootstrap script and full replicate table are supplied with the source data.')
supp.write_text(s,encoding='utf-8')
sl=man/'Supplementary_Figure_Legends_consolidated.md'
st=sl.read_text(encoding='utf-8')
st=st.replace('The extreme negative BCL11B position is preserved by the within-cohort sensitivity, but normal-reference extrapolation remains a limit.',
              'The cohort-internal BCL11B position is primary; normal-reference projection is a separately standardized sensitivity limited by extrapolation.')
sl.write_text(st,encoding='utf-8')
print('Revised manuscript, supplementary narrative and figure legends')
