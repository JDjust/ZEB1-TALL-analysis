"""Apply the final journal-style Methods edit without changing results or legends."""
from pathlib import Path

path=Path(__file__).resolve().parents[1]/'manuscript/ZEB1_ZEB2_TALL_manuscript.md'
text=path.read_text(encoding='utf-8')
intro_old='Human thymus atlases now resolve progenitor entry, T-lineage commitment, cortical expansion and selection with increasing precision'
intro_new='Human thymus atlases and mechanistic studies now resolve progenitor entry, T-lineage commitment, cortical expansion and selection with increasing precision'
text=text.replace(intro_old,intro_new)

start=text.index('## Methods\n')
end=text.index('## Data and code availability\n',start)
methods='''## Methods

### Study design and datasets

We analyzed public, de-identified transcriptomic and chromatin datasets using prespecified cohort labels and completed normalization and model outputs. Normal thymopoiesis was represented by GSE142522, GSE195812, GSE206710 and the Park/HTA donor-stage resource [3–5]. GSE142522 contributes sorted bulk libraries by reported stage. The GSE195812 FACS-stage libraries contain cells pooled from six donors; those libraries were treated as stage summaries rather than six independent donor replicates. GSE206710 has three donors and is most informative for the DN-to-DP interval. In Park/HTA, 15 donors with at least 20 cells in both DN and CD4 single-positive populations contributed the paired stage contrast. The GSE195812 T-lineage UMAP was taken from the existing embedding and used to display cell localization, not as a donor-level test. Dataset units, stage labels and marker availability are detailed in Supplementary Figure S1 and Supplementary Table S2.

The primary diagnostic cohort comprised 1,309 unique patients with one RNA sample each and 17 Pölönen molecular subtypes [1]. Supplied molecular subtype and immunophenotype labels were retained. GSE146901 contributed RNA from eight ETP and ten non-ETP leukemias and four healthy peripheral T-cell samples [7]. GSE162280 comprised 12 BCL11B-rearranged immature leukemias with AML, MPAL or ETP-ALL presentations and no ordinary T-ALL controls [8]. The Lim resource provided malignant-cell pseudobulks and 41 eligible Day-0 patient samples with paired ZBTB16-positive and negative states [6]. The HiChIP analysis used four non-ETP, three ETP, two CD34-positive and one thymus library; seven scATAC peak sets formed a separate assay comparison [9]. Case-level and assay-level identifiers accompany the figure source data.

### ZEB balance and developmental residual

In the primary cohort, balance was defined as z(ZEB1) − z(ZEB2), with each gene standardized across patients after TMM-normalized log2 counts per million. The developmental coordinate was z(CD1A) − [z(CD34) + z(LYL1)]/2. These three genes showed the expected early-to-cortical direction in at least three of the four normal resources. ZEB1, ZEB2 and LMO2 were excluded from the coordinate. TMM normalization and digital expression methods follow established work [47–49]. Expression tracks standardized within an atlas were used only for display; absolute units were not equated across protocols.

The normal reference was a natural cubic spline of balance on developmental coordinate with three spline degrees of freedom, fitted to 20 stage units from GSE195812 and GSE206710 (R² = 0.7282). Each patient's developmental residual was observed balance minus the predicted normal-stage balance at the same coordinate. Figure 3 and Supplementary Figure S1 mark the measured normal-stage domain; predictions outside it are extrapolations. A separately specified spline fitted within the leukemia cohort provided a sensitivity residual for that boundary. Subtype medians and patient counts were calculated from the existing patient-level residual tables. The reported variance comparisons use three completed fits—developmental spline, subtype and their combination—without reselecting markers, subtypes or score weights.

### Independent expression and single-cell comparisons

For GSE146901, we retained the sample-level ZEB1 and ZEB2 values, relative balance, log2 expression ratio and original two-sided Mann–Whitney tests. The four healthy peripheral T cells were shown for context and were excluded from both the ETP-versus-non-ETP tests and the normal-thymus reference. Historical TLX and TAL comparisons were restricted to ZEB1, LMO2 and ZEB2, respectively. DerSimonian–Laird and REML/modified Hartung–Knapp estimates and leave-one-cohort-out results are presented together because four cohorts provide limited precision for method-dependent meta-analysis [51,52]. Historical classes were not recoded as the current 17 molecular subtypes.

For the Lim data, malignant cells were aggregated by biological sample and ZBTB16 state before balance comparison. Eligibility required at least 50 malignant cells for a whole-sample pseudobulk and at least 20 cells in each state for the paired analysis. Forty-one Day-0 samples met the paired-state criterion. Technical duplicate L086 R1/R2 was combined before logCPM, and P058 was removed to avoid overlap with another cohort. The paired Wilcoxon test compared ZBTB16-positive and negative balance within each sample. Induction group and molecular-class displays are descriptive patient-level summaries; no outcome model was fitted to these pseudobulks. Twelve eligible Day-0/Day-28 pairs were retained for the longitudinal sensitivity (paired P = 0.569). Patient aggregation avoids treating cells from the same leukemia as independent patients [34–37].

### BCL11B lesions and chromatin measurements

The GSE162280 case manifest identifies five ZEB2-fusion, five ARID1B-enhancer and two CDK6-enhancer cases. Gene-level counts per million and log2(ZEB1 CPM/ZEB2 CPM) were displayed for each case. The 12 lesion-defined cases were compared within the series; no ordinary T-ALL control group was available. Lesion flags for the 18 BCL11B-subtype patients in the primary cohort were taken from their separate annotation table.

HiChIP analyses used sample-level filtered contact counts and depth-normalized gene-body and promoter-window summaries. The leukemia comparison included four non-ETP and three ETP libraries; the exact two-sided Mann–Whitney P for their ZEB2/ZEB1 contact ratio was 0.057. A second thymus archive, THY-2124, failed source quality control and was excluded. GSE146901 pooled loop-call rates and seven scATAC peak-set summaries were retained with their original assay units. In the pooled-loop matrix, within-feature z scores determine color only; printed rates remain calls per 1,000 loops. Pooled calls were not treated as independent patient replicates. HiChIP and loop methods follow the assay literature [38–41].

### Clinical models and figure provenance

The existing clinical models covered induction failure, M2/M3 marrow morphology, MRD ≥0.1%, MRD ≥0.01%, event-free survival and overall survival. Binary models report odds ratios per one-unit higher balance; survival models report hazard ratios. Supplementary Table S1 gives the unadjusted and molecular-subtype-adjusted estimates, confidence intervals, sample sizes, event counts and P values. Subtype-adjusted binary models used the previously specified Firth procedure where sparse outcomes or separation required it [50]. No outcome threshold or prognostic score was optimized in this study.

Figures were produced in R with ggplot2, patchwork and, for the immunophenotype-to-subtype alluvial, ggalluvial. Panel-level source tables and manifests identify the original input, display transformation and plotting script. R session information and vector and raster exports accompany the figures. All reported estimates are traced to the completed analysis tables rather than recalculated from visualizations.

'''
text=text[:start]+methods+text[end:]

before,sep,after=text.partition('## Figure legends')
replacements={
    'frozen diagnostic RNA profiles':'diagnostic RNA profiles',
    'The previously frozen developmental coordinate used':'The developmental coordinate was defined as',
    'The frozen balance was calculated':'Balance was calculated',
    'frozen model comparison P':'model comparison P',
    'at its frozen developmental coordinate':'at its developmental coordinate',
    'The normal-stage reference curve was used as fixed; no new discovery fit was used to choose a different developmental marker set':'The normal-stage reference curve was applied without changing the three marker genes',
    'independently frozen within-cohort spline residual':'predefined within-cohort spline residual',
    'frozen Mann–Whitney P':'Mann–Whitney P',
    'frozen cross-cohort direction':'cross-cohort estimates',
    'frozen paired P':'paired P',
    'frozen contact analysis':'contact analysis',
    'frozen exact two-sided':'exact two-sided',
    'frozen induction-failure odds ratio':'induction-failure odds ratio',
}
for old,new in replacements.items():
    before=before.replace(old,new)
text=before+sep+after
path.write_text(text,encoding='utf-8')
print('Journal-style Methods and Results phrasing updated')
