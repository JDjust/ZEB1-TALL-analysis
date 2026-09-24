# -*- coding: utf-8 -*-
"""
Assemble total/manuscript/ZEB1_TALL_manuscript.docx

A submission-formatted manuscript (double-spaced, 12pt Times New Roman,
page numbers; no continuous line numbers) with the seven main figures
embedded. Supplementary Figures S1–S15 are generated separately.
"""
import os
import json
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import _style as S

FIGDIR = S.PNG_DIR
OUT = os.path.join(S.TOTAL_DIR, "manuscript", "ZEB1_TALL_manuscript.docx")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(os.path.join(S.DATA_DIR, "author_metadata.json"), encoding="utf-8") as stream:
    AUTHOR_METADATA = json.load(stream)
    AUTHOR = dict(AUTHOR_METADATA["verbatim"])
    AUTHOR.update(AUTHOR_METADATA.get("confirmed_display_overrides", {}))

doc = Document()

# ---- base style ------------------------------------------------------------
normal = doc.styles["Normal"]
normal.font.name = "Times New Roman"
normal.font.size = Pt(12)
normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
pf = normal.paragraph_format
pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
pf.space_after = Pt(0)

sec = doc.sections[0]
sec.left_margin = sec.right_margin = Inches(1)
sec.top_margin = sec.bottom_margin = Inches(1)

# page numbers in footer (no continuous line numbers)
footer = sec.footer
fp = footer.paragraphs[0]
fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = fp.add_run()
fld1 = OxmlElement("w:fldSimple"); fld1.set(qn("w:instr"), "PAGE")
fp._p.append(fld1)


def H(text, size=13, after=6, before=10, italic=False):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    p.paragraph_format.keep_with_next = True
    p.paragraph_format.space_before = Pt(before)
    p.paragraph_format.space_after = Pt(after)
    r = p.add_run(text)
    r.bold = True
    r.italic = italic
    r.font.size = Pt(size)
    return p


def body(text, justify=True, spacing="double", first_indent=True):
    p = doc.add_paragraph()
    fmt = p.paragraph_format
    fmt.widow_control = True
    fmt.line_spacing_rule = WD_LINE_SPACING.DOUBLE if spacing == "double" else WD_LINE_SPACING.ONE_POINT_FIVE
    if justify:
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    if first_indent:
        fmt.first_line_indent = Inches(0.3)
    # allow **bold** markup
    parts = text.split("**")
    for i, seg in enumerate(parts):
        r = p.add_run(seg)
        if i % 2 == 1:
            r.bold = True
    return p


def figure(png, width=6.4):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next = True
    p.paragraph_format.line_spacing = 1.0
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(2)
    p.add_run().add_picture(os.path.join(FIGDIR, png), width=Inches(width))


def legend(title, text):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.05
    p.paragraph_format.keep_together = True
    p.paragraph_format.space_after = Pt(8)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    r = p.add_run(title + " ")
    r.bold = True; r.font.size = Pt(10)
    r2 = p.add_run(text); r2.font.size = Pt(10)


# ===========================================================================
# TITLE PAGE
# ===========================================================================
t = doc.add_paragraph()
t.alignment = WD_ALIGN_PARAGRAPH.CENTER
t.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
r = t.add_run("Context-dependent ZEB1 expression and treatment-response associations "
              "in T-cell acute lymphoblastic leukemia")
r.bold = True; r.font.size = Pt(15)

for txt, it in [
    (AUTHOR["authors"], False),
    (AUTHOR["affiliations"], False),
    (AUTHOR["correspondence"], False),
    ("# Qian Zhou and Qingyun Ni contributed equally to this work.", False),
    ("Running title: Context-dependent ZEB1 expression in T-ALL", True),
]:
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    rr = p.add_run(txt); rr.italic = it; rr.font.size = Pt(12)

p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.LEFT
p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
p.add_run("Keywords: ").bold = True
p.add_run("T-cell acute lymphoblastic leukemia; ZEB1; transcriptional state; "
          "T-cell development; single-cell RNA-seq; drug sensitivity; CRISPR dependency; "
          "reproducibility.")
p = doc.add_paragraph()
p.add_run("Main figures: 7  |  Supplementary figures: 15").italic = True

doc.add_page_break()

# ===========================================================================
# ABSTRACT
# ===========================================================================
H("Abstract", size=14)
body(
    'Developmental identity may shape how ZEB1 relates to leukemic cell state and treatment response in T-cell acute lymphoblastic leukemia (T-ALL). We integrated public bulk and single-cell transcriptomes, chromatin profiles, perturbations, drug screens and clinical endpoints to define these relationships across biological contexts. Across three cohorts, ZEB1 expression was higher in T-ALL than B-ALL (pooled Hedges g = 1.44). Normal thymocyte references showed developmental-stage variation, while TLX-high T-ALL showed directionally higher ZEB1 across four cohorts (pooled g = 0.61; the small-study-adjusted interval included zero). In a 108-patient NOPHO cohort, ZEB1 correlated positively with BCL11B and IL7R and inversely with LMO2. Perturbation responses depended on context: TLX1 knockdown increased ZEB1, and Lmo2-knockout single-cell data showed higher Zeb1 at day 3 under two priming conditions. Patient single-cell associations varied with cohort and malignant-cell selection. Higher ZEB1 expression correlated with higher methotrexate IC50 in GDSC2 (nine T-ALL models; Spearman rho = 0.90) and GDSC1 (ten models; rho = 0.67); eight shared models limit independence. CRISPR effects were heterogeneous, with four of eight T-ALL models crossing the conventional gene-effect reference of -0.5. In patient data, none of 17 analyzable ex-vivo drug correlations met FDR < 0.05; methotrexate was not assayed. Adding diagnosis ZEB1 to mid-induction residual disease did not improve internal prediction of later positivity in an exploratory 66-patient subset. These findings connect ZEB1 to developmental and molecular context and identify a treatment-specific association with methotrexate response in model systems. Clinical predictive value remains unresolved, motivating context-defined interpretation of ZEB1 biology.'
)

# ===========================================================================
# INTRODUCTION
# ===========================================================================
H("Introduction", size=14)
body(
    "T-cell acute lymphoblastic leukemia (T-ALL) is an aggressive malignancy of thymocyte "
    "precursors driven by oncogenic transcription factors (TAL1, TLX1/TLX3, HOXA, LMO1/2) "
    "acting on arrested stages of T-cell development. A recurring theme in the field is that "
    "the developmental identity of the leukemic cell shapes its transcriptional state, drug "
    "response and outcome. Within this framework, the zinc-finger E-box binding homeobox factor "
    "ZEB1 has attracted attention because experimental work proposed that an "
    "LMO2/SAP18/HDAC1 complex represses ZEB1 through promoter histone deacetylation. "
    "That study reported an increased leukemia stem cell-like phenotype and reduced "
    "methotrexate sensitivity following ZEB1 downregulation. In LMO2-expressing Jurkat "
    "cells, ZEB1 re-expression restored the MTX response assessed by EdU incorporation; "
    "the study also reported attenuation of MTX "
    "insensitivity after addition of the HDAC inhibitor trichostatin A [1]. "
    "A separate multiomic study reported a treatment-resistant, "
    "bone-marrow-progenitor-like T-ALL population [3]; a developmental association "
    "alone does not establish a new resistance state."
)
body(
    "Public patient-scale genomics can place these experimental observations in a wider "
    "developmental and molecular context. A transcription factor can mark lineage identity "
    "while also exerting a functional effect within a particular cellular background. "
    "The relevant questions are therefore where ZEB1-associated programs recur, how "
    "perturbation responses vary across models, and whether expression carries information "
    "about a defined treatment or clinical endpoint. Addressing these questions can "
    "identify settings in which a mechanistic hypothesis is most relevant without assuming "
    "that one direction of association must hold across all T-ALL."
)
body(
    "This reanalysis evaluates the proposed model across public datasets while treating "
    "ZEB1 primarily as a continuous variable. Analyses were defined during this reanalysis; "
    "they were not prospectively registered. We report effect sizes, uncertainty, multiple-"
    "testing adjustment where applicable, and patient-level rather than cell-level inference "
    "for single-cell associations [5]. Discovery, sensitivity, and external evidence are "
    "identified separately, including results that disagree across cohorts."
)

# ===========================================================================
# RESULTS
# ===========================================================================
H("Results", size=14)

H("ZEB1 is a T-lineage gene, elevated in T-ALL relative to other leukemias", size=12.5, italic=True)
body(
    "In the MILE cohort (GSE13159; 2,096 samples), ZEB1 was highest in T-ALL among all leukemia "
    "classes (Figure 1A). T-ALL exceeded B-ALL (Hedges g = 1.21), AML (g = 1.87), normal bone "
    "marrow (g = 1.96), CML, MDS and, marginally, CLL (Figure 1B). These lineage comparisons "
    "do not test repression within a matched T-cell background. The elevation replicated in two independent "
    "RNA-seq cohorts (TARGET g = 1.47; St. Jude pharmacotype g = 1.66), with a random-effects "
    "pooled effect of g = 1.44 (95% CI 1.18-1.70; I2 = 80%; Figure 1C) that was robust to "
    "leave-one-cohort-out analysis. ZEB1 and ZEB2 moved in opposite directions between T-ALL and "
    "B-ALL, and single-gene classifiers discriminated T-ALL well (ZEB1 AUC 0.86; ZEB2-low AUC "
    "0.89; Figure 1D,E). Profiling of 29 healthy immune populations (GSE107011) showed that "
    "ZEB1 is intrinsically a T/lymphoid-lineage gene - highest in T-cell subsets, lowest in "
    "monocytes (T versus B +0.94 log2, P = 0.002; myeloid versus B -1.90, P = 2x10-7; "
    "Figure 1F). The T-ALL-versus-B-ALL difference therefore recapitulates normal lineage "
    "identity; this comparison does not isolate leukemia-acquired changes within T-lineage cells."
)
figure("Figure1_disease_specificity.png")
legend("Figure 1.", "ZEB1 is a T-lineage gene stably elevated in T-ALL. (A) Disease-spectrum "
       "ranking of ZEB1 across MILE leukemias. (B) T-ALL versus each comparator (Hedges g, 95% "
       "CI). (C) Three-cohort replication and random-effects meta-analysis (T-ALL vs B-ALL). "
       "(D) ZEB1/ZEB2/LMO2 effect sizes. (E) Single-gene ROC for T-ALL discrimination. "
       "(F) Annotated expression atlas of 14 selected regulatory/lineage genes across "
       "29 healthy immune populations plus a mixed PBMC reference (GSE107011). Donor-celltype means are averaged "
       "equally across available donors, then standardized within each gene across "
       "populations. Columns are grouped by lineage and hierarchically ordered within "
       "lineage using average-linkage Euclidean distances over the displayed genes. "
       "The ZEB1 row is outlined; P, progenitor; O, other/mixed PBMC. Colors saturate "
       "at z = +/-2.5; full values and absolute expression means are in Source Data. "
       "This relative expression display does not establish a causal lineage program.")

H("ZEB1 is broadly expressed across thymopoiesis, without a simple monotonic maturation program", size=12.5, italic=True)
body(
    "We examined normal thymocyte stages to assess the developmental context of ZEB1 expression. "
    "In sorted human thymocytes (GSE142522, CD34+ through "
    "SP8), ZEB1 varied over a range of 1.34 log2 units, with no clear monotonic stage association "
    "(Spearman rho = -0.43, P = 0.40; Figure 2A), unlike ZEB2 and LMO2. FACS-sorted thymic "
    "single cells (GSE195812, 8 developmental libraries) confirmed that ZEB1 peaks at "
    "DN3/immature DP and does not significantly track stage order (rho = -0.60, P = 0.12; "
    "Figure 2B), even though canonical early transcription factors fall and mature T-cell genes "
    "rise along the same axis (Figure 2F). Although adult T-ALL had lower ZEB1 than sorted "
    "thymocytes (g = -1.28), the disease coefficient depended on the maturity score: "
    "unadjusted -0.889 (95% CI -1.491 to -0.288; P = 0.0046), literature-score adjusted "
    "-0.312 (-0.947 to 0.323; P = 0.328), and empirical-score adjusted -0.728 "
    "(-1.366 to -0.090; P = 0.026; Figure 2C). The six normal-stage libraries are not six "
    "independent donor replicates; cross-platform and sorting differences limit causal "
    "interpretation of these models. "
    "The ZEB1-maturation correlation differed between the St. Jude pharmacotype cohort (rho = "
    "0.60) and TARGET (rho = 0.06; Figure 2D). Exploratory nearest-stage mapping distributed "
    "leukemic cells across two predominant normal-reference stages (Figure 2E); these assignments "
    "do not establish the developmental origin of malignant cells. Donor-level means from the Park/HTA "
    "Human T-cell atlas (CELLxGENE; Figure 2G; expanded in Supplementary Figure 2) placed ZEB1 "
    "high in DN/DP thymocytes and lower in SP; an independent three-donor thymus (GSE206710) "
    "confirmed that ZEB1 remains detectable in CD3+ cells and does not rise monotonically from "
    "DN-early to DP-mature."
)
figure("Figure2_development.png")
legend("Figure 2.", "ZEB1 is broadly expressed across thymopoiesis. (A) Bulk sorted-thymocyte trajectories. "
       "(B) ZEB1 by stage in FACS single-cell thymus. (C) T-ALL disease coefficient in "
       "unadjusted and two maturity-score-adjusted OLS models (41 T-ALL samples and six "
       "normal-stage libraries). (D) ZEB1-maturation correlations in 41 adult T-ALL samples, "
       "TARGET and Pharmacotype; cohort differences are retained. (E) Nearest normal stage for adult T-ALL. "
       "(F) Stage-resolved expression of 12 selected genes across eight FACS-defined stages. "
       "Colors show within-gene z scores of the exported stage mean log1p expression "
       "(population standard deviation; color saturation at -2 and +2); they do not compare "
       "absolute abundance between genes. Rows follow the existing stage-order Spearman "
       "correlations, shown at right, and columns retain developmental order. The ZEB1 row "
       "label is highlighted in red. (G) Independent Park/HTA thymus atlas, "
       "donor-level ZEB1 (donors with n_cells >= 20).")

H("Within T-ALL, ZEB1 varies by subtype; evaluated ZEB1 lesions were uncommon", size=12.5, italic=True)
body(
    "Across TARGET subtypes (Figure 3A), ZEB1 was significantly higher in TLX (g = 0.66) and "
    "lower in the pooled LMO-associated group (LMO1/2 plus LMO2_LYL1; g = -0.80), "
    "with no clear differences for TAL, HOXA and NKX2-1 (Figure 3B). Group estimates "
    "differed on a ZEB1-LMO2 plane: TLX was ZEB1-high/LMO2-low while the pooled LMO group was "
    "ZEB1-low/LMO2-high (Figure 3C), supporting subtype-dependent transcriptional contexts. "
    "The original fine labels reveal that LMO2_LYL1 cases had a median ZEB1 of 4.88, "
    "compared with 5.46 in LMO1/2; TLX1 and TLX3 medians were 6.12 and 6.13, "
    "respectively (log2[TPM + 1]; Figure 3G). These are separate descriptions within "
    "the same TARGET cohort. "
    "Testing four cohorts (Figure 3D,E), the TLX-high effect "
    "was directionally consistent (random-effects g = 0.61, 95% CI 0.22-1.00; I2 = 60%). "
    "A small-study sensitivity analysis using REML and modified Hartung-Knapp inference "
    "retained g = 0.61 with a wider 95% CI of -0.09 to 1.31 (P = 0.069); its exploratory "
    "prediction interval was -1.10 to 2.32. Thus, direction was consistent across the four "
    "cohorts while the pooled magnitude and transportability remained imprecise. "
    "The study-defined immature groups use non-equivalent molecular or phenotype labels; "
    "their effects are shown separately and are not treated as one pooled class. In the "
    "pharmacotype cohort, ETP T-ALL was ZEB1-low relative to conventional T-ALL after age "
    "adjustment (Figure 3F). In the evaluated records, no ZEB1 coding mutation or fusion "
    "was annotated among 265 RNA-linked records; variant callability was not established "
    "for every record, so this denominator does not estimate genetic lesion prevalence. "
    "Locus copy number (rho = 0.03) and global CNV burden "
    "(rho = 0.02) had weak rank associations with ZEB1 expression (Supplementary Fig. S3). "
    "These observations do not exclude other "
    "cis-acting mechanisms or rare lesions."
)
figure("Figure3_subtype.png")
legend("Figure 3.", "ZEB1 varies with molecular subtype. (A) TARGET grouped subtype composition. "
       "In A-C, LMO-associated combines LMO1/2 (n=10) and LMO2_LYL1 (n=18); "
       "TAL combines TAL1/TAL2, and TLX combines TLX1/TLX3. "
       "(B) ZEB1 by subtype (Hedges g). (C) ZEB1-LMO2 axis. (D) Cross-cohort ZEB1 effect "
       "matrix; study-defined immature labels are not equivalent. (E) Random-effects TLX "
       "meta-analysis with DerSimonian-Laird and REML/modified Hartung-Knapp intervals; "
       "the latter accounts for small study count. Sensitivity and "
       "leave-one-cohort-out results are provided in data/validation/tlx_meta_*.tsv. "
       "(F) ETP comparison after age adjustment. (G) Original nine TARGET subtype "
       "labels, shown separately using all 265 patient values. Horizontal black marks "
       "are medians and vertical segments are interquartile ranges from the existing "
       "subtype summary; jitter is for display only. This is the same cohort as A-C. "
       "Evaluated mutation, fusion and CNV "
       "results are shown separately in Supplementary Fig. S3 to preserve their distinct effect units.")

H("ZEB1 covaries with T-lineage genes after subtype adjustment", size=12.5, italic=True)
body(
    "Genome-wide in TARGET (n = 265), ZEB1 correlated positively with the T-cell program "
    "(GATA3, TCF7, BCL11B, IL7R, RAG1) and, weakly, negatively with LYL1 and LMO2; CD34, MYC, "
    "TAL1 and NOTCH1 showed weak associations in this analysis (Figure 4A). "
    "Modeling ZEB1 against the whole transcriptome while adjusting for Liu subtype and age left "
    "the genome-wide association broadly similar (t-statistic Spearman = 0.987 between "
    "unadjusted and adjusted; Figure 4C). However, 10,450/13,090 genes had adjusted "
    "FDR<0.05 and 11,290/13,090 coefficients were positive. This unusually broad same-sign "
    "pattern makes a global expression, QC, or composition axis a plausible competing "
    "explanation; it cannot by itself define a specific ZEB1 program. "
    "A raw STAR count audit of all 265 samples found ZEB1 associated with the unstranded "
    "N_noFeature count (rho = 0.754) and assigned-read fraction (rho = -0.532). The "
    "N_noFeature association persisted after rank adjustment for subtype, age and "
    "sample type in 242 cases (partial rho = 0.770). Adding log1p(N_noFeature) to "
    "the subtype, age and sample-type limma-trend model reduced FDR-significant genes "
    "from 10,442 to 6,839 (10,450 in the original subtype-plus-age model), while "
    "GATA3, TCF7 and BCL11B retained positive associations. LMO2 and ZEB2 estimates "
    "changed materially. These are sensitivity analyses: N_noFeature may carry "
    "biological as well as technical variation (Supplementary Fig. S14). The key exception was "
    "LMO2, whose inverse association attenuated after subtype and age adjustment in TARGET "
    "(FDR = 0.18). A separate St. Jude Pharmacotype cohort (n = 116; 12 ETP) showed "
    "positive ZEB1 associations for GATA3, TCF7 and BCL11B after ETP and age adjustment "
    "(FDR = 0.038, 3.5x10-5 and 0.0033, respectively); LMO2 was inversely associated "
    "(FDR = 4.4x10-6). In an additional independent NOPHO RNA-seq cohort "
    "(GSE272023; 108 tumor patients), BCL11B and IL7R correlated positively with ZEB1 "
    "(rho = 0.267 and 0.313), while LMO2 correlated inversely (rho = -0.307); "
    "these directions persisted after CIMP-class rank adjustment. TCF7 was weakly positive "
    "(rho = 0.137) and GATA3 was near zero (rho = -0.058). Figure 4B compares the same "
    "unadjusted patient-level rank correlations across three independent bulk cohorts; "
    "the distinct adjusted sensitivity models are reported in Source Data and S14. "
    "The fixed GATA3/TCF7/BCL11B mean-z score correlated positively with ZEB1 in "
    "TARGET (rho = 0.479) and Pharmacotype (rho = 0.231). In competitive comparisons "
    "against 10,000 random three-gene sets matched on expression mean and detection "
    "fraction, with or without expression SD, upper-tail empirical P values were "
    "0.26-0.34 and 0.20-0.21, respectively (Supplementary Fig. S14E,F). Thus, the "
    "observed positive associations did not distinguish this selected score from "
    "expression-matched background programs under the tested schemes. "
    "Subtype-adjusted LMO2-ZEB1 associations varied by cohort "
    "(Supplementary Figure 4); adjustment changes the estimand and does not rule out regulation. "
    "Functional signatures showed cohort-dependent directions: several changed "
    "sign between cohorts (e.g. ETP and stemness; Figure 4E), with 0/50 Hallmark "
    "pathways at FDR<0.05; Figure 4F). The subtype-adjusted volcano recovers GATA3, TCF7, BCL11B "
    "and IL7R among the ZEB1-high tail and LMO2/LYL1 among the ZEB1-low tail (Figure 4D). "
    "Gene-level direction did reproduce across cohorts "
    "(rank Spearman = 0.32, P < 1x10-258). Direct residualization of LMO2 versus ZEB1 across "
    "four T-ALL cohorts showed the TARGET inverse correlation (n = 242, rho = -0.23) attenuating "
    "after subtype adjustment (rho = -0.05, P = 0.45; Figure 4G). In GSE110636, however, "
    "the adjusted inverse association persisted (rho = -0.555, P = 0.00395). These observational "
    "results establish heterogeneity, not the absence or direction of a regulatory edge. "
    "The existing TARGET subgroup analysis retained an inverse point estimate in "
    "28 LMO2-associated cases (10 LMO1/2 and 18 LMO2_LYL1; rho = -0.22, P = 0.257), "
    "compared with rho = -0.08 (P = 0.221) among the other 237 cases "
    "(Supplementary Figure 4D,E). These pooled categories contain distinct molecular "
    "subtypes; their correlations do not establish a difference between subtypes "
    "or exclude regulation within a particular cellular background."
)
body(
    "A targeted analysis retaining the available original molecular subtype labels "
    "identified an inverse LMO2-ZEB1 association within TARGET HOXA cases "
    "(n = 33, rho = -0.54, bootstrap 95% CI -0.77 to -0.21; permutation P = 0.0008, "
    "BH q = 0.011 across 14 eligible cohort-subtype tests; Supplementary Table 1). "
    "The LMO2_LYL1 estimate was also negative but imprecise (n = 18, rho = -0.32, "
    "95% CI -0.81 to 0.22). An exploratory HC3 omnibus rank-slope interaction test "
    "indicated heterogeneity across seven eligible TARGET subtypes (n = 234; "
    "P = 0.010, BH q = 0.031 across three evaluable cohorts). These within-cohort "
    "results support subtype-dependent associations; they are not independent "
    "replication of the HOXA finding or evidence of direct repression."
)
figure("Figure4_transcriptome.png")
legend("Figure 4.", "ZEB1-associated T-lineage expression and subtype sensitivity. (A) Selected gene "
       "correlations with ZEB1. (B) Heatmap of unadjusted patient-level Spearman rho for "
       "six selected genes in TARGET (n=265), St. Jude Pharmacotype (n=116) and "
       "GSE272023/NOPHO tumor samples (n=108). Color encodes "
       "direction and strength; cell numbers give rho rounded to two decimals. "
       "GSE272023 provides additional positive BCL11B/IL7R and "
       "inverse LMO2 evidence, while GATA3 is near zero there. P and FDR, with their "
       "different testing families, are in Source Data. "
       "(C) Genome-wide t-statistics remain similar after adjustment. (D) Subtype-adjusted "
       "volcano with selected genes labelled. (E) All eleven shared program correlations "
       "in TARGET and Pharmacotype, with numerical annotations. Panels B and E use "
       "the same -1 to +1 Spearman correlation color scale, without row standardization. "
       "(F) Hallmark enrichment estimates; 0/50 pathways met FDR<0.05. (G) LMO2-ZEB1 Spearman rho before "
       "and after subtype residualization; these are conditional associations, not causal tests.")

H("Malignant-cell restriction changes the apparent single-cell validation", size=12.5, italic=True)
body(
    "In diagnostic single-cell RNA-seq (GSE227122; Figure 5A), ZEB1 was sparsely detected. "
    "Exploratory patient-level associations were imprecise in this discovery cohort "
    "(n = 10; all FDR > 0.8; Figure 5B). The prior GSE248287 analysis pooled malignant and "
    "normal cells before calculating patient means. We reprocessed the 15 diagnosis samples "
    "using the authors' Celltypes_all = Malignant annotation, verified metadata-to-matrix row "
    "alignment from per-cell raw UMI totals, summed malignant-cell UMIs for each patient, and "
    "calculated log2(1 + CPM). All 15 patients had at least 30 malignant cells. None of 12 "
    "analysis-defined ZEB1 correlations passed BH FDR < 0.05 (Figure 5C). Several effect estimates "
    "changed sign relative to the mixed-cell analysis (Figure 5G), so the earlier two-cohort "
    "mini-meta-analysis is invalid and is withdrawn. "
    "A fixed GATA3/TCF7/BCL11B score computed from raw-UMI patient pseudobulk was "
    "positively associated with ZEB1 in the ten GSE227122 candidate-malignant patients "
    "(rho = 0.455, permutation P = 0.193), but not in the 15 author-annotated malignant "
    "GSE248287 patients (rho = -0.289, P = 0.297; Supplementary Fig. S15). The patient "
    "cohorts differ in malignant annotation and precision; neither result is a definitive "
    "test of the context-dependent bulk association. Across minimum-cell thresholds of "
    "30-1,000, the score correlation remained positive in GSE227122 (rho 0.38-0.45; "
    "9-10 patients) and negative in GSE248287 (rho -0.31 to -0.22; 9-15 patients). "
    "ZEB1 detection fractions ranged from 17.8% to 65.1% and 17.2% to 58.2%, "
    "respectively. Descriptive adjustment for median UMI retained both directions; "
    "conditioning on ZEB1 detection fraction attenuated the GSE248287 correlation "
    "to 0.04. Because detection and expression are coupled, this does not identify "
    "a technical cause. Extending Harmony to a 50-iteration cap reached its "
    "convergence criterion at iteration 27. Across two computational runs and "
    "Leiden resolutions 0.6/1.0, the discovery score correlation remained positive "
    "(rho 0.455-0.515; ten patients), while patient malignant-label Jaccard overlap "
    "with the original labels ranged from 0.706 to 0.994 (Supplementary Fig. S6B). "
    "This supports within-discovery association stability under the tested "
    "computational settings, without independently validating the heuristic labels. "
    "A complementary Park/HTA donor-held-out reference classified 68/70 normal "
    "donor-stage profiles correctly (20 donors; balanced accuracy 97.5%) using "
    "18 developmental markers that excluded ZEB1 and the tested core genes. Seven "
    "of ten leukemia patient mean profiles lay beyond the held-out normal "
    "95th-percentile nearest-centroid distance (Supplementary Fig. S2E,F), supporting "
    "caution when translating normal stage similarity into leukemia stage identity. "
    "This distance flag may reflect biology, population or measurement differences. "
    "Mapping the discovery cells onto the "
    "normal thymic reference showed heterogeneous, not "
    "single-stage, assignments (Figure 5D), and higher-ZEB1 patients were descriptively less "
    "CD8SP-like (rho = -0.75, P = 0.013, FDR ns; Figure 5E). In paired diagnosis/end-of-induction "
    "samples, blasts were cleared without a consistent ZEB1 shift (Figure 5F). These data "
    "do not establish a reproducible ZEB1-low functional state or a cross-cohort program."
)
figure("Figure5_singlecell.png")
legend("Figure 5.", "Single-cell analysis. (A) Per-patient lineage composition at diagnosis. "
       "(B) Discovery patient-level ZEB1 association estimates (n=10). (C) Independent validation "
       "restricted to author-annotated malignant cells (GSE248287, n=15; all BH FDR nonsignificant). "
       "(D) Patient-by-stage heatmap of exploratory normal-thymus mapping "
       "probabilities; each column is one diagnosis patient. (E) ZEB1 versus CD8SP-like "
       "mapping (95% data ellipse). (F) Diagnosis-to-EOI blast clearance. (G) Mixed-cell versus "
       "malignant-only patient-level estimates; connected values are different analyses of the "
       "same cohort and are not independent validation.")

H("Chromatin context and perturbation responses vary across developmental models", size=12.5, italic=True)
body(
    "H3K27ac HiChIP (GSE243915) showed higher ZEB1-locus contact frequency in non-ETP T-ALL "
    "than in ETP, CD34+ progenitors or thymus, consistent with the ETP-low association but only "
    "at borderline significance (P = 0.057; Figure 6A) and with no ZEB1-LMO2 trans contacts. "
    "The assayed ZEB1 promoter CpGs were hypomethylated in both CIMP groups (GSE69954; "
    "P = 0.23; Figure 6B); this DNA methylation comparison does not test histone "
    "deacetylation. LMO2 and TAL1 ChIP tracks showed signal in the "
    "ZEB1/ZEB1-AS1 overlapping promoter region in T-ALL lines (GSE154675). "
    "The 12-kb window includes both genes; canonical ZEB1 TSS-centered windows also "
    "overlap ZEB1-AS1. Without input controls these regional signals cannot establish "
    "ZEB1-specific occupancy or repression (Figure 6C; Supplementary Figure 7D). "
    "TLX1 knockdown in ALL-SIL (GSE110635) raised ZEB1 "
    "(log2FC +0.33; Figure 6D), but TLX1 manipulation does not directly test LMO2 repression. "
    "In Lmo2-transgenic mouse contrasts, Zeb1 decreased slightly in young DN2 cells "
    "(log2FC -0.137, FDR 0.423) and increased in young and old DN3a cells and DN thymocytes "
    "(GSE188225 and GSE49164; Figure 6E). These stage-specific contrasts include several "
    "from one study and are not independent replications. In GSE188225, five young mice "
    "contribute both DN2 and DN3a libraries (17 libraries from 12 mouse IDs overall; "
    "Supplementary Information, perturbation sample units). In the reversible Lmo2 model "
    "GSE144035 [31], sample metadata supported six model-matched Dox comparisons: "
    "three evolving and three Lmo2-independent tumors. Lmo2 and Zeb1 both decreased "
    "after Dox in these six pairs; Zeb1 changes ranged from -1.230 to -0.033 on the "
    "log2(1+CPM) scale (Supplementary Figure 12E-F). Two dependent tumors lacked +Dox "
    "libraries. These descriptive model-level changes cannot establish a direct signed "
    "regulatory effect. DN1 Lmo2 knockout (GSE287751; two libraries per arm) "
    "lowered Lmo2 and raised Zeb1 descriptively (log2FC +0.335), a direction compatible "
    "with release from repression but insufficient to establish a direct mechanism. "
    "Resolving the same study's author-supplied single-cell condition labels provided "
    "additional context-specific direction evidence: Lmo2-KO increased Zeb1 at day 3 "
    "in both unprimed and Flt3L-primed cells (KO-minus-control log2(1+CPM) differences "
    "+0.283 and +0.129), while pre-Notch differences were -0.008 and +0.168. "
    "Lmo2 decreased in all four contexts. These eight conditions were analyzed "
    "descriptively and were not treated as independent biological replicates "
    "(Supplementary Fig. S12C). "
    "None of these experiments isolates a direct early effect "
    "at the human ZEB1 promoter. CollecTRI ULM on the TLX1 knockdown "
    "(GSE110635; 3 control, 6 KD) left "
    "no clear ZEB1 or TLX1 TF activity change (P = 0.89 and 0.67); no analysis-defined T-lineage TF "
    "passed FDR < 0.05 (Figure 6F). Re-analysis of the archived network reproduced all 737 TF "
    "effect estimates. Of nine TLX1 targets, seven were detected in at least one sample and six "
    "in every sample; corresponding ZEB1 coverage was 57/63 and 34/63. Separate siRNAs gave "
    "small ZEB1 activity changes of opposite signs (Supplementary Figure 10C). The weak activity "
    "response of the directly knocked-down TF limits interpretation of ZEB1 activity as a negative "
    "mechanistic result. The patient associations and these perturbations do not "
    "establish a context-general signed regulatory edge."
)
figure("Figure6_mechanism.png")
legend("Figure 6.", "Chromatin context and transcriptional responses to perturbation. (A) H3K27ac HiChIP ZEB1-locus contact "
       "rate. (B) Promoter CpGs remain hypomethylated in CIMP+ and CIMP- T-ALL. (C) LMO2/TAL1/"
       "LDB1/GATA2 ChIP signal heatmap over the shared ZEB1/ZEB1-AS1 12-kb region (no input). "
       "Values are exact means over covered bases; coverage fractions and narrower-window "
       "sensitivities are in Source Data. (D) TLX1 knockdown raises "
       "ZEB1. (E) ZEB1/Zeb1 log2FC by perturbation and developmental context. "
       "Rows show GSE110635 (TLX1 KD), GSE188225 (young DN2, young DN3a and old DN3a), "
       "GSE49164 (DN thymocytes) and GSE287751 (DN1 KO), respectively. NA denotes "
       "unavailable FDR, not a nonsignificant test. "
       "The DN1 Lmo2 knockout is shown separately from Lmo2 ON contrasts. "
       "Red and blue denote FDR-significant positive and negative changes; grey denotes "
       "a nonsignificant result or unavailable FDR. Exact values and missing intervals "
       "are in Source Data. "
       "(F) Per-sample CollecTRI ULM activities centered on the control mean for each TF, "
       "without row scaling: three controls and three samples for each siRNA sequence. "
       "Labels preserve the source replicate index. This heatmap displays variability; "
       "complete tests and target coverage are reported in Supplementary Figure 10 and Source Data.")

H("ZEB1 associations with drug response and genetic dependency vary by context", size=12.5, italic=True)
body(
    "Finally, we tested whether ZEB1 state carries a general pharmacological or genetic "
    "vulnerability. "
    "In a 116-patient St. Jude pharmacotype cohort, 18 ex-vivo drugs were recorded; "
    "17 had at least eight paired observations for ZEB1-LC50 analysis (per-drug n=14-94), "
    "and no analyzed drug "
    "reached FDR<0.05 for a ZEB1-LC50 correlation; the strongest nominal signals (panobinostat, "
    "nelarabine) were small-n exploratory tails, and methotrexate was absent from the St. Jude "
    "panel (Figure 7A). Filling that gap in GDSC2, nine T-ALL lines with both ZEB1 RNA and MTX "
    "LN_IC50 showed a positive Spearman correlation (rho = 0.90, exploratory permutation "
    "P = 0.0020). GDSC1 showed the same direction in ten lines (rho = 0.673, permutation "
    "P = 0.0397); eight ModelIDs were shared and their MTX responses correlated across "
    "screens (rho = 0.905; Supplementary Figure 11). Thus the result is cross-screen "
    "consistency among largely the same models, not independent biological replication. "
    "The T-ALL-restricted signal remains power-limited and is not claimed as a biomarker; "
    "the pan-cancer matched GDSC2 set (n = 96, rho = -0.07, P = 0.48) is retained in "
    "Source Data as context, not as a T-ALL validation cohort. Of 92 patients with both "
    "MRD measurements in the St. Jude source table, 26 were below 0.01% at both time points, "
    "48 moved from ≥0.01% at mid induction to below 0.01% later, and 18 remained ≥0.01%; "
    "there were no below-to-above transitions (Figure 7G). "
    "An exploratory low-dimensional model adding diagnosis ZEB1 to continuous mid-induction MRD "
    "reduced leave-one-out AUC from 0.737 to 0.708 and increased Brier score from 0.172 to "
    "0.177; the adjusted ZEB1 OR per SD was 1.16 (95% CI 0.57-2.38; Figure 7G). The "
    "0.01% boundary is confirmed by the source study [2]. Paired MRD records came from "
    "Total XV (23 patients) and XVI (69); all 19 Total XVII patients in our 116-patient "
    "table lacked both MRD measurements. Adding protocol to both comparator models "
    "yielded leave-one-out AUCs of 0.784 without and 0.767 with ZEB1 (Brier 0.164 and "
    "0.169). These internal sensitivities do not establish incremental prediction. "
    "TARGET follow-up on RNA-seq cases "
    "(n = 262) gave event-limited, non-significant Cox models for continuous ZEB1: OS HR per SD "
    "= 1.15 (95% CI 0.71-1.86, 14 events, P = 0.58) and EFS HR = 1.02 (27 events, P = 0.90; "
    "Supplementary Figure 12). ZEB1 showed only a weak, unstable association with mid-induction "
    "minimal residual disease that was not evident at end induction (Figure 7B). "
    "In a single locked DepMap Public 26Q1 release, the eight T-ALL lines with CRISPR data "
    "had mean ZEB1 gene effect -0.54; four crossed the conventional -0.5 reference and four "
    "did not (Figure 7C). Compared with 13 B-ALL lines, the mean T-minus-B effect was -0.55 "
    "(bootstrap 95% CI -0.94 to -0.21; Mann-Whitney P = 0.00596, BH FDR = 0.0318 across "
    "32 selected genes; Figure 7E). This is a model-lineage difference in a small panel, "
    "not uniform dependency, a patient biomarker, or a therapeutic effect. "
    "Selected target gene effects (BCL2L1, CDK6, PSMB5) are shown for context; these models "
    "do not establish a ZEB1-expression-specific vulnerability (Figure 7D)."
)
figure("Figure7_therapy.png")
legend("Figure 7.", "Drug response, genetic dependency and residual disease. (A) ZEB1 versus 17 analyzable ex-vivo drugs (per-drug n=14-94); "
       "CHZ868 had six paired observations and is displayed only in the missingness "
       "summary in Supplementary Fig. S8. "
       "(B) ZEB1 versus minimal residual disease. (C) ZEB1 CRISPR gene effect across T-ALL lines "
       "(heterogeneous; 4/8 below -0.5). "
       "(D) Paired T-ALL and B-ALL mean gene effects for selected targets from the "
       "same DepMap release. (E) T-ALL versus B-ALL gene-effect mean differences "
       "with bootstrap intervals; FDR across 32 genes. "
       "(F) GDSC1 and GDSC2 methotrexate versus ZEB1 in T-ALL (n=10 and n=9; eight "
       "ModelIDs shared, linked by grey lines). (G) Patient-level mid- to end-induction MRD flow (n=92 paired), "
       "with the internally validated AUC change after adding diagnosis ZEB1 among the "
       "66 mid-induction-positive patients. The source-defined 0.01% boundary is used; "
       "mid-induction is day 19 in Total XV and day 15 in Total XVI. Independent "
       "clinical validation was not available.")

# ===========================================================================
# DISCUSSION
# ===========================================================================
H("Discussion", size=14)
body(
    'ZEB1 expression in T-ALL is embedded in developmental and molecular context. Its higher expression relative to B-ALL, variation across normal thymocyte stages and positive TLX-high direction across cohorts provide a framework for interpreting heterogeneous patient and model-system findings. The positive BCL11B and IL7R associations and inverse LMO2 association in NOPHO, together with context-specific perturbation responses, support a biological relationship worth resolving at the level of defined cell states. Bulk associations and selected model systems do not by themselves identify a universal regulatory mechanism; patient single-cell analyses further show that malignant-cell selection can alter the apparent relationship.'
)
body(
    'The clinically relevant question is whether this context helps distinguish treatment response. The methotrexate association is a focused lead: higher ZEB1 tracked higher IC50 in both GDSC screens, although eight overlapping models prevent an independent-replication claim. This direction differs from a general ZEB1-low chemoresistance model and suggests that drug identity, developmental state and model background should be considered together. The patient ex-vivo panel did not contain methotrexate, so its results cannot directly adjudicate this signal. Separately, diagnosis ZEB1 did not improve internal prediction of persistent MRD beyond mid-induction MRD in the available subset. Drug sensitivity in cultured models and incremental patient-level prediction address different questions; the present data connect the former to a specific treatment while leaving clinical translation unresolved.'
)
body(
    "The original mechanistic evidence addresses a distinct experimental comparison. "
    "Wu et al. introduced ZEB1 into LMO2-expressing Jurkat cells and measured EdU-positive "
    "cells after 24 h of MTX exposure (150 nM); they reported restoration of drug "
    "responsiveness after ZEB1 re-expression and attenuation of MTX insensitivity with "
    "TSA (1.5 micromolar; their Figure 5) [1]. They also reported LMO2-associated "
    "promoter histone changes and Sap18-dependent ZEB1 repression (their Figures 2 and 3). "
    "Our GDSC analysis relates baseline expression to IC50 across different cell lines, "
    "without matched LMO2 perturbation or ZEB1 rescue. Its positive association therefore "
    "does not directly retest that within-model rescue experiment. The public-data "
    "findings motivate investigation of which cellular backgrounds preserve the reported "
    "mechanism and how assay endpoints shape the observed drug-response relationship."
)
body(
    "The TLX-high ZEB1 "
    "association is directionally consistent across four cohorts, although heterogeneity "
    "and the small number of studies limit its precision. An exploratory split of the "
    "study-supplied TLX labels in two independent array cohorts found a positive TLX3-versus-"
    "non-TLX contrast in both GSE26713 (g=0.42; 22 versus 88 patients) and GSE62156 "
    "(g=1.11; 12 versus 47), whereas TLX1 had opposite directions (g=-0.49; 7 versus 88, "
    "and g=1.66; 5 versus 47, respectively). These small post-hoc strata do not establish "
    "a TLX3-specific mechanism or independent clinical biomarker (Source Data: "
    "tlx1_tlx3_split.tsv). ETP status and the study-defined "
    "immature molecular classes differ, so their effects cannot be pooled into a single claim. "
    "The value of the study is threefold: it provides a quantitative, multi-omic map of ZEB1 "
    "biology in T-ALL; it quantifies how subtype adjustment changes observed gene associations "
    "across cohorts; and it supplies analysis code and source tables for the displayed results."
)
body(
    "A further methodological safeguard underlies these conclusions. We identified a discordant microarray probeset "
    "(208078_s_at), whose deposited representative sequence reads as ZEB1 but whose "
    "oligonucleotides match SIK1 and an overlapping Ensembl model, but not ZEB1, in the "
    "tested transcriptome. These matches do not establish gene-level SIK1 specificity. "
    "Its signal was "
    "lower in T-ALL than several comparator classes in GSE13159 (Supplementary Fig. S13). By "
    "quantifying ZEB1 exclusively from the four validated ZEB1-specific probesets and from "
    "gene-level RNA-seq, we ensure that the T-lineage, developmental-state-linked behavior reported here "
    "is less vulnerable to that specific annotation error. We have not established whether "
    "the 2018 mechanistic study used this probeset."
)
body(
    "Our study has limitations inherent to public-data reanalysis. Several orthogonal layers are "
    "small (H3K27ac HiChIP n = 7; DepMap T-ALL n = 8; GDSC2 MTX T-ALL n = 9; paired single-cell "
    "n = 5), so we report them as bounded, hypothesis-level evidence rather than definitive "
    "negatives, and we flag each with effect sizes and confidence intervals. The St. Jude "
    "pharmacotype panel lacks methotrexate; GDSC2 fills that gap but the T-ALL n = 9 correlation "
    "is power-limited and is not treated as an MTX biomarker. The pan-cancer GDSC2 set is "
    "context in Source Data, not an independent T-ALL validation cohort. PRISM 24Q2 was "
    "not obtained. Platform and age differences between cohorts contribute to heterogeneity. "
    "The larger childhood T-ALL genomic resource [6] is a candidate for future expression "
    "validation only after access rights, endpoint definitions, and patient overlap with "
    "TARGET are audited. Larger public ETP-enriched malignant-cell series or additional "
    "non-overlapping MTX-screened T-ALL models would strengthen biological validation. "
    "TARGET OS/EFS Cox models on continuous ZEB1 "
    "(n = 262) are event-limited (OS 14 events, P = 0.58; EFS 27 events, P = 0.90); "
    "two subtypes have zero OS events and the subtype-adjusted Cox model shows separation "
    "symptoms. MRD associations are weak and unstable; the 66-patient persistence analysis "
    "did not show incremental value for diagnosis ZEB1 beyond mid-induction MRD. These analyses "
    "do not establish an independent prognostic "
    "association (Supplementary Figure 12). Overall, the current evidence supports bounded "
    "developmental and subtype associations while leaving cell-intrinsic validation unresolved."
)

# ===========================================================================
# METHODS
# ===========================================================================
H("Materials and Methods", size=14)


def mh(t):
    H(t, size=12, before=6, after=2, italic=True)


mh("Datasets")
body(
    "Public datasets analyzed included: GSE13159 (MILE); TARGET-ALL-P2 (STAR TPM); the St. Jude "
    "pharmacotype cohort [2] (RNA-seq FPKM and 18-drug ex-vivo LC50); GSE107011 (healthy immune "
    "profiling); GSE142522 (sorted thymocytes and adult T-ALL); GSE195812 and GSE206710 "
    "(FACS-sorted thymic single-cell RNA-seq); CELLxGENE Human T cells (Park et al. Science "
    "2020); GSE227122 [7] and GSE248287 [4] (T-ALL single-cell/CITE "
    "RNA-seq); GSE62156, GSE26713 and GSE110636 (subtype replication); GSE272023 (CIMP); "
    "GSE110637, GSE243915, GSE69954, GSE154675, GSE70734, GSE280728 and GSE287757 "
    "(chromatin/methylation); GSE110635, GSE188225, GSE186943 and GSE144035 (TLX1/LMO2 "
    "perturbation); GDSC1/GDSC2 methotrexate LN_IC50; TARGET-ALL-P2 clinical follow-up; "
    "and DepMap Public 26Q1 (expression, copy number, CRISPR gene effect; Chronos method [8]). "
    "Accessions and per-sample usage are listed in the Source Data and report accompanying this "
    "manuscript.", first_indent=False)

with open(os.path.join(S.DATA_DIR, "reference_catalog.json"), encoding="utf-8") as ref_file:
    reference_catalog = json.load(ref_file)
body("Dataset-associated source publications: " + "; ".join(
    entry["dataset"] + " [" + ",".join(map(str, entry["reference_numbers"])) + "]"
    for entry in reference_catalog["dataset_citations"]) + ".",
    first_indent=False)

mh("ZEB1 quantification and probe validation")
body(
    "On Affymetrix HG-U133 Plus 2.0 (GPL570) arrays, gene-level ZEB1 was quantified as the mean "
    "of the four ZEB1-specific probesets 210875_s_at, 212758_s_at, 212764_at and 239952_at, "
    "whose 25-mer oligonucleotides map uniquely to ZEB1 on chromosome 10. We explicitly excluded "
    "the frequently reused probeset 208078_s_at: although its GEO representative sequence is a "
    "ZEB1 transcript (NM_030751), all eleven of its physical oligonucleotides match SIK1 "
    "and an overlapping Ensembl transcript model, while 0/11 match ZEB1, in an exact-match "
    "GENCODE v47 transcriptome audit. This does not establish that its array signal "
    "quantitatively measures gene-level SIK1. "
    "RNA-seq ZEB1 used the gene-level quantification of each resource. This checks the "
    "identity of the selected microarray probesets; it does not establish the probes used "
    "in earlier publications.", first_indent=False)
mh("TARGET raw-count QC and expression-wide sensitivity")
body(
    "For 265 TARGET samples, we parsed unstranded assigned and STAR summary rows from "
    "the original STAR ReadsPerGene count files. N_noFeature, assigned-read fraction and "
    "mitochondrial assigned-count fraction were linked to the processed log2(TPM+1) matrix "
    "by sample identifier. Expression PC1 was calculated after excluding ZEB1. Rank-based "
    "partial associations residualized ZEB1 and each QC metric on Liu subtype, age and "
    "bone marrow versus peripheral blood sample type in the 242 complete cases. "
    "Genome-wide sensitivity used the same 242 cases and limma-trend robust fitting, "
    "comparing the subtype-plus-age model against additions of sample type, assigned "
    "fraction, log1p(N_noFeature) or mitochondrial fraction. FDR used Benjamini-Hochberg. "
    "For competitive gene-set controls, each GATA3/TCF7/BCL11B gene was matched to "
    "100 nearest genes by Euclidean distance of percentile-ranked mean log expression "
    "and detection fraction, with a second scheme adding expression SD. ZEB1 and "
    "the three selected genes were excluded from candidate pools. We drew 10,000 "
    "three-gene sets per scheme without repeated genes within a set, standardized "
    "each gene across cohort patients, and correlated the mean score with ZEB1. "
    "The empirical positive-tail probability used (1 + number of random-set rho values "
    "at least as large as observed)/(10,001). Matching did not use ZEB1 correlations. "
    "These are post hoc competitive controls on selected genes, not independent "
    "confirmatory tests; all pools, draws and seeds are supplied. "
    "Raw count metrics are diagnostic covariates, not validated measures of technical "
    "quality or proof of confounding. Scripts, sample-level tables and job provenance are "
    "in total/code and total/data/validation.", first_indent=False)
body(
    "Probe identity was evaluated by exact matching of physical Affymetrix 25-mers, on "
    "both strands, against GENCODE v47 transcript sequences. All transcript/gene hits were "
    "retained, including ambiguous matches. Figure S13B displays each of the 11 physical "
    "208078_s_at oligonucleotides against ZEB1, SIK1 and the overlapping Ensembl model. "
    "The companion audit files, reference checksum, matching code and complete hit tables "
    "are supplied under data/probe_audit. The scope of the manuscript is probe identity "
    "for ZEB1 measurement; attribution to probes used by earlier studies requires their "
    "actual expression inputs.", first_indent=False)

mh("Statistical framework")
body(
    "For the independent GSE272023/NOPHO check, we matched 108 tumor sample IDs to "
    "the deposited VST expression matrix, excluding three non-tumor controls. "
    "Seven Ensembl gene IDs, including ZEB1, were resolved using the locally archived "
    "NCBI/Ensembl symbol map. We computed patient-level Spearman correlations for six "
    "selected genes and BH-adjusted those six P values within this exploratory family. "
    "A rank-residual sensitivity conditioned on CIMP-high versus CIMP-low; CIMP is "
    "not a complete molecular subtype or batch adjustment. Figure 4B displays the "
    "unadjusted rank correlations across three independent bulk cohorts so the visual "
    "effect unit is consistent. The corresponding adjusted TARGET and Pharmacotype "
    "models are separate sensitivity analyses.", first_indent=False)
body(
    "ZEB1 was analyzed as a continuous variable; dichotomizations (median, quartiles) were used "
    "only for display. Tests were two-sided; Benjamini-Hochberg FDR families are specified "
    "in the source tables, with exploratory nominal P values identified separately. Group "
    "differences were quantified with Hedges g and 95% confidence intervals and Wilcoxon tests. "
    "Multi-cohort effects were combined with DerSimonian-Laird random-effects meta-analysis, "
    "reporting Q and I2. For the four-cohort TLX contrast, sensitivity analyses used REML "
    "heterogeneity estimation and modified Hartung-Knapp inference (variance scale bounded "
    "below by one; t distribution with k-1 degrees of freedom), an exploratory prediction "
    "interval with k-2 degrees of freedom, and leave-one-cohort-out fits. Genome-wide "
    "associations used Spearman correlation and, for adjustment, "
    "limma with the model Y ~ ZEB1 + subtype + age; gene-set enrichment used fgsea on the "
    "Spearman-ranked list. Single-cell results were summarized at patient level; GSE248287 "
    "used author-annotated malignant cells and raw-UMI pseudobulk, whereas the GSE227122 "
    "malignant-cell assignment remains heuristic. Paired comparisons used exact "
    "Wilcoxon tests. Cox proportional-hazards models on TARGET OS/EFS used continuous ZEB1 "
    "(per SD), unadjusted and age-adjusted, without an optimal cutoff. An exploratory "
    "subtype-adjusted fit was flagged for zero-event strata and nonfinite/zero confidence "
    "limits and was not used for an independence claim. Endpoint fields were checked against "
    "the TARGET Phase II Validation clinical workbook dated 20230727, preserving the literal "
    "First Event category None as no event and retaining genuinely blank fields as missing. "
    "Independent refits used R survival with Breslow and Efron ties; proportional hazards "
    "were assessed using scaled Schoenfeld residuals with the Kaplan-Meier time transform, "
    "and follow-up was summarized using reverse Kaplan-Meier estimates. CollecTRI "
    "TF activity after TLX1 knockdown used decoupler 2.1.4 ULM (tmin = 5), applied to "
    "natural-log(1 + CPM) expression with the archived OmniPath CollecTRI network. "
    "Inhibitory edges received weight -1, other edges +1, and duplicate source-target "
    "weights were averaged. We exported the frozen network, per-sample activities and "
    "target detection counts, and compared each siRNA and pooled KD against shared controls "
    "using Welch tests with BH adjustment across returned TFs within each contrast. "
    "Evidence was interpreted according "
    "to the cohort, assay and unit of replication for each analysis. Cross-cohort associations "
    "and orthogonal assays address different sources of uncertainty; neither establishes "
    "causality. MRD prediction and nearest-stage mapping remain exploratory and lack "
    "independent clinical or malignant-cell stage validation.", first_indent=False)

mh("Exploratory MRD persistence analysis")
body(
    "In the source pharmacotyping study [2], primary ALL cells were tested at six "
    "concentrations per drug using a four-day MTT assay for twelve agents and a "
    "96-hour mesenchymal stromal co-culture/flow-cytometry assay for the remaining "
    "six. Responses outside the tested range were assigned half the lowest or twice "
    "the highest tested concentration, respectively. LC50 values therefore include "
    "boundary assignments and are not uniformly uncensored measurements. The current "
    "drug-wise rank analyses retain the deposited values and missingness; no drug "
    "imputation is used for the displayed ZEB1 correlations. The panel does not "
    "include methotrexate.", first_indent=False)
body(
    "From the 116-patient St. Jude source table, we first counted transitions among 92 "
    "patients with both mid- and end-induction MRD measurements. We then retained the 66 "
    "patients with mid-induction MRD ≥0.01% for the exploratory persistence model. "
    "Source Methods [2] define negative bone-marrow flow-cytometric MRD as <0.01%. "
    "The source table's day-15 field represents day 19 in Total XV and day 15 in "
    "Total XVI/XVII; end-induction measurements occur at days 42-46. The source assigns "
    "values below detection to 0.005% before logarithmic analysis. Our binary outcome "
    "uses the reported threshold, and the eligible early-positive predictor values need "
    "no below-limit substitution. Reported zeros do not establish biological eradication. "
    "Individual missingness reasons remain unavailable. The binary endpoint was later "
    "MRD ≥0.01%. We compared "
    "a fixed logistic model using log10(mid-induction MRD percent) with the same model plus "
    "diagnosis ZEB1_log. A fixed L2 penalty (C=1) and training-fold standardization were "
    "used for leave-one-out probability predictions. We report ROC AUC and Brier score from "
    "those predictions; a separate unpenalized logistic model supplies a descriptive "
    "ZEB1 odds ratio per SD with Wald confidence interval. This is internal, exploratory "
    "assessment rather than an externally validated prediction model. A separate "
    "sensitivity added a Total XV versus XVI indicator to both comparator models, "
    "retaining the same penalty and fold-specific scaling. All four model results and "
    "patient predictions are reported; no model was selected for significance. Source "
    "MRD values were checked against the 116-patient metadata by ID.", first_indent=False)

mh("Single-cell processing")
body(
    "For GSE227122, the executed log records 31,297 cells after the gene-count and "
    "mitochondrial filters (200 <= nFeature < 8,000; mitochondrial fraction <20%), "
    "followed by exclusion of 121 Scrublet-predicted doublets. Omicverse Harmony was "
    "used, but stopped before convergence at its 10-iteration limit; igraph Leiden "
    "clustering completed. This limits confidence in integration-derived labels and "
    "embeddings. Patient count aggregation uses the unintegrated raw-count layer. "
    "GSE227122 candidate blasts were assigned heuristically: cluster median T-lineage "
    "score had to exceed the myeloid and B scores, with patient entropy below 1.6 bits "
    "or median T score above 0.4; only T/Unknown-lineage cells within qualifying "
    "clusters were retained. An extended-Harmony sensitivity recreated the original "
    "Omicverse default scaled-PCA path, used a 50-iteration cap and random state 0, "
    "and evaluated Leiden resolutions 0.6 and 1.0 with the same label rule. Original "
    "labels were reproduced exactly before this comparison. Patient score sensitivity "
    "used the original three-gene means/SDs, with 10,000 permutations and 2,000 patient "
    "bootstrap resamples; no resolution was chosen by its ZEB1 association. "
    "Both computational runs and all four run/resolution combinations are retained "
    "to report variability in the clustering outputs. The corrected run converged "
    "after 27 Harmony iterations. Execution logs, including the first run's metadata "
    "export error, accompany the code. These are technical sensitivity runs, "
    "not additional biological replicates. "
    "For GSE248287, we used the original authors' malignant-cell "
    "annotation, verified metadata order against raw per-cell UMI totals, required at least "
    "30 malignant cells per diagnosis patient, and summed raw counts before log2(1 + CPM). "
    "The GSE248287 study sorted immune and leukemia fractions and pooled them to "
    "increase immune-cell representation. Observed cell fractions therefore describe "
    "the sequenced mixture and do not estimate the original marrow or blood blast burden. "
    "Patient-level Spearman correlations were BH-adjusted across 12 features defined for this reanalysis. "
    "Per-patient retained cell counts, median UMI depth, detected genes and ZEB1 detection "
    "fractions were audited against the summed-count tables. Minimum malignant-cell "
    "thresholds of 30, 100, 250, 500 and 1,000 were evaluated with the full-cohort "
    "three-gene scaling frozen; all thresholds are reported. These exploratory "
    "sensitivities used 10,000 permutations and 2,000 patient bootstrap resamples "
    "(seed 20260924, deterministically offset by retained patient set). Separate "
    "descriptive partial-rank correlations evaluated one depth/detection covariate at "
    "a time; detection is coupled to expression and these are not causal adjustments. "
    "For the separate fixed three-gene check, we also summed raw counts for "
    "GATA3, TCF7, BCL11B and ZEB1 by diagnosis patient in both cohorts, calculated "
    "log2(1 + CPM) from each retained library, standardized each of the three signature "
    "genes across patients within cohort, and averaged them. GSE227122 uses the prior "
    "project's heuristic candidate-malignant calls; GSE248287 uses author malignant "
    "calls. Two-sided Spearman permutation P values used 100,000 label permutations; "
    "BH correction covered the score and three individual genes within each cohort. "
    "The 10 and 15 patients were analyzed separately without cell-level pooling. "
    "For a complementary reference validation, raw Park/HTA counts and discovery "
    "malignant-cell counts were normalized per cell to log1p(10,000 counts per "
    "retained library total). Eighteen fixed developmental markers excluded ZEB1, "
    "ZEB2, LMO2, GATA3, TCF7, BCL11B and IL7R. Donor-stage mean profiles with at "
    "least 20 cells were retained (70 profiles; six smaller groups excluded). In each "
    "leave-one-donor-out fold, feature scaling and stage centroids used only other "
    "donors; classification used nearest cosine distance. The full reference was "
    "then applied to diagnosis patient mean profiles. Similarity margins, softmax "
    "affinity entropy (temperature 0.2), and an exploratory flag above the held-out "
    "normal 95th-percentile distance were reported. Affinities are uncalibrated, "
    "and no cross-cohort ZEB1 expression residual was interpreted. These profile-level "
    "checks do not validate the separate per-cell mapping. "
    "Discovery cells were mapped to normal thymic stage centroids by cosine distance; "
    "mapping entropy quantified assignment uncertainty.", first_indent=False)
body(
    "For GSE287751, the deposited author cell metadata workbook used Strict OOXML. "
    "We parsed its original worksheet XML and shared-string table, recovered all "
    "11,661 author-retained barcodes, and matched them to the raw UMI matrix. "
    "Author labels defined eight combinations of timepoint, Flt3L priming and "
    "control/Lmo2-KO. Counts were summed within each condition, normalized by total "
    "UMIs to CPM and transformed as log2(1+CPM). We display KO-minus-control "
    "differences within each of four timepoint/priming contexts. These are "
    "descriptive contrasts; no cell-level P values or pseudo-replicate confidence "
    "intervals were calculated.", first_indent=False)

mh("Chromatin, methylation and perturbation analyses")
body(
    "H3K27ac HiChIP contacts (GSE243915) were counted in a ZEB1-centered locus window "
    "(hg38 chr10:31,318,447-31,529,814; 211 kb), corresponding to the annotated ZEB1 gene body "
    "rather than a classical promoter, and normalized per million. We therefore report these "
    "counts as ZEB1-locus contacts. A separate 12 kb promoter-anchored window "
    "(hg38 chr10:31,316,447-31,328,447; 2 kb upstream of the annotated gene start through 10 kb "
    "into the locus) was used only for ChIP occupancy summaries, not for HiChIP counting. "
    "The Ensembl GRCh38 annotation retrieved on 2026-09-24 places ZEB1 on the positive "
    "strand, with canonical transcript ENST00000424869.6 starting at chr10:31,319,216 "
    "(one-based). The legacy 12-kb window overlaps ZEB1-AS1 on the opposite strand. "
    "We re-extracted exact per-base bigWig values in that window, TSS +/-1 kb, TSS +/-2 kb "
    "and the gene body, recording the fraction of covered bases and both covered-base "
    "and missing-as-zero means. The main heatmap uses covered-base means; the locus "
    "display uses 200-bp means with uncovered positions displayed as zero, transformed "
    "by natural-log(1+signal). All 16 GEO sample records identify human samples and "
    "report hg38 (GRCh38) for human read mapping; the deposited processing description "
    "reports HISAT2 mapping and MACS2 peak calling. Track chromosome-10 lengths also "
    "match GRCh38. These records do not specify a track-level normalization formula "
    "or establish matched input controls for the displayed signals. "
    "450K methylation used Illumina TSS200/TSS1500/5'UTR/1stExon promoter probes, which are "
    "independent of these HiChIP windows. Perturbation RNA-seq was analyzed with "
    "DESeq2 (human) or limma-trend (mouse); only knockdowns with a verified on-target effect "
    "were interpreted. ChIP occupancy without input controls was reported descriptively and not "
    "as evidence of repression.", first_indent=False)

mh("Figure generation and reproducibility")
body(
    "All main figures and Supplementary Figures S1-S15 were generated in Python "
    "(matplotlib/seaborn) using the supplied source-data tables, module outputs and "
    "validation exports. The figure scripts retain the panel-specific input paths and "
    "display transformations; panels may combine multiple inputs. S1-S12 are produced by figS_all.py and "
    "the probe-validation panel S13 by figS13_probe_annotation.py and TARGET QC panel "
    "S14 by figS14_target_qc.py, and the fixed patient-level check S15 by "
    "figS15_three_gene_patient.py. "
    "Figures are provided as vector PDFs with embedded, editable text (TrueType) for downstream "
    "editing in Adobe Illustrator. Supplementary PDFs are in figures/Supplementary/.", first_indent=False)

# ===========================================================================
# DECLARATIONS
# ===========================================================================
H("Data and code availability", size=13)
body("The analyzed public matrices and summary tables are identified by accession or resource "
     "in the Methods. Raw-data access requirements and redistribution rights vary by resource "
     "and require final author review. The current working package contains source-data tables "
     "and figure-generation code (total/data and total/code). "
     + AUTHOR_METADATA.get("code_availability", "No public code repository or DOI is currently available."),
     first_indent=False)

H("Author contributions", size=13)
body(AUTHOR["contributions"].split("; (VII)")[0] +
     ". [AUTHOR CONFIRMATION: final approval of the newly revised manuscript.]", first_indent=False)
H("Funding", size=13)
body(AUTHOR["funding"], first_indent=False)
H("Competing interests", size=13)
body(AUTHOR["competing_interests"], first_indent=False)
H("Acknowledgments", size=13)
body(AUTHOR["acknowledgments"], first_indent=False)
H("Ethics statement", size=13)
body(AUTHOR["ethics"], first_indent=False)
H("Use of generative AI", size=13)
body("Codex was used to assist code revision, analysis execution, figure preparation, and "
     "manuscript drafting.", first_indent=False)

H("References", size=13)
ref_path = os.path.join(S.DATA_DIR, "reference_catalog.json")
with open(ref_path, encoding="utf-8") as ref_file:
    core_refs = json.load(ref_file)["records"]
for number, ref in enumerate(core_refs, 1):
    if ref.get("record_type") == "dataset":
        body(f"[{number}] {', '.join(ref['authors_family'])}. {ref['title']} [dataset]. "
             f"NCBI Gene Expression Omnibus; {ref['year']}. {ref['accession']}. {ref['url']}.",
             spacing="onehalf", first_indent=False)
        continue
    first_author = ref["authors_family"][0] if ref["authors_family"] else "[AUTHOR NEEDED]"
    locator = ref["page"] or ref["article_number"] or ""
    citation = (f"[{number}] {first_author} et al. {ref['title']}. "
                f"{ref['journal']}. {ref['year']};{ref['volume'] or ''}:"
                f"{locator}. doi:{ref['doi']}.")
    body(citation, spacing="onehalf", first_indent=False)
body("GSE144035 is cited directly as a GEO dataset; no unverified primary-publication "
     "attribution is assumed. Bibliographic identities were checked against cached "
     "Crossref, PubMed or GEO records; assay-specific supporting passages and integrity "
     "status require source-level review.", first_indent=False)

doc.save(OUT)
print("saved", OUT)
