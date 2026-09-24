# -*- coding: utf-8 -*-
"""
Assemble total/manuscript/ZEB1_TALL_Supplementary.docx

Nature/Cell-like supplementary PDF companion: no continuous line numbers,
1.15 line spacing, 11 pt Times captions, 180 mm figures, inventory table.
"""
import os
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor, Twips, Emu
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import _style as S

PNG = S.SUPP_PNG_DIR
OUT = os.path.join(S.MS_DIR, "ZEB1_TALL_Supplementary.docx")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

FIGS = [
    ("FigureS1_QC_sensitivity.png",
     "Supplementary Figure 1.",
     "Cohort QC and sensitivity analyses for the T-ALL versus comparator ZEB1 elevation. "
     "(A) MILE (GSE13159) sample composition by disease class. "
     "(B) Probe sensitivity: Hedges g from the four-probe ZEB1 mean versus the canonical "
     "212764_at probe alone. Direction and magnitude are preserved. "
     "(C) Leave-one-cohort-out random-effects meta-analysis of T-ALL versus B-ALL. "
     "(D) Platform replication (microarray versus RNA-seq)."),
    ("FigureS2_thymus_atlas.png",
     "Supplementary Figure 2.",
     "Park/HTA Human T-cell atlas (CELLxGENE; donor-level inference) and an independent "
     "three-donor thymus (GSE206710). "
     "(A) Donor-mean ZEB1 in DN, DP, CD4SP and CD8SP thymocytes (donors with ≥20 cells). "
     "ZEB1 is high in DN/DP and lower in SP. "
     "(B) Paired DN → CD4SP ZEB1 within donors (n = 15; Wilcoxon P = 6.1 × 10−5). "
     "(C) GSE206710 stage-mean ZEB1, LMO2 and ZEB2 (DN early → DP mature); ZEB1 does not "
     "rise monotonically (ρ = −0.50, P = 0.67). "
     "(D) Paired CD3− versus CD3+ ZEB1 within GSE206710 donors (Wilcoxon P = 0.25). "
     "(E) Leave-one-donor-out classification of 70 donor-stage mean profiles from 20 "
     "Park/HTA donors using 18 fixed developmental markers, excluding ZEB1, ZEB2, "
     "LMO2, GATA3, TCF7, BCL11B and IL7R. Numbers are profile counts, colors are "
     "row-normalized fractions; 68/70 profiles were correctly classified. "
     "(F) Affinity of ten diagnosis patient mean profiles to the full normal reference. "
     "Asterisks denote nearest-centroid distances above the 95th percentile of the "
     "held-out normal profiles (7/10 patients). Affinity uses softmax cosine similarity "
     "with temperature 0.2 and is not a calibrated cell-state probability. Patient "
     "profiles average potentially heterogeneous malignant cells; this does not "
     "validate per-cell assignment fractions. The unit of reference validation is the donor."),
    ("FigureS3_subtype_CIMP_cis.png",
     "Supplementary Figure 3.",
     "Subtype extras, CIMP and cis-genetic tests. "
     "(A) GSE272023 CIMP-high versus CIMP-low mean ZEB1 expression "
     "(65 versus 43 samples; Hedges g = 0.25, 95% CI -0.13 to 0.64, P = 0.43). "
     "(B) TARGET driver mutations versus ZEB1; none pass FDR < 0.05. "
     "(C) Fusion class versus no-fusion. "
     "(D) No ZEB1 coding mutation/fusion was annotated among 265 RNA-linked records; "
     "sample-level variant callability is unverified, so this is not a lesion prevalence "
     "estimate. Locus copy-number and CNV burden have weak rank associations "
     "with ZEB1 (ρ ≈ 0.03 and 0.02)."),
    ("FigureS4_LMO2_ZEB1_false_edge.png",
     "Supplementary Figure 4.",
     "LMO2–ZEB1 associations before and after subtype adjustment. "
     "(A) Observed LMO2–ZEB1 Spearman correlation heatmap, before and after subtype "
     "residualization across four cohorts. "
     "(B) Spearman ρ before (grey) and after (red) residualizing both genes on molecular subtype "
     "in TARGET (n = 242) and three independent GEO cohorts. "
     "TARGET ρ = −0.23 (P = 3.2 × 10−4) → −0.049 (P = 0.45). GSE110636 (n = 25) remains "
     "significant after adjustment and is reported as cohort heterogeneity. "
     "(C) TARGET limma: selected T-lineage genes retain associations after subtype + age adjustment; "
     "the LMO2 association attenuates "
     "(FDR = 0.18); blue adjusted points indicate FDR >= 0.05. "
     "(D,E) Patient-level TARGET expression from the previously completed subgroup "
     "analysis: LMO2-associated cases (LMO1/2, n = 10; LMO2_LYL1, n = 18) and "
     "the remaining 237 cases. Colors retain the original molecular labels; axes "
     "are log2(TPM + 1), as defined in the original analysis script. Correlations "
     "and nominal P values are reused from that analysis, without model refitting. "
     "These panels include all 265 patients, including 23 with Unknown subtype; "
     "the subtype-residual analysis in A/B includes 242 classified patients. "
     "Neither pooled subgroup is a single molecular subtype, and these displays "
     "are not independent validation or tests of subtype interaction."),
    ("FigureS5_complete_GSEA.png",
     "Supplementary Figure 5.",
     "Gene-set enrichment overview; full tested lists are in the source tables. "
     "(A) Hallmark NES (0/50 pathways FDR < 0.05). "
     "(B,C) Reactome and GO-BP sets selected by absolute NES; color gives -log10 BH FDR. "
     "Histone/ribosome genes occur in the leading edges, so pathway names alone do not "
     "establish a specific methylation or immune mechanism. "
     "(D) Four discovery co-expression modules correlate with ZEB1 (ρ = 0.84–0.87); "
     "these exploratory associations do not establish the presence or absence of a "
     "distinct ZEB1-low state."),
    ("FigureS6_scrna_composition_sensitivity.png",
     "Supplementary Figure 6.",
     "Single-cell QC and cell-composition sensitivity. "
     "(A) Lineage composition of GSE227122 blasts split by ZEB1 bin. "
     "(B) Patient-level Jaccard overlap of the original and extended-Harmony "
     "heuristic malignant labels at Leiden resolutions 0.6 and 1.0 in two computational "
     "runs. Both reached the integration convergence criterion at iteration 27, "
     "but clustering labels differed; all four results are retained. The original "
     "label rule was reproduced exactly before rerunning integration with a "
     "50-iteration cap. This evaluates label stability, not accuracy against "
     "independent malignant-cell ground truth. Original discovery association "
     "tables remain supplied; fixed-score associations are in Figure S15 and "
     "the extended-label sensitivity tables. "
     "(C) Malignant-cell fraction in the 15 GSE248287 diagnosis patients according to the "
     "original authors' Celltypes_all annotation. These are fractions of the sequenced "
     "mixture after sorting and pooling, not estimates of native patient blast burden. "
     "(D) Mixed-cell patient means versus malignant-cell raw-UMI pseudobulk produce materially "
     "different correlations. The earlier two-cohort meta-analysis was withdrawn."),
    ("FigureS7_CIMP_occupancy.png",
     "Supplementary Figure 7.",
     "CIMP, promoter methylation and occupancy extras. "
     "(A) GSE69954 NOPHO 450k: ZEB1 promoter probes remain hypomethylated in CIMP+ and CIMP−. "
     "(B) Gene-body CIMP+ minus CIMP− differences in mean beta values; colored points "
     "indicate BH FDR < 0.05. These group comparisons do not establish an individual-level "
     "methylation-expression relationship or test histone deacetylation. "
     "(C) GSE154675 LMO2/TAL1/LDB1/GATA2 exact mean signal over covered bases in the "
     "legacy 12-kb ZEB1/ZEB1-AS1 shared region (no input control). "
     "(D) Sixteen local ChIP tracks in 200-bp bins around the canonical ZEB1 TSS; "
     "natural-log(1+mean signal), with uncovered bases displayed as zero. Gray shading "
     "in the locus schematic marks the legacy 12-kb window; dashed lines mark TSS +/-1 kb. "
     "Arrows show annotated gene directions and spans clipped to the displayed locus, "
     "not exon structure. The overlapping antisense gene prevents ZEB1-exclusive "
     "attribution even with narrower windows. Raw bins, coverage and +/-1/2-kb "
     "sensitivity summaries are supplied in data/validation/chromatin_windows."),
    ("FigureS8_complete_18drug.png",
     "Supplementary Figure 8.",
     "St. Jude pharmacotype panel: 18 drugs collected, 17 analyzable for "
     "ZEB1 versus LC50 (CHZ868 has only six measurements and is excluded). "
     "(A) Spearman rho for 17 drugs (per-drug n=14-94); 0/17 pass FDR < 0.05 "
     "(dot size reflects n). Methotrexate is absent. "
     "(B) Missingness; HDAC/nelarabine tails have >80% missing LC50. "
     "(C) Drug-by-gene context matrix for ZEB1, ZEB2 and LMO2, reusing the existing "
     "51 gene-drug correlation estimates. Rows are ordered by the source-table drug "
     "class and drug name. Numbers and colors represent Spearman rho on a common "
     "-1 to +1 scale without row standardization; color does not indicate significance. "
     "The source table retains n, P and FDR for each estimate. "
     "(D) Univariate OLS P for the ZEB1 term; FDR remains non-significant."),
    ("FigureS9_complete_DepMap.png",
     "Supplementary Figure 9.",
     "Complete DepMap 26Q1 T-ALL view. "
     "(A) ZEB1 CRISPR gene effect in every T-ALL line with CRISPR data (n = 8), "
     "with a row-aligned annotation of the release-provided dependency probability "
     "(P(dep.), color scale 0-1; displayed values rounded to two decimals). "
     "Probability is a separate DepMap measure, not a P value or FDR. "
     "(B) ZEB1 RNA versus self-dependency. "
     "(C) Heatmap of the full selected druggable/lineage panel across T-ALL (8), "
     "B-ALL (13) and myeloid (45) cell lines, ordered by T-ALL mean gene effect. "
     "All means use the same 26Q1 release; lineage differences do not establish "
     "selectivity in patients."),
    ("FigureS10_perturbation_TF.png",
     "Supplementary Figure 10.",
     "Perturbation extras and CollecTRI ULM TF activity (decoupler 2.1.4; GSE110635). "
     "(A) TLX1 knockdown (on-target log2FC = −1.17) raises ZEB1 (log2FC = +0.33). "
     "(B) Both siRNA sequences show ZEB1 up; they share the same controls. "
     "(C) Heatmap of mean inferred CollecTRI ULM activity changes for each siRNA "
     "(3 samples each) and their pooled comparison (6 KD) against 3 shared controls. "
     "Row labels give targets detected in at least one sample / all network targets. "
     "These are raw activity changes, without row standardization. Pooled ZEB1 activity Δ = −0.007 "
     "(P = 0.89) and TLX1 activity P = 0.67; no analysis-defined T-lineage TF passes FDR < 0.05 "
     "(MEF2C nominal P = 6.7 × 10−4, FDR = 0.14). TLX1 target coverage was 7/9 and "
     "ZEB1 coverage 57/63; only 6 and 34 targets, respectively, were detected in all samples. "
     "The archived 62,411-edge network reproduced all 737 returned TF effects; "
     "separate-siRNA tests use BH across 737 TFs within each contrast and share controls. "
     "Coverage counts do not establish assay sensitivity; the weak on-target TLX1 activity "
     "response precludes interpreting ZEB1 non-significance as a strong negative mechanism test. "
     "(D) Mouse Lmo2 transgenic RNA-seq (GSE188225): Zeb1 effects and 95% confidence "
     "intervals for young DN2, young DN3a and old DN3a contrasts. All contrasts derive "
     "from one study and are not independent cohort replications."),
    ("FigureS11_MTX_PRISM.png",
     "Supplementary Figure 11.",
     "Direct methotrexate sensitivity across two GDSC screens. "
     "(A) GDSC1: ten T-ALL models, ZEB1 versus MTX LN_IC50 Spearman ρ = 0.673; "
     "200,000-permutation P = 0.0397; bootstrap 95% CI 0.013–0.925. "
     "(B) GDSC2: nine T-ALL models, ρ = 0.900; permutation P = 0.0020; "
     "bootstrap 95% CI 0.402–1.000. These P values are exploratory. "
     "(C) Eight shared ModelIDs have cross-screen response ρ = 0.905. "
     "(D) Leave-one-cell-line-out ρ values; full ModelID overlap and omissions are in "
     "total/data/validation/. The shared lines make this cross-screen consistency, "
     "not independent biological validation or a patient treatment biomarker."),
    ("FigureS12_PDX_survival.png",
     "Supplementary Figure 12.",
     "PDX L-IC / mouse perturbation and TARGET OS/EFS Cox. "
     "(A) GSE260697 REL13 and REL14 PDXs: individual FPKM libraries for ZEB1, LMO2 and "
     "CD1A, plotted as log2(1+FPKM). Colors identify PDX model and circles/triangles "
     "identify L-IC/DP compartments; REL13 has three libraries per compartment and "
     "REL14 has two. "
     "Points are libraries, not independent patients; no sample pairing is inferred. "
     "The standardized cycle signature is excluded from this shared expression axis. "
     "(B) GSE287751 DN1 Lmo2 knockout lowers Lmo2 and raises Zeb1 descriptively "
     "(log2FC +0.335; two libraries per arm), compatible in direction with release from repression. "
     "(C) GSE287751 author-labelled raw-UMI condition pseudobulk (11,661 retained cells "
     "in eight experimental conditions). U denotes unprimed, F Flt3L-primed, preN "
     "pre-Notch, and d3 day 3 of OP9-DLL1 culture. Cells were assigned using the "
     "author's Strict OOXML metadata, which explicitly identifies control and Lmo2-KO. "
     "Zeb1 KO-minus-control log2(1+CPM) differences were +0.283 (d3 U), +0.129 (d3 F), "
     "-0.008 (preN U) and +0.168 (preN F); Lmo2 decreased in all four contexts. "
     "These are descriptive condition contrasts, without cell-level significance tests "
     "or a claim that the four contexts are independent biological replicates. "
     "(D) TARGET Cox on continuous ZEB1 (n = 262 RNA-seq cases): OS 14 events, HR per SD = 1.15 "
     "(95% CI 0.71–1.86, P = 0.58); EFS 27 events, HR = 1.02 (P = 0.90). Unadjusted and "
     "age-adjusted models are shown. The subtype-adjusted model has zero-event strata and "
     "unstable intervals; see Source Data. No optimal-cutoff survival. ZEB1 is not claimed as "
     "an independent prognostic factor. "
     "(E-F) GSE144035 model-matched Zeb1 and Lmo2 expression before/after Dox: "
     "three evolving (red) and three independent (blue) tumor models, one library per "
     "model and condition. Expression is log2(1+CPM), normalized over all supplied count rows. "
     "Lines join the same VTL model, without claiming the same recipient mouse or independent "
     "biological replication. Both genes decline in all six paired models. Two dependent "
     "models have only -Dox libraries and are not included in the paired display. "
     "The author's Independent -Dox coefficient is not an OFF/ON perturbation effect and "
     "has been removed from Figure 6E; source metadata and the withdrawn row remain archived."),
    ("FigureS13_probe_annotation.png",
     "Supplementary Figure 13.",
     "Probe-annotation validation of ZEB1 quantification (companion GPL570 oligonucleotide "
     "audit). "
     "(A) Four-layer identity of the candidate ZEB1 probesets: 208078_s_at carries a ZEB1 GEO "
     "representative sequence (NM_030751) but all 11 oligonucleotides match SIK1 and an "
     "overlapping Ensembl transcript model, with 0/11 matching ZEB1. This does not establish "
     "gene-level SIK1 specificity. The four probesets used throughout this study "
     "(210875_s_at, 212758_s_at, "
     "212764_at, 239952_at) map uniquely to ZEB1 on chromosome 10. "
     "(B) Exact-match matrix for each of the 11 physical 208078_s_at 25-mers against "
     "GENCODE v47, retaining both SIK1 and overlapping ENSG00000275993 hits. "
     "(C) Individual-sample signal distributions in GSE13159 classes with at least 30 "
     "samples; boxes show the median and interquartile range, whiskers extend to the "
     "most extreme observations within 1.5 interquartile ranges, and outliers are omitted "
     "from display. This is the signal of the non-ZEB1 probeset. "
     "(D) Complete five-probeset Spearman correlation matrix in 117 GSE26713 samples, "
     "including correlations among the validated ZEB1 probes. Probe agreement alone is "
     "not a sequence-identity test."),
    ("FigureS14_TARGET_global_axis_QC.png",
     "Supplementary Figure 14.",
     "TARGET expression-wide association and raw STAR count sensitivity. "
     "(A) ZEB1 associations with expression PC1 and raw-count summary metrics, including "
     "subtype, age and sample-type adjusted rank correlations in 242 cases. "
     "(B) ZEB1 versus unstranded STAR N_noFeature count in 265 samples, coloured by subtype. "
     "(C) Selected gene coefficients from limma-trend models with and without log1p(N_noFeature); "
     "error bars are approximate plus/minus 1.96 standard errors. "
     "(D) Genome-wide FDR-significant gene counts under alternative sample-level covariates. "
     "(E,F) Empirical cumulative distributions for 10,000 random three-gene sets per "
     "matching scheme in TARGET and Pharmacotype, respectively. Genes were matched "
     "on mean log expression and detection fraction, or these variables plus expression "
     "SD. The dashed line is the observed GATA3/TCF7/BCL11B score correlation. "
     "Competitive P denotes the unadjusted empirical positive-tail probability, "
     "not a patient-permutation P value. The two matching schemes are sensitivity "
     "analyses of the same selected score, not independent replications. "
     "N_noFeature can reflect technical and biological variation; these diagnostics do not "
     "identify a causal source of the global axis."),
    ("FigureS15_three_gene_patient.png",
     "Supplementary Figure 15.",
     "Patient-level malignant-cell check of a fixed GATA3/TCF7/BCL11B score. "
     "(A) GSE227122: raw-UMI patient pseudobulk from diagnosis cells assigned candidate "
     "malignant by the project heuristic (n=10). Spearman rho=0.455; permutation P=0.193. "
     "(B) GSE248287: diagnosis malignant cells selected by the original authors' "
     "Celltypes_all annotation (n=15). Spearman rho=-0.289; permutation P=0.297. "
     "The score is the mean of within-cohort patient z scores for the same three genes. "
     "(C) Individual gene and score correlations. No association passes within-cohort "
     "BH FDR<0.05. (D) Per-patient median raw UMI depth and fraction of retained "
     "malignant cells with at least one ZEB1 UMI. (E) Fixed-score correlations across "
     "five minimum-cell thresholds; n denotes retained patients. All threshold-specific "
     "permutation tests, bootstrap intervals and exclusions are supplied in "
     "data/validation/scrna_patient_qc/cell_threshold_sensitivity.tsv. No threshold was "
     "selected for significance. The cohorts have different malignant-call provenance and are not "
     "pooled. These exploratory results describe context and precision; they do not "
     "negate the independently replicated bulk associations."),
]

doc = Document()
normal = doc.styles["Normal"]
normal.font.name = "Times New Roman"
normal.font.size = Pt(11)
normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
pf = normal.paragraph_format
pf.line_spacing = 1.15
pf.space_after = Pt(0)

sec = doc.sections[0]
sec.page_width = Inches(8.5)
sec.page_height = Inches(11)
sec.left_margin = sec.right_margin = Inches(0.85)
sec.top_margin = sec.bottom_margin = Inches(0.85)
# explicitly NO line numbers

footer = sec.footer
fp = footer.paragraphs[0]
fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = fp.add_run("Supplementary Information  ·  ")
r.font.size = Pt(9)
r.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
fld = OxmlElement("w:fldSimple")
fld.set(qn("w:instr"), "PAGE")
fp._p.append(fld)


def shade_header(cell, hex_color="1A1A1A"):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), hex_color)
    shd.set(qn("w:val"), "clear")
    tcPr.append(shd)


def set_run_font(run, size=11, bold=False, italic=False, color=None, name="Times New Roman"):
    run.font.name = name
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    run._element.rPr.rFonts.set(qn("w:eastAsia"), name)
    if color:
        run.font.color.rgb = RGBColor(*color)


def H(text, size=14, before=0, after=8):
    p = doc.add_paragraph()
    p.paragraph_format.keep_with_next = True
    p.paragraph_format.space_before = Pt(before)
    p.paragraph_format.space_after = Pt(after)
    p.paragraph_format.line_spacing = 1.15
    r = p.add_run(text)
    set_run_font(r, size=size, bold=True)
    return p


def body(text, size=11, justify=True, first=False, space_after=8):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.15
    p.paragraph_format.space_after = Pt(space_after)
    if justify:
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    if first:
        p.paragraph_format.first_line_indent = Inches(0.25)
    r = p.add_run(text)
    set_run_font(r, size=size)
    return p


def caption(title, text):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.0
    p.paragraph_format.keep_together = True
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(10)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    r = p.add_run(title + " ")
    set_run_font(r, size=10, bold=True)
    r2 = p.add_run(text)
    set_run_font(r2, size=10)


# ----- cover ----------------------------------------------------------------
t = doc.add_paragraph()
t.alignment = WD_ALIGN_PARAGRAPH.CENTER
t.paragraph_format.space_after = Pt(4)
r = t.add_run("SUPPLEMENTARY INFORMATION")
set_run_font(r, size=11, bold=True, color=(0xC1, 0x27, 0x2D))

t = doc.add_paragraph()
t.alignment = WD_ALIGN_PARAGRAPH.CENTER
t.paragraph_format.space_after = Pt(10)
r = t.add_run("Context-dependent ZEB1 expression and treatment-response associations "
              "in T-cell acute lymphoblastic leukemia")
set_run_font(r, size=14, bold=True)

t = doc.add_paragraph()
t.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = t.add_run("Supplementary Figures S1-S15  |  15 vector PDFs in figures/Supplementary/")
set_run_font(r, size=11, italic=True, color=(0x55, 0x55, 0x55))

body(
    "This file accompanies the main manuscript. Figures follow the same colour-blind-safe "
    "palette, Arial labelling and vector-PDF export (TrueType-42) as the seven main figures. "
    "Source Data tables live in total/data/ (SourceData_FigS*). The post-hoc TLX1/TLX3 "
    "split uses two independent array cohorts and is reported in "
    "total/data/validation/tlx1_tlx3_split.tsv with its input and executable script. "
    "Additional analyses are mapped onto S2 (Park/HTA + GSE206710), S4 (LMO2–ZEB1 heterogeneity), "
    "S6 (scRNA composition sensitivity), S10 (CollecTRI ULM), S11 (GDSC MTX) and S12 (PDX + TARGET Cox). "
    "S13 documents the GPL570 probe-annotation audit; S14 reports TARGET raw-count QC "
    "sensitivity; S15 checks a fixed three-gene score in patient malignant-cell pseudobulk.",
    first=False)

H("Inventory", size=13, before=8, after=6)

table = doc.add_table(rows=1, cols=3)
table.style = "Table Grid"
table.alignment = WD_TABLE_ALIGNMENT.CENTER
hdr = ["Figure", "Maps to", "Evidence summary"]
for i, h in enumerate(hdr):
    cell = table.rows[0].cells[i]
    cell.text = ""
    p = cell.paragraphs[0]
    r = p.add_run(h)
    set_run_font(r, size=9, bold=True, color=(255, 255, 255))
    shade_header(cell, "1A1A1A")

inventory = [
    ("S1", "M1 QC / sensitivity", "T-ALL ZEB1-high is probe- and platform-robust"),
    ("S2", "§3.2 Park/HTA + GSE206710", "Donor-level ZEB1 high in DN/DP, lower in SP"),
    ("S3", "M3 CIMP / cis-genetic", "Observed CIMP/copy-number associations and queried lesion annotations"),
    ("S4", "§3.1 / M3–M4", "LMO2–ZEB1 conditional associations differ by cohort"),
    ("S5", "M4 complete GSEA", "Pathway enrichment and sensitivity to the broad expression axis"),
    ("S6", "§3.5 / M5", "Malignant-cell sensitivity; prior mixed-cell validation withdrawn"),
    ("S7", "M6 CIMP / occupancy", "Promoter and gene-body methylation by CIMP; occupancy without input"),
    ("S8", "M7 drug panel", "18 collected, 17 analyzable; 0/17 FDR<0.05; MTX absent"),
    ("S9", "M8 complete DepMap", "ZEB1 effects are heterogeneous across eight T-ALL lines; T–B mean differs"),
    ("S10", "§3.4 / M9 + decoupler", "TLX1 KD expression response, inferred activity and age/stage-specific Lmo2 contrasts"),
    ("S11", "3.3 GDSC MTX", "GDSC1/2 show concordant MTX associations in overlapping T-ALL models"),
    ("S12", "§3.6 / M10 + TARGET", "PDX libraries, context-specific mouse perturbations and event-limited survival"),
    ("S13", "GPL570 probe audit", "Study uses ZEB1-specific probes; 208078_s_at does not map to ZEB1"),
    ("S14", "TARGET raw STAR QC", "Global axis is sample-metric sensitive; selected T-lineage genes persist"),
    ("S15", "Fixed three-gene scRNA check", "Patient-level candidate/author malignant associations differ by cohort"),
]
for fig, maps, claim in inventory:
    row = table.add_row().cells
    for i, val in enumerate((fig, maps, claim)):
        row[i].text = ""
        p = row[i].paragraphs[0]
        r = p.add_run(val)
        set_run_font(r, size=9, bold=(i == 0))

for row in table.rows:
    for cell in row.cells:
        for p in cell.paragraphs:
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.space_before = Pt(2)

# Keep full S12 methods/provenance in the text, with a concise on-page caption.
for index, entry in enumerate(FIGS):
    if entry[0].startswith("FigureS12_"):
        H("Additional details for Supplementary Figure 12", size=13, before=12, after=6)
        body(entry[2])
        FIGS[index] = (entry[0], entry[1],
            "PDX expression, mouse perturbations and TARGET survival. "
            "(A) ZEB1, LMO2 and CD1A in REL13/REL14 L-IC and DP libraries "
            "(GSE260697; log2[1+FPKM]); points are libraries, not independent patients. "
            "(B) GSE287751 DN1 knockout contrasts, two libraries per arm. "
            "(C) Author-labelled GSE287751 condition pseudobulk contrasts "
            "(log2[1+CPM]); U, unprimed; F, Flt3L-primed; preN, pre-Notch; d3, "
            "day 3. The four contexts are descriptive, not four biological replicates. "
            "(D) Unadjusted and age-adjusted ZEB1 Cox estimates, HR per SD with "
            "95% intervals (262 patients; OS 14 events, EFS 27). "
            "(E,F) GSE144035 model-matched Zeb1/Lmo2 before and after Dox, "
            "log2(1+CPM): three evolving and three independent models, one library "
            "per condition; lines indicate model identity, not recipient-mouse pairing. "
            "Full sample definitions, numerical contrasts and source qualifications "
            "appear in Additional details for Supplementary Figure 12.")
        break

H("Cohort reuse, overlap and data access", size=13, before=12, after=6)
body(
    "Exact identifier intersections are recorded in data/validation/cohort_overlap/"
    "overlap_edges.tsv, with per-analysis memberships and hashed source paths. "
    "The 242 molecularly classified TARGET patients are a subset of the 265-patient "
    "cohort; the 234 patients in the fine-subtype interaction analysis exclude the "
    "23 Unknown and eight TAL2 cases. All 28 LMO-associated cases are included in "
    "these nested analyses. These panels reuse one cohort rather than providing "
    "additional independent patient validation. Of 116 Pharmacotype patients, 92 "
    "have paired MRD and 66 meet the persistence-model entry criteria; both are "
    "nested subsets of the same study.")
body(
    "The MTX analyses include 10 GDSC1 and nine GDSC2 models, with eight shared "
    "ModelIDs (two GDSC1-only and one GDSC2-only). The eight-model DepMap ZEB1 CRISPR "
    "set shares five ModelIDs with GDSC1 MTX and six with GDSC2 MTX. These assays "
    "measure different responses but partly reuse the same biological models. "
    "Five young GSE188225 mouse IDs also occur in both DN2 and DN3a. Exact matching "
    "was restricted to compatible ID namespaces; different accessions, papers or "
    "unmatched identifiers do not establish that patients or donors are independent.")
body(
    "The GEO accession registry (data/dataset_registry.tsv) records deposited "
    "designs, supplementary-file links, parent-series relations and available "
    "raw-access statements from cached source records. A listed processed-file "
    "URL or SRA relation is not verification of access to every raw file. "
    "Controlled raw-data access is not inferred from availability of processed "
    "matrices. The present figure package uses the archived exports and local "
    "validation tables documented alongside each analysis; the exact-ID overlap "
    "inventory is limited to the enumerated sets and is not a universal "
    "cross-study patient linkage.")

H("Perturbation sample units and shared specimens", size=13, before=12, after=6)
body(
    "The accompanying perturbation_unit_inventory.tsv in data/validation/perturbation_units "
    "maps the displayed perturbations to their sample manifests and independence limits. "
    "For GSE188225, all 17 library-to-GSM mappings, mouse IDs, ages, genotype and sex were "
    "checked against the deposited GEO SOFT records. These libraries represent 12 mouse IDs. "
    "Five young mice (two WT and three TG) each contribute both DN2 and DN3a libraries; "
    "young-stage contrasts therefore share biological specimens. Young DN2 has two WT "
    "and three TG mice, young DN3a three per genotype, and old DN3a three per genotype. "
    "The stage contrasts are displayed separately and are not pooled as independent studies. "
    "For GSE49164, deposited sample titles identify six WT and three Lmo2-TG arrays, plus "
    "three Lmo2-TG/Lyl1-knockout arrays in a separate arm; mouse identities are not resolved "
    "by this metadata check. Figure 6 uses the TG-versus-WT comparison. GSE110635 siRNA "
    "comparisons share three controls. GSE144035 pairs are model-matched, with recipient-animal "
    "independence unresolved. GSE287751 single-cell conditions are descriptive experimental "
    "contexts, not cell-level biological replicates. This inventory adds provenance without "
    "refitting the existing expression models.", first=False)

H("Notes on analyses 3.1–3.6", size=13, before=12, after=6)
body(
    "3.1 Conditional LMO2–ZEB1 association. In TARGET (n = 242 with assigned subtype) the LMO2–ZEB1 Spearman "
    "correlation is ρ = −0.23 (P = 3.2 × 10−4) and collapses to ρ = −0.049 (P = 0.45) after "
    "residualizing both genes on Liu subtype. GSE62156 behaves the same way. GSE26713 never "
    "shows a negative association. GSE110636 (n = 25) remains inversely associated after "
    "adjustment (ρ = −0.555, P = 0.00395). This heterogeneity does not resolve regulation "
    "or causality.", first=False)
body(
    "3.2 Second normal atlas. Park et al. Science 2020 Human T cells were obtained from "
    "CELLxGENE (donor-level means; n_cells ≥ 20). ZEB1 is high in DN (mean 1.45, 20 donors) "
    "and DP (1.37) and lower in CD4SP/CD8SP (0.61/0.65). GSE206710 (three donors) remains "
    "the independent second atlas: ZEB1 is detectable in CD3+ cells and is not a monotonic "
    "stage program.", first=False)
body(
    "3.3 MTX / GDSC. St. Jude pharmacotype lacks methotrexate. In GDSC1, ten T-ALL "
    "models with ZEB1 RNA and MTX LN_IC50 gave Spearman rho = 0.673 (200,000-permutation "
    "P = 0.0397); GDSC2 gave rho = 0.900 in nine models (P = 0.0020). Eight ModelIDs "
    "overlap, so the screens show assay consistency rather than independent biological "
    "validation. Cross-screen MTX response rho in the shared eight models was 0.905. "
    "The association is exploratory and does not establish a patient biomarker. "
    "PRISM 24Q2 was not obtained (figshare 403); pharmacotype drug proxies and DepMap "
    "DHFR gene effect do not test the same MTX exposure.", first=False)
body(
    "3.5 scRNA reanalysis. The earlier GSE248287 correlations used all cell types and the "
    "resulting Fisher-z mini-meta is withdrawn. We used the authors' malignant-cell annotation "
    "for all 15 diagnosis patients, verified metadata order against per-cell RNA totals, summed "
    "raw UMIs by patient and computed log2(1 + CPM). All patients passed the minimum of 30 "
    "malignant cells; none of 12 analysis-defined correlations passed BH FDR < 0.05. "
    "GSE227122 blast selection is heuristic and remains exploratory.", first=False)
body(
    "3.6 Survival. TARGET-ALL-P2 clinical follow-up was merged with RNA-seq cases (n = 262). "
    "Cox models used continuous ZEB1 (per SD), not an optimal cutoff. OS: 14 events, HR = 1.15 "
    "(P = 0.58); EFS: 27 events, HR = 1.02 (P = 0.90). Subtype-adjusted models have "
    "zero-event strata (NKX2_1 and Unknown for OS), with three zero/infinite confidence "
    "limits across subtype coefficients. They are not used to infer independence. "
    "The 265 linked patients were checked one-to-one against the original TARGET Phase II "
    "Validation clinical workbook (20230727); 262 had positive follow-up and complete model "
    "fields. First Event contained 230 literal None entries (no event) and three genuinely "
    "blank entries; blanks were not recoded as censored. Independent R survival refits "
    "reproduced the original Breslow hazard ratios within 5e-9; Efron estimates were similar. "
    "Reverse Kaplan-Meier median follow-up was 2289 days for OS (95% CI 2191-2325) and "
    "2302 days for EFS (2203-2334). Under Breslow ties, scaled Schoenfeld-residual tests "
    "for ZEB1 gave P = 0.738 (OS) and 0.286 (EFS) in unadjusted models, and 0.733 and "
    "0.294 after age adjustment. The age-adjusted EFS global test gave P = 0.099; "
    "limited event counts constrain detection of departures from proportional hazards. "
    "All diagnostics and patient inclusion flags are in data/validation/target_survival. "
    "Mid-/end-induction MRD analyses do not establish ZEB1 as an independent "
    "prognostic factor. Figure 7G uses 92 paired patients for its MRD flow, with "
    "the source-confirmed 0.01% boundary and counts in "
    "total/data/validation/mrd_transitions.tsv. In the exploratory paired-MRD subset, "
    "66 patients were mid-induction positive and 18 remained positive at end induction. "
    "Source Methods (Lee et al., PMID 36604538) define mid induction as day 19 for "
    "Total XV and day 15 for XVI/XVII, with end induction at days 42-46. The paired "
    "subset contains 23 Total XV and 69 Total XVI patients; all 19 Total XVII patients "
    "in the analysis table lack both MRD fields. Adding diagnosis ZEB1 reduced "
    "leave-one-out AUC from 0.737 to 0.708. With protocol included in both models, "
    "AUC was 0.784 without and 0.767 with ZEB1. Patient manifests, model predictions "
    "and inspected source-method passages are in data/validation/clinical_endpoint_audit. "
    "These models remain internally assessed exploratory analyses.",
    first=False)

doc.add_page_break()

# ----- targeted within-subtype supplement ----------------------------------
H("Within-subtype LMO2–ZEB1 associations", size=12)
body("Original molecular labels were retained separately within each cohort; "
     "TLX1/TLX3 and LMO1/2/LMO2_LYL1 were not merged where separately available. "
     "Strata with at least 10 complete patients and variable expression underwent "
     "two-sided Spearman permutation tests (9,999 permutations, plus-one correction) "
     "and paired patient bootstrap percentile intervals (2,000 resamples). Unknown "
     "labels and smaller strata remained descriptive. BH adjustment covered all 14 "
     "eligible cohort–subtype tests. The sample threshold and tests were defined for "
     "this targeted reanalysis; results are exploratory. All strata, sample IDs and "
     "input hashes are in data/validation/lmo2_within_subtype. Seed: 20260924, "
     "offset by sorted cohort–subtype index.")
body("To assess heterogeneity within each evaluable cohort, expression ranks were "
     "centered and scaled to unit population standard deviation within subtype. "
     "OLS models included subtype intercepts, LMO2 rank and subtype-by-LMO2 "
     "interactions. An HC3 Wald chi-square test jointly assessed the interaction "
     "terms, with BH adjustment across three evaluable cohorts. This tests differences "
     "in standardized rank slopes and uses an asymptotic approximation; small strata "
     "limit precision. GSE110636 had only one eligible subtype and was not tested. "
     "Expression-derived subtype definitions can affect conditional associations; "
     "these observational tests do not establish a regulatory mechanism.")
caption("Supplementary Table 1.", "Within-subtype estimates. Intervals are unadjusted "
        "95% percentile bootstrap intervals; q values adjust permutation P values "
        "across the 14 eligible strata. Sparse strata are retained in the source table.")
within = pd.read_csv(os.path.join(S.DATA_DIR, "validation/lmo2_within_subtype/within_subtype.tsv"), sep="\t")
table = doc.add_table(rows=1, cols=5)
table.style = "Table Grid"
for cell, title in zip(table.rows[0].cells, ["Cohort", "Original subtype", "n", "rho [95% CI]", "BH q"]):
    cell.text = title
for r in within[within.inference].itertuples():
    values = [r.cohort, r.subtype, str(r.n),
              f"{r.rho:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}]", f"{r.q_bh:.3f}"]
    for cell, value in zip(table.add_row().cells, values):
        cell.text = value
for row in table.rows:
    for cell in row.cells:
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                set_run_font(run, size=8)
interaction = pd.read_csv(os.path.join(S.DATA_DIR, "validation/lmo2_within_subtype/cohort_interactions.tsv"), sep="\t")
body("Omnibus rank-slope interaction results: " + "; ".join(
    f"{r.cohort}, n={r.n}, {r.subtypes} subtypes, P={r.p_hc3:.4f}, BH q={r.q_bh:.4f}"
    for r in interaction[interaction.p_hc3.notna()].itertuples()) + ".", size=10)
doc.add_page_break()

# ----- figures --------------------------------------------------------------
for i, (png, title, legend) in enumerate(FIGS):
    path = os.path.join(PNG, png)
    if not os.path.exists(path):
        body(f"[Missing {png}]", size=10)
        continue
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.keep_with_next = True
    p.paragraph_format.line_spacing = 1.0
    picture_width = 6.2 if png.startswith(("FigureS12_", "FigureS14_")) else 6.65
    p.add_run().add_picture(path, width=Inches(picture_width))
    caption(title, legend)
    if i < len(FIGS) - 1:
        doc.add_page_break()

doc.save(OUT)
print("saved", OUT)
