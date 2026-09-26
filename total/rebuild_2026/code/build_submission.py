"""Assemble the journal submission package in ZEB1/submission.

Inputs: polished manuscript Markdown, Supplementary Information Markdown,
finalized figures (finalize_journal_figures.R) and supplementary tables.
Outputs: manuscript DOCX (legends only and with embedded figures),
Supplementary Information DOCX, Supplementary Tables XLSX, cover letter DOCX.
"""
from pathlib import Path
import csv
import re
import shutil
import subprocess

import openpyxl
from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt
from openpyxl.styles import Alignment, Font

CODE = Path(__file__).resolve().parent
REBUILD = CODE.parent
ROOT = REBUILD.parents[1]
MS_DIR = REBUILD / "manuscript"
TABLES = REBUILD / "data/source_data_rebuilt/supplementary_tables"
CLIN_SENS = REBUILD / "data/source_data_rebuilt/revision2_statistics/clinical_age_sex_wbc_sensitivity.tsv"
SUB = ROOT / "submission"
FIG = SUB / "figures"
OUT_MS = SUB / "manuscript"
OUT_SUPP = SUB / "supplementary"
OUT_CL = SUB / "cover_letter"
BUILD = SUB / "_build"
for d in (OUT_MS, OUT_SUPP, OUT_CL, BUILD):
    d.mkdir(parents=True, exist_ok=True)

PANDOC = shutil.which("pandoc")
if not PANDOC:
    raise SystemExit("pandoc not found")


def superscripts(text: str) -> str:
    return re.sub(r"<sup>(.*?)</sup>", lambda m: "^" + m.group(1).replace(" ", "\\ ") + "^", text)


def reference_docx() -> Path:
    ref = BUILD / "reference.docx"
    subprocess.run([PANDOC, "-o", str(ref), "--print-default-data-file", "reference.docx"],
                   check=True)
    doc = Document(ref)
    for name in ("Normal", "Body Text", "First Paragraph", "Compact"):
        if name not in [s.name for s in doc.styles]:
            continue
        st = doc.styles[name]
        st.font.name = "Times New Roman"
        st.font.size = Pt(12)
        st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
        pf = st.paragraph_format
        pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        pf.space_after = Pt(0)
        pf.space_before = Pt(0)
    for name, size in (("Title", 16), ("Heading 1", 14), ("Heading 2", 13), ("Heading 3", 12)):
        if name not in [s.name for s in doc.styles]:
            continue
        st = doc.styles[name]
        st.font.name = "Arial"
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = None
        st.element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
        st.paragraph_format.space_before = Pt(12)
        st.paragraph_format.space_after = Pt(6)
    if "Caption" in [s.name for s in doc.styles]:
        doc.styles["Caption"].font.italic = False
    doc.save(ref)
    return ref


def add_line_numbers(path: Path) -> None:
    doc = Document(path)
    for section in doc.sections:
        sect = section._sectPr
        ln = OxmlElement("w:lnNumType")
        ln.set(qn("w:countBy"), "1")
        ln.set(qn("w:restart"), "continuous")
        sect.append(ln)
        footer = section.footer.paragraphs[0]
        footer.alignment = 1
        run = footer.add_run()
        for tag, text in (("begin", None), (None, "PAGE"), ("end", None)):
            if tag:
                fld = OxmlElement("w:fldChar")
                fld.set(qn("w:fldCharType"), tag)
                run._r.append(fld)
            else:
                instr = OxmlElement("w:instrText")
                instr.set(qn("xml:space"), "preserve")
                instr.text = text
                run._r.append(instr)
    doc.save(path)


def to_docx(md_text: str, target: Path, ref: Path, line_numbers: bool) -> None:
    src = BUILD / (target.stem + ".md")
    src.write_text(md_text, encoding="utf-8")
    subprocess.run([PANDOC, str(src), "-f", "markdown+pipe_tables+superscript-raw_html",
                    "-t", "docx", "--reference-doc", str(ref),
                    "--resource-path", str(SUB), "-o", str(target)], check=True)
    if line_numbers:
        add_line_numbers(target)
    print(f"  {target.relative_to(ROOT)}  ({target.stat().st_size // 1024} KB)")


def word_counts(md: str) -> tuple[int, int]:
    abstract = md.split("## Abstract", 1)[1].split("\n## ", 1)[0]
    body = md.split("## Introduction", 1)[1].split("## Methods", 1)[0]
    clean = lambda s: len(re.sub(r"[#*\[\]()]", " ", s).split())
    return clean(abstract), clean(body)


# ── Main manuscript ─────────────────────────────────────────────────────────
main_md = (MS_DIR / "ZEB1_ZEB2_TALL_manuscript_polished.md").read_text(encoding="utf-8")
main_md = superscripts(main_md)
n_abs, n_body = word_counts(main_md)
n_refs = len(re.findall(r"^\d+\. ", main_md.split("## References", 1)[1], flags=re.M))
counts = (f"**Word count:** abstract {n_abs}; main text (Introduction to Discussion) {n_body}.  \n"
          f"**Figures:** 7. **References:** {n_refs}. **Supplementary material:** "
          f"7 figures and 4 tables.\n")
main_md = main_md.replace("\n---\n\n## Abstract", "\n" + counts + "\n---\n\n## Abstract", 1)

ref = reference_docx()
print("[manuscript]")
to_docx(main_md, OUT_MS / "ZEB1_ZEB2_TALL_manuscript.docx", ref, line_numbers=True)

with_figs = main_md
for i in range(1, 8):
    pat = re.compile(rf"(\*\*Figure {i}\..*?)(\n\n)", re.S)
    with_figs = pat.sub(lambda m, i=i: m.group(1) + f"\n\n![](figures/main/Figure{i}.png){{width=6.5in}}\n\n",
                        with_figs, count=1)
to_docx(with_figs, OUT_MS / "ZEB1_ZEB2_TALL_manuscript_with_figures.docx", ref, line_numbers=False)
(OUT_MS / "ZEB1_ZEB2_TALL_manuscript.md").write_text(main_md, encoding="utf-8")

# ── Supplementary Information ───────────────────────────────────────────────
si = (MS_DIR / "Supplementary_Information.md").read_text(encoding="utf-8")
si = si.split("## Supplementary Table S1.", 1)[0].rstrip() + "\n"
si = si.replace("Table S4", "Table S3")
si = re.sub(r" ?Table S3 records excluded exploratory lines[^.]*\.", "", si)
si = re.sub(r"the full-precision values, event counts and Cox warnings are in `[^`]+`",
            "full-precision values, event counts and Cox warnings are provided in "
            "Supplementary Table S4", si)
si = re.sub(r"\bfrozen\s+", "", si)
si = re.sub(r"\bFrozen\s+(\w)", lambda m: m.group(1).upper(), si)
for i in range(1, 8):
    pat = re.compile(rf"(## Supplementary Figure S{i}\..*?\n\n.*?)(\n\n|\Z)", re.S)
    si = pat.sub(lambda m, i=i: m.group(1) +
                 f"\n\n![](figures/supplementary/FigureS{i}.png){{width=6.5in}}\n\n",
                 si, count=1)
si += """
## Supplementary Tables

Supplementary Tables S1–S4 are provided as separate sheets of the file *Supplementary_Tables_S1-S4.xlsx*.

**Supplementary Table S1. Clinical models and subtype attenuation.** Six clinical endpoints with unadjusted and molecular-subtype-adjusted estimates, 95% confidence intervals, denominators, event counts, P values and Benjamini–Hochberg adjustment across endpoints within each adjustment level.

**Supplementary Table S2. Datasets, biological units and evidentiary roles.** Accession, assay, biological unit, sample size and the role of each dataset in the analysis.

**Supplementary Table S3. Full 17-subtype developmental residual summary.** Within-cohort residual medians with 3,000-replicate stratified patient-bootstrap 95% intervals, joint-model development-adjusted subtype effects with HC3 robust tests and Benjamini–Hochberg FDR, and the normal-stage projection sensitivity.

**Supplementary Table S4. Clinical sensitivity models additionally adjusted for age, sex and white-cell count.** Firth logistic (binary endpoints) and Efron-ties Cox (survival endpoints) estimates per one-unit higher ZEB balance, with endpoint-specific complete cases and model warnings.
"""
si = superscripts(si)
print("[supplementary]")
to_docx(si, OUT_SUPP / "ZEB1_ZEB2_TALL_Supplementary_Information.docx", ref, line_numbers=False)

# ── Supplementary tables workbook ───────────────────────────────────────────
wb = openpyxl.Workbook()
wb.remove(wb.active)
sheets = [
    ("Table S1", "Clinical models and subtype attenuation", TABLES / "Table_S1_clinical_models.tsv"),
    ("Table S2", "Datasets, biological units and evidentiary roles", TABLES / "Table_S2_dataset_units.tsv"),
    ("Table S3", "Full 17-subtype developmental residual summary", TABLES / "Table_S4_full_subtype_residuals.tsv"),
    ("Table S4", "Clinical sensitivity adjusted for age, sex and white-cell count", CLIN_SENS),
]
drop_cols = {"source", "source_path", "model_id", "source_file"}
for name, title, path in sheets:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    keep = [i for i, h in enumerate(rows[0]) if h.strip().lower() not in drop_cols]
    ws = wb.create_sheet(name)
    ws.append([f"Supplementary {name}. {title}"])
    ws["A1"].font = Font(bold=True, size=12)
    ws.append([])
    for r, row in enumerate(rows):
        vals = []
        for i in keep:
            v = row[i] if i < len(row) else ""
            try:
                v = float(v) if re.fullmatch(r"-?\d+(\.\d+)?([eE][+-]?\d+)?", v) else v
            except ValueError:
                pass
            vals.append(v)
        ws.append(vals)
        if r == 0:
            for c in ws[ws.max_row]:
                c.font = Font(bold=True)
                c.alignment = Alignment(wrap_text=True, vertical="top")
    for col in ws.columns:
        width = max(len(str(c.value)) if c.value is not None else 0 for c in col[2:])
        ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 48)
    ws.freeze_panes = "A4"
xlsx = OUT_SUPP / "Supplementary_Tables_S1-S4.xlsx"
wb.save(xlsx)
print(f"  {xlsx.relative_to(ROOT)}")

# ── Cover letter ────────────────────────────────────────────────────────────
cover = """**Date:** [date of submission]

To the Editor-in-Chief, *Haematologica*

Dear Editor,

We are pleased to submit our manuscript entitled **"Molecular subtypes of T-cell acute lymphoblastic leukemia differentially reconfigure a developmentally patterned ZEB1–ZEB2 axis"** for consideration as an Original Article in *Haematologica*.

Contemporary genomic studies now resolve T-cell acute lymphoblastic leukemia (T-ALL) into at least 17 driver-defined molecular subtypes. How much of a subtype's gene-expression configuration is inherited from its normal developmental position, and how much reflects the molecular subtype itself, has not been addressed systematically. We examined this question for the ZEB1–ZEB2 axis by combining four normal human thymus references with a 1,309-patient molecularly annotated diagnostic cohort, and tested the resulting framework in independent bulk RNA, patient-matched single-cell, lesion-defined and chromatin-contact data.

We believe three findings will interest your readership:

1. A ZEB-independent developmental coordinate explains 26% of ZEB1–ZEB2 balance variation in T-ALL, and molecular subtype contributes a further 14 percentage points after development (P = 7.8 × 10⁻⁴⁸). ETP-like leukemia, despite its low raw balance, lies close to its developmental expectation, whereas BCL11B, SPI1 and LMO2 γδ-like subtypes deviate strongly.
2. All 12 BCL11B-rearranged leukemias in a lesion-defined series converge on ZEB2-dominant expression, including seven cases in which ZEB2 is not the rearrangement partner.
3. Unadjusted clinical associations of ZEB balance with induction failure and minimal residual disease largely attenuate after molecular-subtype adjustment, indicating that subtype composition, rather than ZEB expression itself, carries most of the apparent prognostic signal.

The within-cohort developmental adjustment we describe separates inherited developmental position from subtype-associated reconfiguration without requiring cross-platform normalization, and is applicable to other lineage-associated expression programs in hematologic malignancies.

This manuscript has not been published and is not under consideration elsewhere. All authors have read and approved the submission and declare no competing interests. The study used only publicly available, de-identified data. Analysis code is available at https://github.com/JDjust/ZEB1-TALL-analysis.

Suggested reviewers: [name, institution, e-mail — three experts in T-ALL genomics or human thymopoiesis, without recent collaboration with the authors].

Thank you for considering our work.

Sincerely,

Chao Wu, on behalf of all authors  
Department of Pulmonary Oncology, Tianjin Medical University Cancer Institute & Hospital, Tianjin, China  
chaowutjmuch@163.com

Limei Li  
Department of Hematology, The Second Affiliated Hospital of Hainan Medical University, Haikou, China  
lilimei116@126.com
"""
print("[cover letter]")
to_docx(cover, OUT_CL / "Cover_letter_Haematologica.docx", ref, line_numbers=False)
print("Word counts: abstract", n_abs, "main text", n_body, "references", n_refs)
