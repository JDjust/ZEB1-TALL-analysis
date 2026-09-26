"""Typeset the frozen subtype residual table as a one-page landscape supplement."""
from pathlib import Path
import csv
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer

BASE = Path(__file__).resolve().parents[1]
src = BASE / "data/source_data_rebuilt/supplementary_tables/Table_S4_full_subtype_residuals.tsv"
dst = BASE / "manuscript/Supplementary_Table_S4.pdf"
with src.open(encoding="utf-8-sig", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

def f(value):
    return f"{float(value):.2f}"

def p(value):
    x = float(value)
    return f"{x:.1e}" if x < 0.001 else f"{x:.3f}"

headers = ["Subtype", "n", "Dev. median", "Balance median", "Normal residual", "Cohort residual", "Outside range", "P", "FDR"]
data = [headers]
for r in rows:
    subtype = r["subtype"].replace("γδ", "gamma-delta").replace("αβ", "alpha-beta")
    data.append([subtype, r["n"], f(r["dev_median"]), f(r["balance_median"]),
                 f(r["residual_normal_median"]), f(r["residual_within_median"]),
                 r["n_outside_normal_dev_range"], p(r["wilcox_p_residual_normal"]),
                 p(r["wilcox_fdr"])])

doc = SimpleDocTemplate(str(dst), pagesize=landscape(A4), leftMargin=28, rightMargin=28,
                        topMargin=26, bottomMargin=25)
title_style = ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=13, leading=16, spaceAfter=5)
note_style = ParagraphStyle("note", fontName="Helvetica", fontSize=7.5, leading=10, spaceBefore=8)
table = Table(data, colWidths=[119, 25, 67, 75, 79, 79, 72, 60, 60], rowHeights=[25] + [22]*len(rows), repeatRows=1)
table.setStyle(TableStyle([
    ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#21384b")),
    ("TEXTCOLOR", (0,0), (-1,0), colors.white),
    ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
    ("FONTNAME", (0,1), (-1,-1), "Helvetica"),
    ("FONTSIZE", (0,0), (-1,-1), 8),
    ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ("ALIGN", (1,0), (-1,-1), "CENTER"),
    ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f1f6f8")]),
    ("LINEBELOW", (0,0), (-1,0), .7, colors.HexColor("#21384b")),
    ("BOTTOMPADDING", (0,0), (-1,0), 7),
    ("TOPPADDING", (0,0), (-1,0), 7),
]))
story = [Paragraph("Supplementary Table S4. Developmental residuals across 17 T-ALL molecular subtypes", title_style), table,
         Paragraph("Values are frozen subtype medians. Normal residual: observed balance minus the normal-stage spline expectation. Cohort residual: sensitivity from the within-cohort spline. Outside range counts patients beyond the measured normal-stage coordinate domain. P and FDR are frozen Wilcoxon results for the normal-reference residual. Statistical significance alone does not imply a large subtype-specific displacement: ETP-like remains predominantly developmentally inherited, whereas BCL11B is extreme but substantially extrapolated (14/18 outside range). Full-precision values and frozen call labels are retained in the accompanying TSV source table.", note_style)]
doc.build(story)
print(dst)
