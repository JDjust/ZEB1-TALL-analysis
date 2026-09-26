"""Typeset the frozen clinical models as a one-page landscape supplement."""
from pathlib import Path
import csv
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph

BASE = Path(__file__).resolve().parents[1]
src = BASE / "data/source_data_rebuilt/supplementary_tables/Table_S1_clinical_models.tsv"
dst = BASE / "manuscript/Supplementary_Table_S1.pdf"
with src.open(encoding="utf-8-sig", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

def fmt_p(v):
    x = float(v)
    return f"{x:.1e}" if x < .001 else f"{x:.3f}"

data = [["Endpoint", "Scale", "Adjustment", "n", "Events", "Estimate", "95% CI", "P"]]
for r in rows:
    data.append([r["endpoint"].replace(">=", "≥"), r["measure"], r["adjustment"], r["n"], r["events"],
                 f'{float(r["estimate"]):.2f}',
                 f'{float(r["ci_low"]):.2f}–{float(r["ci_high"]):.2f}', fmt_p(r["p"])])

doc = SimpleDocTemplate(str(dst), pagesize=landscape(A4), leftMargin=35, rightMargin=35,
                        topMargin=35, bottomMargin=30)
title = ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=14, leading=17, spaceAfter=12)
note = ParagraphStyle("note", fontName="Helvetica", fontSize=8.2, leading=11, spaceBefore=11)
table = Table(data, colWidths=[145, 38, 115, 45, 52, 68, 95, 65], rowHeights=[29]+[28]*len(rows), repeatRows=1)
table.setStyle(TableStyle([
    ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#21384b")),
    ("TEXTCOLOR", (0,0), (-1,0), colors.white),
    ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
    ("FONTNAME", (0,1), (-1,-1), "Helvetica"),
    ("FONTSIZE", (0,0), (-1,-1), 9),
    ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ("ALIGN", (1,0), (-1,-1), "CENTER"),
    ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f1f6f8")]),
    ("LINEBELOW", (0,0), (-1,0), .7, colors.HexColor("#21384b")),
]))
story = [Paragraph("Supplementary Table S1. Frozen clinical associations before and after molecular-subtype adjustment", title),
         table, Paragraph("Estimates are odds ratios (OR) or hazard ratios (HR) per unit balance from the frozen model files. Molecular-subtype adjustment attenuates the induction-failure, MRD ≥0.1%, event-free-survival and overall-survival associations. M2/M3 morphology and MRD ≥0.01% retain adjusted associations; no result is presented as an independent prognostic tool. Full model identifiers, exact values and source paths are provided in the accompanying TSV table.", note)]
doc.build(story)
print(dst)
