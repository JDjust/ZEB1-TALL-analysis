"""Render the reader-facing supplementary workbook as readable PDF table parts.

Wide source tables are split by columns, repeating identifying columns. The
workbook remains the full-precision source; PDF numbers are rounded for display.
"""
from html import escape
from pathlib import Path

import openpyxl
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "submission/supplementary/Supplementary_Tables_S1-S8.xlsx"
OUTPUT = ROOT / "submission/first_submission/Supplementary_Tables_S1-S8_readable.pdf"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
pdfmetrics.registerFont(TTFont("Arial", "C:/Windows/Fonts/arial.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Bold", "C:/Windows/Fonts/arialbd.ttf"))

# Zero-based source columns. Repeated columns identify the same records.
PARTS = {
    "Table S1": [[0, 1, 2, 3, 4, 5], [0, 1, 2, 6, 7, 8], [0, 1, 2, 9, 10, 11]],
    "Table S2": [[0, 1, 2, 3, 4, 5]],
    "Table S3": [[0, 1, 2, 3, 4, 5], [0, 1, 6, 7, 8, 9]],
    "Table S4": [[0, 1, 2, 3, 4, 5, 6, 7, 8]],
    "Table S5": [[0, 1, 2, 3, 4], [0, 5, 6], [0, 7, 8, 9, 10], [0, 11]],
    "Table S6": [[0, 1, 2, 3, 4], [0, 5, 6, 7, 8, 9]],
    "Table S7 HiChIP": [[0, 1, 2, 3, 4, 5, 6, 7], [0, 1, 8, 9]],
    "Table S7 scATAC": [[0, 1, 2, 3, 4, 5], [0, 1, 6, 7, 8, 9]],
    "Table S8 Crossfit": [[0, 1, 2, 3, 4]],
    "Table S8 Scores": [[0, 1, 2, 3, 4, 5]],
    "Table S8 Models": [[0, 1, 2, 3, 4, 5], [0, 1, 6, 7, 8]],
    "Table S8 Contrasts": [[0, 1, 2, 3, 4, 5]],
}


def display(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value == 0:
            return "0"
        if abs(value) < 0.001:
            return f"{value:.2e}"
        return f"{value:.3f}"
    return str(value)


styles = {
    "title": ParagraphStyle("title", fontName="Arial-Bold", fontSize=11,
                            leading=14, textColor=colors.HexColor("#203744")),
    "note": ParagraphStyle("note", fontName="Arial", fontSize=7.5,
                           leading=9.2, textColor=colors.HexColor("#51636b")),
    "header": ParagraphStyle("header", fontName="Arial-Bold", fontSize=8,
                             leading=9.5, alignment=TA_CENTER,
                             textColor=colors.white),
    "cell": ParagraphStyle("cell", fontName="Arial", fontSize=8,
                           leading=10, alignment=TA_LEFT,
                           textColor=colors.HexColor("#203744")),
}


def para(value, style, header=False):
    value = display(value)
    if header:
        value = value.replace("_", "_<br/>")
    return Paragraph(escape(value).replace("&lt;br/&gt;", "<br/>"), style)


def widths(headers, rows, columns, available):
    weights = []
    for col in columns:
        samples = [display(r[col]) for r in rows]
        longest = max([len(headers[col])] + [len(x) for x in samples])
        weight = min(26, max(7, longest))
        if col == 0:
            weight = max(weight, 13)
        weights.append(weight)
    total = sum(weights)
    return [available * x / total for x in weights]


book = openpyxl.load_workbook(SOURCE, data_only=True, read_only=True)
page_w, page_h = landscape(A4)
doc = SimpleDocTemplate(str(OUTPUT), pagesize=(page_w, page_h),
                        leftMargin=20 * mm, rightMargin=20 * mm,
                        topMargin=15 * mm, bottomMargin=16 * mm)
available = page_w - 40 * mm
story = []

for sheet in book:
    raw = list(sheet.values)
    title = str(raw[0][0])
    note = str(raw[1][0])
    headers = list(raw[2])
    rows = [list(r) for r in raw[3:] if any(v is not None for v in r)]
    groups = PARTS[sheet.title]
    for part_no, columns in enumerate(groups, start=1):
        if story:
            story.append(Spacer(1, 8 * mm))
        suffix = f" (part {part_no} of {len(groups)})" if len(groups) > 1 else ""
        title_flow = para(title + suffix, styles["title"])
        note_flow = para(note + " Exact values remain in source TSV files.", styles["note"])
        data = [[para(headers[c], styles["header"], header=True) for c in columns]]
        data += [[para(r[c], styles["cell"]) for c in columns] for r in rows]
        table = Table(data, colWidths=widths(headers, rows, columns, available),
                      repeatRows=1, hAlign="LEFT")
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#244455")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1),
             [colors.white, colors.HexColor("#F3F7F8")]),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ("LINEBELOW", (0, 0), (-1, 0), .45, colors.HexColor("#244455")),
        ]))
        story.append(KeepTogether([title_flow, note_flow, table]))

doc.build(story)
print(OUTPUT)
