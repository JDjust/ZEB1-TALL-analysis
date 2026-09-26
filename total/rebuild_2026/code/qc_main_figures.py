"""Build a review-only contact sheet and verify canonical main-figure artifacts."""
from pathlib import Path
import csv
import fitz
from PIL import Image, ImageDraw, ImageFont, ImageOps

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "edition_20260925" / "figures" / "main"
REVIEW = ROOT / "figures" / "exports_for_review"
LOG = ROOT / "logs"
REVIEW.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

rows = []
canvas = Image.new("RGB", (2600, 3900), "white")
draw = ImageDraw.Draw(canvas)
font = ImageFont.truetype("arial.ttf", 36)
for i in range(1, 7):
    tag = f"F{i}"
    png, pdf, svg, tiff = (MAIN / f"{tag}.{ext}" for ext in
                           ("png", "pdf", "svg", "tiff"))
    for path in (png, pdf, svg, tiff):
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty main figure artifact: {path}")
    doc = fitz.open(pdf)
    if len(doc) != 1:
        raise RuntimeError(f"Expected one PDF page for {tag}: {len(doc)}")
    with Image.open(png) as original:
        tile = ImageOps.contain(original.convert("RGB"), (1260, 1190))
    col, row = (i - 1) % 2, (i - 1) // 2
    x, y = 30 + 1300 * col, 65 + 1300 * row
    canvas.paste(tile, (x + (1260 - tile.width) // 2,
                        y + (1190 - tile.height) // 2))
    draw.text((x, y - 55), tag, font=font, fill="#152b39")
    rows.append({"figure": tag, "pdf_pages": len(doc),
                 "pdf_width_pt": round(doc[0].rect.width, 1),
                 "pdf_height_pt": round(doc[0].rect.height, 1),
                 "png_width_px": original.width,
                 "png_height_px": original.height,
                 "pdf_text_words": len(doc[0].get_text().split()),
                 "source_data_files": len(list((ROOT / "data" /
                                         "source_data_rebuilt" / tag).glob("*.tsv"))),
                 "legend_present": (ROOT / "manuscript" /
                                    f"Figure_{i}_legend.md").exists()})
    doc.close()

sheet = REVIEW / "main_figures_contact_sheet.png"
canvas.save(sheet, optimize=True)
with (LOG / "main_figure_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
print(sheet)
