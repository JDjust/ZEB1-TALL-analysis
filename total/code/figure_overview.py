"""Build naturally ordered contact sheets for whole-package visual review."""
from pathlib import Path
import re
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]


def sheet(folder, pattern, output, columns=3):
    files = sorted(folder.glob(pattern), key=lambda p: int(re.search(r"\d+", p.name)[0]))
    width, height = 600, 530
    canvas = Image.new("RGB", (width * columns, height * ((len(files) + columns - 1) // columns)), "white")
    draw = ImageDraw.Draw(canvas)
    for i, path in enumerate(files):
        x, y = (i % columns) * width, (i // columns) * height
        with Image.open(path) as original:
            thumb = original.convert("RGB")
            thumb.thumbnail((width - 12, height - 32))
            canvas.paste(thumb, (x + (width - thumb.width) // 2, y + 25))
        draw.text((x + 8, y + 6), path.stem, fill="black")
    canvas.save(output)
    print(f"{output.name}: {len(files)} figures")


if __name__ == "__main__":
    sheet(ROOT / "figures/png", "Figure*.png", ROOT / "figures/main_overview.png")
    sheet(ROOT / "figures/Supplementary/png", "FigureS*.png", ROOT / "figures/supplementary_overview.png")
