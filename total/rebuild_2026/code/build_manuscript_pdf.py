"""Build printable manuscript PDFs using local Pandoc and headless Chrome."""
from pathlib import Path
import shutil
import subprocess
import time

base = Path(__file__).resolve().parents[1]
man = base / "manuscript"
pandoc = shutil.which("pandoc")
chrome = shutil.which("chrome") or r"C:\Program Files\Google\Chrome\Application\chrome.exe"
if not pandoc or not Path(chrome).exists():
    raise SystemExit("Pandoc and Google Chrome are required for PDF export")

for stem, md_name in (("ZEB1_ZEB2_TALL_manuscript", "ZEB1_ZEB2_TALL_manuscript.md"),
                      ("Supplementary_Information", "Supplementary_Information.md")):
    html = man / f"{stem}.html"
    pdf = man / f"{stem}.pdf"
    css = "supplement_print.css" if stem == "Supplementary_Information" else "print.css"
    subprocess.run([pandoc, str(man / md_name), "-f", "gfm", "-t", "html", "-s",
                    "--css", css, "-o", str(html)], check=True)
    profile = base / "logs" / f"chrome_pdf_profile_{stem}"
    subprocess.run([str(chrome), "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                    f"--user-data-dir={profile}", f"--print-to-pdf={pdf}", html.as_uri()],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(20):
        if pdf.exists() and pdf.stat().st_size > 10000:
            break
        time.sleep(.25)
    print(pdf.name, pdf.stat().st_size)
    html.unlink(missing_ok=True)
