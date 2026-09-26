"""Print manuscript text and append native vector figures and reader tables."""
from pathlib import Path
import re
import subprocess
import tempfile

import fitz

ROOT=Path(__file__).resolve().parents[3]
SUB=ROOT/"submission"
FIRST=SUB/"first_submission"
BUILD=SUB/"_build"
PANDOC=r"C:\Users\Ls180\AppData\Local\Microsoft\WinGet\Packages\JohnMacFarlane.Pandoc_Microsoft.Winget.Source_8wekyb3d8bbwe\pandoc-3.10\pandoc.exe"
CHROME=r"C:\Program Files\Google\Chrome\Application\chrome.exe"
CSS="""
@page { size: A4; margin: 19mm 20mm 18mm 20mm; }
body { font-family: 'Times New Roman', serif; font-size: 11pt; line-height: 1.48; color: #17252d; }
h1,h2,h3 { font-family: Arial, sans-serif; break-after: avoid; color: #203744; }
h1 { font-size: 17pt; } h2 { font-size: 13pt; margin-top: 1.25em; }
h3 { font-size: 11.2pt; margin-top: 1em; }
p { margin: .22em 0 .56em; orphans: 3; widows: 3; }
p:has(+ table) { break-before: page; break-after: avoid; }
p:has(+ p > math[display="block"]) { break-after: avoid; }
p:has(> math[display="block"]) { break-inside: avoid; }
table { width: 100%; border-collapse: collapse; font-family: Arial, sans-serif;
        font-size: 7.7pt; line-height: 1.27; break-inside: avoid; }
thead { display: table-header-group; background: #244455; color: white; }
tr { break-inside: avoid; }
.math.display { display: block; text-align: center; margin: .75em 0 1em; }
math[display="block"] { display: block; text-align: center; margin: .75em 0 1em; }
th,td { padding: 4px 5px; border-bottom: .4px solid #cbd5d8; vertical-align: top; }
td:first-child { font-weight: 600; }
blockquote { margin-left: 0; border-left: 2px solid #9fb5bd; padding-left: 9px; }
"""

def print_md(md,out,stem):
    source=BUILD/(stem+"_vector_print.md")
    html=BUILD/(stem+"_vector_print.html")
    source.write_text(md,encoding="utf-8")
    subprocess.run([PANDOC,str(source),"-f","markdown+pipe_tables+superscript-raw_html",
                    "-t","html5","-s","--mathml","-o",str(html)],check=True)
    s=html.read_text(encoding="utf-8")
    s=s.replace("</head>","<style>"+CSS+"</style></head>",1)
    html.write_text(s,encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="zeb_preview_",dir=BUILD) as profile:
        subprocess.run([CHROME,"--headless=new","--disable-gpu","--no-first-run",
                        "--no-pdf-header-footer",f"--user-data-dir={profile}",
                        f"--print-to-pdf={out}",html.as_uri()],
                       check=True,timeout=180,
                       stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    if not out.exists() or out.stat().st_size<10000:
        raise RuntimeError(f"Chrome did not render {out}")

def combine(parts,out):
    result=fitz.open()
    for p in parts:
        with fitz.open(p) as one:
            result.insert_pdf(one)
    result.save(out,garbage=4,deflate=True)
    print(out,len(result),"pages",out.stat().st_size,"bytes")
    result.close()

FIRST.mkdir(parents=True,exist_ok=True)
main=(SUB/"manuscript/ZEB1_ZEB2_TALL_manuscript.md").read_text(encoding="utf-8")
main_text=BUILD/"Main_Article_text_vector.pdf"
print_md(main,main_text,"main")
combine([main_text]+[SUB/f"figures/main/Figure{i}.pdf" for i in range(1,8)],
        FIRST/"Main_Article_with_Figures.pdf")

si=(BUILD/"ZEB1_ZEB2_TALL_Supplementary_Information_editorial_rebuild.md").read_text(encoding="utf-8")
si=re.sub(r"^!\[\]\(figures/supplementary/[^\n]+\)\{[^\n]+\}\s*$","",si,flags=re.M)
si_text=BUILD/"Supplementary_Information_text_vector.pdf"
print_md(si,si_text,"supplement")
combine([si_text]+[SUB/f"figures/supplementary/FigureS{i}.pdf" for i in range(1,10)]
        +[FIRST/"Supplementary_Tables_S1-S8_readable.pdf"],
        FIRST/"Supplementary_Information_complete.pdf")
