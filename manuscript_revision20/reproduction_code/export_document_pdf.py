"""Export final, citation-resolved HTML using installed Chromium and native MathML.

This is a documented fallback after bounded Word PDF-export failures; DOCX
remains the editable manuscript. No scientific values are recomputed.
"""
from pathlib import Path
import json, hashlib, importlib.metadata, argparse
from playwright.sync_api import sync_playwright

B=Path(__file__).resolve().parent.parent
reports=[]
parser=argparse.ArgumentParser()
parser.add_argument('--documents',nargs='+',choices=['MANUSCRIPT','SUPPLEMENTARY_INFORMATION'],default=['MANUSCRIPT','SUPPLEMENTARY_INFORMATION'])
selected=parser.parse_args().documents
with sync_playwright() as p:
    browser=p.chromium.launch(channel='msedge',headless=True)
    for name in selected:
        source=B/'manuscript'/f'{name}.html'
        page=browser.new_page()
        page.goto(source.as_uri(),wait_until='load',timeout=60000)
        page.add_style_tag(content='''
          @page { size:A4; margin:20mm 21mm 20mm 21mm; }
          html { background:white; }
          body { max-width:none; margin:0; padding:0; color:#000;
                 font-family:Cambria,"Times New Roman",serif; font-size:11pt; line-height:1.22; }
          h1 { font-size:16pt; } h2 { font-size:13pt; } h3 { font-size:11pt; }
          h1,h2,h3 { break-after:avoid; }
          h2#abstract { break-before:page; page-break-before:always; }
          p,li { orphans:3; widows:3; } p { margin:0 0 6pt; }
          math { font-family:"Cambria Math",serif; } math[display=block] { margin:12pt 0; }
          a { color:inherit; text-decoration:none; } pre { white-space:pre-wrap; }
          .csl-entry { break-inside:avoid; margin-bottom:6pt; }
          table { border-collapse:collapse; font-size:9pt; } img { max-width:100%; }
        ''')
        page.evaluate('document.fonts.ready')
        pdf=B/'manuscript'/f'{name}.pdf'
        page.pdf(path=str(pdf),format='A4',print_background=True,prefer_css_page_size=True,
                 display_header_footer=True,header_template='<span></span>',
                 footer_template='<div style="font-size:8pt;width:100%;text-align:right;padding-right:21mm;color:black"><span class="pageNumber"></span></div>')
        reports.append(dict(document=name,engine='Microsoft Edge Chromium print',
            browser_version=browser.version,playwright=importlib.metadata.version('playwright'),
            html_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
            reason='Word opened actual DOCX but PDF export stalled; two bounded export attempts terminated only identified own processes.',
            scope='HTML and native MathML print fallback; this PDF is not represented as Word-exported'))
        page.close()
    browser.close()
receipt=B/'QA/DOCUMENT_PDF_ENGINE.json'
if receipt.exists():
    prior=json.loads(receipt.read_text(encoding='utf-8'))
    reports += [r for r in prior if r['document'] not in selected]
receipt.write_text(json.dumps(reports,indent=2),encoding='utf-8')
print(json.dumps([dict(document=r['document'],engine=r['engine']) for r in reports]))
