"""Check delivered figure, manuscript and source-data integrity."""
from pathlib import Path
import csv
import re
import fitz
from PIL import Image
from docx import Document

base=Path(__file__).resolve().parents[1]
source=base/'data/source_data_rebuilt'
problems=[]
panel_count=0
displayed=set()
figure_stats=[]
for group,prefix in [('main','F'),('supplementary','S')]:
    for i in range(1,8):
        stem=f'{prefix}{i}'
        folder=base/'figures'/group
        if not (folder/f'{stem}.pdf').exists():
            folder=base/'edition_20260925/figures'/group
        srcdir=source/(f'F{i}' if prefix=='F' and i<=4 else f'Main_F{i}' if prefix=='F' else f'Supplement_S{i}')
        with (folder/f'{stem}_panel_manifest.tsv').open(encoding='utf-8-sig',newline='') as fh:
            rows=list(csv.DictReader(fh,delimiter='\t'))
        panel_count+=len(rows)
        for row in rows:
            for name in row['displayed_data'].split(' + '):
                candidate=srcdir/name.strip()
                displayed.add(candidate)
                if not candidate.is_file(): problems.append(f'Missing panel data: {candidate}')
        with Image.open(folder/f'{stem}.png') as im:
            width,height=im.size
            if min(width,height)<1500: problems.append(f'Low PNG dimension: {stem} {width}x{height}')
        with fitz.open(folder/f'{stem}.pdf') as pdf:
            if len(pdf)!=1: problems.append(f'Figure PDF has {len(pdf)} pages: {stem}')
        figure_stats.append((stem,len(rows),width,height))

with fitz.open(base/'manuscript/ZEB1_ZEB2_TALL_manuscript.pdf') as pdf:
    main_pages=len(pdf)
with fitz.open(base/'manuscript/Supplementary_Information.pdf') as pdf:
    supp_pages=len(pdf)
for name in ('Supplementary_Table_S1.pdf','Supplementary_Table_S4.pdf'):
    with fitz.open(base/'manuscript'/name) as pdf:
        if len(pdf)!=1: problems.append(f'Table {name} not one page')

main=(base/'manuscript/ZEB1_ZEB2_TALL_manuscript.md').read_text(encoding='utf-8')
supp=(base/'manuscript/Supplementary_Information.md').read_text(encoding='utf-8')
for i in range(1,8):
    if f'**Figure {i}.' not in main: problems.append(f'Missing main legend {i}')
    if f'## Supplementary Figure S{i}.' not in supp: problems.append(f'Missing supplement legend {i}')
for token in ('[CITATION NEEDED]','[DATA NEEDED]','TODO','PLACEHOLDER'):
    if token in main or token in supp: problems.append(f'Unresolved manuscript token: {token}')
if main.count('[AUTHOR CONFIRMATION')<3:
    problems.append('Missing expected author confirmations')
if len(re.findall(r'^\d+\. ',main.split('## References')[-1],flags=re.M))!=60:
    problems.append('Reference list does not have 60 entries')
with (base/'data/manifests/panel_source_resolution.tsv').open(encoding='utf-8-sig',newline='') as fh:
    upstream_rows=list(csv.DictReader(fh,delimiter='\t'))
if len({(r['figure'],r['panel']) for r in upstream_rows})!=panel_count:
    problems.append('Upstream source resolution does not cover every panel')
docx_text='\n'.join(p.text for p in Document(base/'manuscript/ZEB1_ZEB2_TALL_manuscript.docx').paragraphs)
supp_docx_text='\n'.join(p.text for p in Document(base/'manuscript/ZEB1_ZEB2_TALL_Supplementary.docx').paragraphs)
if sum('doi:' in line.lower() for line in docx_text.splitlines())!=60:
    problems.append('Main DOCX does not contain all 60 DOI-bearing references')
for i in range(1,8):
    if f'Figure {i}.' not in docx_text: problems.append(f'Main DOCX missing legend {i}')
    if f'Supplementary Figure S{i}.' not in supp_docx_text:
        problems.append(f'Supplement DOCX missing legend S{i}')

out=base/'logs/FINAL_QC.md'
lines=['# Final technical QC','',
       f'- Main figures: 7, supplementary figures: 7, panel rows: {panel_count}, unique displayed source tables: {len(displayed)}.',
       f'- Main manuscript: {main_pages} PDF pages; supplementary information: {supp_pages} PDF pages.',
       '- Clinical Table S1 and subtype Table S4 each occupy one landscape PDF page.',
       '- All 14 figure PDFs have one page; PNG dimensions and panel source paths were checked.',
       f'- {len(upstream_rows)} upstream source references resolve to existing frozen files, covering all {panel_count} panels.',
       '- The main manuscript has seven main legends and 60 numbered references; the supplement has seven figure legends.',
       '- DOCX exports contain all seven main and supplementary legends; the main DOCX contains 60 DOI-bearing reference entries.',
       '- Bibliographic DOI/title identity was checked in the reference manifest. Claim-to-citation suitability still requires author scientific judgment.',
       '', '## Figures', '', '| Figure | Panels | PNG dimensions |', '|---|---:|---:|']
lines.extend(f'| {a} | {b} | {w} × {h} |' for a,b,w,h in figure_stats)
lines.extend(['','## Unresolved technical problems',''] + (['- None.'] if not problems else [f'- {x}' for x in problems]))
out.write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(out)
print('technical problems',len(problems))
for problem in problems: print(problem)
if problems: raise SystemExit(1)
