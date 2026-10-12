"""Read actual Office PDFs and refresh final visual inspection objects."""
from pathlib import Path
import json, hashlib
import fitz
from PIL import Image, ImageDraw, ImageChops
import sys

B = Path(__file__).resolve().parent.parent
Q = B / 'QA'
D = Q / 'document_pages'
D.mkdir(exist_ok=True)
reports = []
for name in ['MANUSCRIPT', 'SUPPLEMENTARY_INFORMATION']:
    path = B / 'manuscript' / (name + '.pdf')
    with fitz.open(path) as doc:
        pages = []
        thumbs = []
        for i, page in enumerate(doc):
            text = page.get_text()
            spans = [s for b in page.get_text('dict')['blocks'] if b['type'] == 0
                     for line in b['lines'] for s in line['spans']]
            outside = [s['text'] for s in spans if s['bbox'][0] < -.5 or s['bbox'][1] < -.5
                       or s['bbox'][2] > page.rect.width + .5 or s['bbox'][3] > page.rect.height + .5]
            unresolved = bool(__import__('re').search(r'@[A-Za-z][A-Za-z0-9]+', text))
            eq_terms = [k for k in ['Φ', 'ρ', 'β', 'Δ', 'epsilon', 'ε', 'α'] if k in text]
            rec = dict(page=i+1, characters=len(text), outside_page_text=outside,
                       unresolved_citation_key=unresolved, equation_symbols=eq_terms)
            pages.append(rec)
            # 144 dpi page images support actual equation/reference inspection.
            pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
            pix.save(D / f'{name}_{i+1:02d}.png')
            image = Image.frombytes('RGB', [pix.width, pix.height], pix.samples)
            image.thumbnail((400, 535))
            thumbs.append(image)
        cols, cellw, cellh = 4, 420, 565
        sheet = Image.new('RGB', (cols*cellw, ((len(thumbs)+cols-1)//cols)*cellh), '#EAEAEA')
        draw = ImageDraw.Draw(sheet)
        for i, im in enumerate(thumbs):
            x, y = i%cols*cellw, i//cols*cellh
            sheet.paste(im, (x+(cellw-im.width)//2, y+25))
            draw.text((x+12, y+7), f'{name} p{i+1}', fill='black')
        sheet.save(Q / (name + '_PAGES_CONTACT_SHEET.png'))
        status = 'PASS' if all(p['characters']>80 and not p['outside_page_text']
                               and not p['unresolved_citation_key'] for p in pages) else 'FAIL'
        reports.append(dict(document=name, status=status, pages=len(doc),
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest(), page_records=pages,
                            scope='PDF bounds, text/citation readback; human visual review recorded separately'))
(Q/'DOCUMENT_PDF_QA.json').write_text(json.dumps(reports, indent=2), encoding='utf-8')

if '--documents-only' in sys.argv:
    print(json.dumps([dict(document=r['document'],status=r['status'],pages=r['pages']) for r in reports]))
    assert all(r['status']=='PASS' for r in reports), 'Document PDF check failed'
    sys.exit(0)

names = [f'Fig{i}' for i in range(1,7)] + [f'SFig{i}' for i in range(1,7)]
edits = []
sheet = Image.new('RGB', (2700, 2520), '#EAEAEA')
draw = ImageDraw.Draw(sheet)
for i, name in enumerate(names):
    orig, edited = [Q / (name+suffix) for suffix in ['_office_original.png','_office_edit.png']]
    a, b = [Image.open(p).convert('RGB') for p in [orig,edited]]
    box = ImageChops.difference(a,b).getbbox()
    assert box, name
    edits.append(dict(figure=name, persisted_edit_changes_pixels=True, difference_bbox=box,
                      original_export_sha256=hashlib.sha256(orig.read_bytes()).hexdigest(),
                      edited_export_sha256=hashlib.sha256(edited.read_bytes()).hexdigest()))
    a.thumbnail((884,592))
    x, y = i%3*900, i//3*630
    sheet.paste(a, (x+(900-a.width)//2, y+28))
    draw.text((x+12,y+8),name,fill='black')
sheet.save(Q/'OFFICE_COMPLETE_CONTACT_SHEET.png')
(Q/'OFFICE_PIXEL_EDIT_QA.json').write_text(json.dumps(dict(status='PASS', figures=edits),indent=2),encoding='utf-8')
print(json.dumps([dict(document=r['document'],status=r['status'],pages=r['pages']) for r in reports]))
assert all(r['status']=='PASS' for r in reports), 'Document PDF check failed'
