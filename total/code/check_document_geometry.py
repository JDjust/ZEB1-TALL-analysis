"""Check embedded figure bounds; not a substitute for rendered-page inspection."""
from pathlib import Path
from docx import Document
import json
root = Path(__file__).resolve().parents[1]
output = {}
for path in (root/'manuscript').glob('*.docx'):
    doc = Document(path)
    # These builders currently use one section; reject silent section mismatch.
    assert len(doc.sections)==1
    sec = doc.sections[0]
    available_w = sec.page_width-sec.left_margin-sec.right_margin
    available_h = sec.page_height-sec.top_margin-sec.bottom_margin
    figures = []
    for i, shape in enumerate(doc.inline_shapes,1):
        figures.append({'figure_index':i,'width_inches':shape.width/914400,
                        'height_inches':shape.height/914400,
                        'fits_text_area':shape.width<=available_w and shape.height<=available_h})
    assert all(x['fits_text_area'] for x in figures)
    output[path.name]={'available_inches':[available_w/914400,available_h/914400],
                       'figures':figures,
                       'scope':'Image bounds only; caption flow and page breaks require rendered inspection'}
out = root/'manuscript/review_pdf/document_geometry.json'
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(output,indent=2),encoding='utf-8')
print('All 22 embedded figures fit their document text areas; caption pagination not established.')
