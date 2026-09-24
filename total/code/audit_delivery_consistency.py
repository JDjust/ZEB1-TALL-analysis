"""Check current delivery identity and numbering; not a scientific validity certificate."""
from pathlib import Path
import hashlib,json,re
from collections import Counter
from docx import Document
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reproducibility/current_delivery_checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
result={'scope':'Current image identity, figure/reference numbering, source hashes and Word pagination freshness', 'documents':{}}
for name,folder,n in [('ZEB1_TALL_manuscript','figures/png',7),('ZEB1_TALL_Supplementary','figures/Supplementary/png',15)]:
    path=ROOT/'manuscript'/(name+'.docx'); doc=Document(path)
    embedded=Counter(hashlib.sha256(doc.part.related_parts[s._inline.graphic.graphicData.pic.blipFill.blip.embed].blob).hexdigest() for s in doc.inline_shapes)
    files=list((ROOT/folder).glob('Figure*.png'))
    expected=Counter(sha(p) for p in files)
    assert len(files)==len(doc.inline_shapes)==n
    assert embedded==expected, 'Embedded images differ from current PNGs: '+name
    text='\n'.join(p.text for p in doc.paragraphs)
    prefix='Figure' if n==7 else 'Supplementary Figure'
    captions=[int(x) for x in re.findall(r'^'+prefix+r' (\d+)\.',text,re.M)]
    assert captions==list(range(1,n+1)),captions
    markers=re.findall(r'\[(?:AUTHOR CONFIRMATION|AUTHOR NEEDED|DATA NEEDED|CITATION NEEDED|METHOD CONFIRMATION)[^\]]*\]',text)
    assert not any(x.startswith(('[DATA NEEDED','[CITATION NEEDED','[AUTHOR NEEDED')) for x in markers)
    result['documents'][name]={'sha256':sha(path),'embedded_figures_match':n,'caption_numbers':captions,'remaining_markers':markers}
    if n==7:
        body, refs=text.split('\nReferences\n',1)
        numbers=[int(x) for x in re.findall(r'^\[(\d+)\]',refs,re.M)]
        catalog=json.loads((ROOT/'data/reference_catalog.json').read_text(encoding='utf-8'))['records']
        assert numbers==list(range(1,len(catalog)+1))
        used=set()
        for token in re.findall(r'\[([\d,;\s\-–]+)\]',body):
            for part in re.split(r'[,;]',token):
                limits=re.split('[-–]',part.strip())
                if len(limits)==2:used.update(range(int(limits[0]),int(limits[1])+1))
                elif limits[0]:used.add(int(limits[0]))
        assert used<=set(numbers)
        for ref in catalog:
            identifier=ref.get('doi') or ref.get('accession')
            assert identifier and identifier in refs
        result['references']={'numbered':len(numbers),'cited':sorted(used),'uncited':sorted(set(numbers)-used),'scope':'Catalog identifiers and numbered-text consistency only, not renewed identity/full-text verification'}
claims=pd.read_csv(ROOT/'claims.tsv',sep='\t')
assert claims.claim_id.is_unique
for path,group in claims.groupby('source_table'):
    assert group.source_sha256.eq(sha(ROOT/path)).all(),path
result['claims']={'rows':len(claims),'all_source_hashes_current':True}
pagination=json.loads((ROOT/'manuscript/review_pdf/caption_pagination_checks.json').read_text(encoding='utf-8'))
for name,record in result['documents'].items():
    assert pagination['documents'][name]==record['sha256'], 'Word pagination check stale: '+name
assert len(pagination['checks'])==22 and all(x['same_page'] for x in pagination['checks'])
result['word_caption_pagination_current']=True
OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'figures':22,'claims':len(claims),'references':result['references'],'remaining_author_markers':sum(len(v['remaining_markers']) for v in result['documents'].values())},ensure_ascii=False))
