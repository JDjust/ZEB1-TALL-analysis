"""Bounded source-record integrity check; absence of flags is not a guarantee."""
from pathlib import Path
import json, hashlib, re, csv, html
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import requests

root = Path(__file__).resolve().parents[1]
out = root/'data/reference_records/integrity'
out.mkdir(exist_ok=True)
catalog = json.loads((root/'data/reference_catalog.json').read_text(encoding='utf-8'))['records']
cached = {}
for path in (root/'data/reference_records').glob('pubmed*.xml'):
    for a in ET.fromstring(path.read_bytes()).findall('PubmedArticle'):
        doi = next((x.text for x in a.findall('.//ArticleId') if x.get('IdType')=='doi'), '')
        if doi: cached[doi.lower()] = (a, str(path.relative_to(root)))
selected = [1,2,3,4,6,28]
now = datetime.now(timezone.utc).isoformat()

def crossref(number):
    ref = catalog[number-1]
    url = 'https://api.crossref.org/works/'+ref['doi']
    path = out/f'crossref_ref{number}.json'
    try:
        if path.exists():
            raw=path.read_bytes()
        else:
            response = requests.get(url, timeout=(8,25))
            response.raise_for_status()
            raw=response.content; path.write_bytes(raw)
        m=json.loads(raw)['message']
        assert m['DOI'].lower()==ref['doi'].lower()
        normalized=lambda s: re.sub(r'\W','',html.unescape(re.sub(r'<[^>]+>','',s)).lower())
        assert normalized(m['title'][0])==normalized(ref['title'])
        return {'reference':number,'doi':ref['doi'],'checked_utc':now,'url':url,
                'status':'retrieved; identity matched','record_saved_utc':datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat(),'update_to':m.get('update-to',[]),
                'relation':m.get('relation',{}),'source':str(path.relative_to(root)),
                'sha256':hashlib.sha256(raw).hexdigest()}
    except Exception as e:
        return {'reference':number,'doi':ref['doi'],'checked_utc':now,'url':url,'status':'failed','error':str(e)}

with ThreadPoolExecutor(max_workers=3) as pool:
    live=list(pool.map(crossref,selected))
pmids=[]
for i in selected:
    doi=catalog[i-1]['doi'].lower()
    if doi in cached: pmids.append(cached[doi][0].findtext('MedlineCitation/PMID'))
if '29778661' not in pmids: pmids.append('29778661')
for i in [2,3,6]:
    searchpath=out/f'pubmed_search_ref{i}.json'
    if searchpath.exists():
        result=json.loads(searchpath.read_bytes())
    else:
        response=requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'pubmed','term':catalog[i-1]['doi']+'[AID]','retmode':'json'},timeout=(8,25))
        response.raise_for_status(); result=response.json(); searchpath.write_bytes(response.content)
    ids=result['esearchresult']['idlist']
    if len(ids)==1 and ids[0] not in pmids: pmids.append(ids[0])
url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
fetch={'url':url,'pmids':pmids,'checked_utc':now}
try:
    r=requests.get(url,params={'db':'pubmed','id':','.join(pmids),'retmode':'xml'},timeout=(8,30))
    r.raise_for_status(); tree=ET.fromstring(r.content)
    assert tree.findall('PubmedArticle')
    path=out/'live_core_pubmed.xml'; path.write_bytes(r.content)
    for a in tree.findall('PubmedArticle'):
        doi=next((x.text for x in a.findall('.//ArticleId') if x.get('IdType')=='doi'),'')
        if doi: cached[doi.lower()]=(a,str(path.relative_to(root)))
    fetch['status']='retrieved'; fetch['sha256']=hashlib.sha256(r.content).hexdigest()
except Exception as e: fetch.update(status='failed',error=str(e))

rows=[]
for i,ref in enumerate(catalog,1):
    doi=ref.get('doi','').lower()
    if not doi: continue
    row={'reference':i,'doi':doi,'pmid':'','record_source':'','publication_types':'','linked_notices':'','scope':'PubMed record unavailable in inspected local set'}
    if doi in cached:
        a,source=cached[doi]
        row.update(pmid=a.findtext('MedlineCitation/PMID'),record_source=source,
                   publication_types='; '.join(''.join(x.itertext()) for x in a.findall('.//PublicationType')),
                   linked_notices=json.dumps([{'type':x.get('RefType'),'pmid':x.findtext('PMID'),'source':x.findtext('RefSource'),'note':x.findtext('Note')} for x in a.findall('.//CommentsCorrections')],ensure_ascii=False),
                   scope='Live retrieved record' if 'live_core_pubmed' in source else 'Previously cached record; not refreshed this pass')
    rows.append(row)
with (out/'pubmed_notice_inventory.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
(out/'run.json').write_text(json.dumps({'checked_utc':now,'crossref':live,'pubmed_fetch':fetch},indent=2),encoding='utf-8')
print(json.dumps({'crossref':[(x['reference'],x['status'],x.get('update_to')) for x in live],'pubmed':fetch,'records_with_notices':[x for x in rows if x['linked_notices'] not in ('','[]')]},ensure_ascii=False,indent=2))
