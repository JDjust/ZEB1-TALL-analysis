"""Retrieve source metadata for the 16 already-plotted tracks; no signal processing."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,re
import requests
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/validation/chromatin_windows/source_metadata'
OUT.mkdir(exist_ok=True)
manifest=pd.read_csv(OUT.parent/'bigwig_manifest.tsv',sep='\t')
def fetch(name):
    gsm=name.split('_')[0]
    url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
    path=OUT/(gsm+'.soft.txt')
    if not path.exists():
        response=requests.get(url,timeout=25);response.raise_for_status()
        if '!Sample_geo_accession = '+gsm not in response.text:
            raise ValueError('Not a GEO sample record: '+gsm)
        path.write_text(response.text,encoding='utf-8')
    text=path.read_text(encoding='utf-8')
    # GEO's text endpoint hard-wraps values, including URLs, across physical lines.
    logical=[]
    for line in text.splitlines():
        if line.startswith(('!','^','#')) or not logical:
            logical.append(line)
        elif line.strip():
            logical[-1] += line
    text='\n'.join(logical)
    def vals(field):return re.findall(r'^!Sample_'+field+r' = (.*)$',text,re.M)
    processing=vals('data_processing')
    return {'file':name,'gsm':gsm,'source_url':url,'title':' | '.join(vals('title')),
        'characteristics':' | '.join(vals(r'characteristics_ch\d+')),
        'organism':' | '.join(vals(r'organism_ch\d+')),
        'processing':' | '.join(processing),
        'supplementary_files':' | '.join(vals(r'supplementary_file(?:_\d+)?')),
        'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
with ThreadPoolExecutor(max_workers=4) as pool:
    rows=list(pool.map(fetch,manifest.file))
pd.DataFrame(rows).to_csv(OUT/'track_processing.tsv',sep='\t',index=False)
for r in rows:
    print(r['gsm'],r['title'],r['processing'])
