"""Resolve GEO sample contrasts and describe matched model +/-Dox changes."""
from pathlib import Path
import gzip, re, json, hashlib
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/validation/gse144035'
soft=OUT/'GSE144035_family.soft.gz'
raw=OUT/'GSE144035_All_counts_data.txt.gz'
records=[]
for block in gzip.open(soft,'rt').read().split('^SAMPLE = ')[1:]:
    title=re.search(r'!Sample_title = (.+)',block)[1].strip()
    fields=dict(re.findall(r'!Sample_characteristics_ch1 = ([^:]+): ([^\r\n]+)',block))
    records.append({'gsm':block.splitlines()[0].strip(),'title':title,
                    'model':re.match(r'(VTL-\d+)',title)[1],
                    'sample':re.search(r'\[([^]]+)\]',title)[1],
                    'status':fields['lmo2_status'],'dox':fields['dox']})
meta=pd.DataFrame(records)
assert len(meta)==14 and meta['sample'].is_unique
meta.to_csv(OUT/'sample_manifest.tsv',sep='\t',index=False)
d=pd.read_csv(raw,sep='\t')
counts=d.set_index('Gene.ID')[meta['sample']].apply(pd.to_numeric)
assert counts.ge(0).all().all() and np.equal(counts,counts.round()).all().all()
log=np.log2(1+counts.div(counts.sum(),axis=1)*1e6)
rows=[]
for model,m in meta.groupby('model'):
    if set(m.dox)!={'Plus','Minus'}: continue
    ids=m.set_index('dox')['sample']
    for gene in ['Zeb1','Lmo2']:
        geneids=d.loc[d['Gene.Name'].eq(gene),'Gene.ID']
        assert len(geneids)==1
        values=log.loc[geneids.iloc[0]]
        rows.append({'model':model,'status':m.status.iloc[0],'gene':gene,
                     'minus_log2_1pCPM':values[ids['Minus']], 'plus_log2_1pCPM':values[ids['Plus']],
                     'delta_plus_minus':values[ids['Plus']]-values[ids['Minus']]})
pairs=pd.DataFrame(rows)
assert len(pairs)==12
pairs.to_csv(OUT/'paired_model_changes.tsv',sep='\t',index=False)
d.loc[d['Gene.Name'].isin(['Zeb1','Lmo2']),['Gene.ID','Gene.Name','Dependent -Dox','Independent -Dox','FDR','P value']].to_csv(OUT/'author_table_rows.tsv',sep='\t',index=False)
record={'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw,soft]},
        'samples':14,'paired_models':6,'unpaired_dependent_models':2,
        'scale':'log2(1+CPM) from all supplied count rows; descriptive model-matched Plus minus Minus Dox',
        'interpretation':'Independent -Dox author column is not a within-model Plus versus Minus Dox contrast; no independent-biological-replicate or inferential P claim for paired descriptive differences.'}
(OUT/'audit.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(pairs.to_string(index=False))
# Preserve the mislabelled legacy row for audit, then exclude it from perturbation source data.
source=ROOT/'data/SourceData_Fig6E_perturbation_forest.tsv'
effects=pd.read_csv(source,sep='\t')
legacy=effects[effects.dataset.eq('GSE144035')]
if len(legacy):
    legacy.to_csv(OUT/'withdrawn_mislabelled_F6E_row.tsv',sep='\t',index=False)
    effects[~effects.dataset.eq('GSE144035')].to_csv(source,sep='\t',index=False)
