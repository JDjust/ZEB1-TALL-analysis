"""Verify existing perturbation sample metadata; no expression reanalysis."""
from pathlib import Path
import csv, gzip, hashlib, json, shutil
import pandas as pd

root = Path(__file__).resolve().parents[1]
out = root / 'data/validation/perturbation_units'
out.mkdir(parents=True, exist_ok=True)
legacy = Path('D:/ZEB1/_new ZEB1')
soft = legacy/'_geo_dl/GSE188225_family.soft.gz'
source = legacy/'LMO2_ZEB1_module_closeout/GSE188225_sample_manifest_corrected.csv'
records = {}
for block in gzip.open(soft, 'rt', encoding='utf-8').read().split('^SAMPLE = ')[1:]:
    gsm = block.splitlines()[0].strip()
    rec = {}
    for line in block.splitlines():
        if line.startswith('!Sample_characteristics_ch1 = '):
            key, value = line.split(' = ', 1)[1].split(': ', 1)
            rec[key] = value
        elif line.startswith('!Sample_description = '):
            rec['count_column'] = line.split(' = ', 1)[1]
    records[gsm] = rec
df = pd.read_csv(source, dtype=str)
for row in df.to_dict('records'):
    rec = records[row['geo_accession']]
    for field, geo in [('count_column','count_column'),('mouse_id','mouse'),('age_weeks','age (weeks)'),('age_category','age category'),('genotype_geo','genotype'),('sex','gender')]:
        assert row[field] == rec[geo], (row['geo_accession'], field)
df.to_csv(out/'GSE188225_verified_samples.tsv', sep='\t', index=False)
summary = []
for (age, stage, genotype), group in df.groupby(['age_category','stage','genotype'], sort=False):
    summary.append(dict(dataset='GSE188225', age=age, stage=stage, genotype=genotype,
                        libraries=len(group), unique_mice=group.mouse_id.nunique(),
                        mouse_ids=';'.join(group.mouse_id), gsm_ids=';'.join(group.geo_accession)))
pd.DataFrame(summary).to_csv(out/'GSE188225_group_units.tsv', sep='\t', index=False)
overlap = df.groupby('mouse_id').filter(lambda x: len(x)>1)
overlap.to_csv(out/'GSE188225_shared_mice.tsv', sep='\t', index=False)
assert len(df)==17 and df.mouse_id.nunique()==12 and overlap.mouse_id.nunique()==5
shutil.copy2(soft, out/soft.name)
shutil.copy2(source, out/source.name)
matrix = legacy/'_geo_dl/GSE49164_series_matrix.txt.gz'
meta = {}
with gzip.open(matrix, 'rt', encoding='utf-8') as stream:
    for line in stream:
        if line.startswith(('!Sample_title\t','!Sample_geo_accession\t')):
            cells = next(csv.reader([line], delimiter='\t'))
            meta[cells[0]] = cells[1:]
array_samples = pd.DataFrame({'gsm':meta['!Sample_geo_accession'],'title':meta['!Sample_title']})
array_samples['group'] = array_samples.title.map(lambda x: 'TG_Lyl1_KO' if 'knockout' in x else ('WT' if x.startswith('Wild-type') else 'TG'))
assert array_samples.group.value_counts().to_dict() == {'WT':6,'TG':3,'TG_Lyl1_KO':3}
array_samples.to_csv(out/'GSE49164_deposited_samples.tsv', sep='\t', index=False)
shutil.copy2(matrix, out/matrix.name)
inventory = [
    ['GSE110635','Human ALL-SIL; TLX1 siRNA, 24 h','3 controls; 3 s46; 3 s47 libraries','Two siRNA contrasts share three controls; replicate index retained, not independent cohorts','../regulon_coverage/gse110635_symbol_coldata.tsv'],
    ['GSE188225','Mouse Lmo2-TG vs WT, young DN2 / young DN3a / old DN3a','17 libraries from 12 mouse IDs; 5 mice contribute both young stages','Young DN2: TG 3 / WT 2; young DN3a: 3 / 3; old DN3a: 3 / 3. Stages not independent studies','GSE188225_verified_samples.tsv'],
    ['GSE49164','Mouse DN thymocytes, 2 months; Lmo2-TG vs WT','6 WT / 3 TG arrays; separate TG+Lyl1-KO arm has 3 arrays','GEO replicate-labelled arrays; individual mouse IDs unresolved. Main Figure 6 uses TG vs WT, not the Lyl1-KO arm','GSE49164_deposited_samples.tsv'],
    ['GSE287751 bulk','Mouse DN1 Lmo2 KO vs control','2 libraries per arm in existing source table','Descriptive TPM-derived change; animal identities not resolved here','../../../../modules/module10/tables/M10_r2_GSE287751_DN1_KO_vs_ctrl.tsv'],
    ['GSE287751 single-cell','Mouse Lmo2 KO vs EV; day 3 / preNotch; unprimed / flt3l','8 author-labelled conditions; 11661 retained cells','Four context contrasts are not four biological replicates; cells are not independent experimental units','../gse287751/analysis_scope.json'],
    ['GSE144035','Mouse VTL model-matched -Dox / +Dox','6 paired models among 14 deposited libraries','Three evolving and three independent models; recipient-animal independence not established','../gse144035/sample_manifest.tsv'],
]
with (out/'perturbation_unit_inventory.tsv').open('w', newline='', encoding='utf-8') as f:
    writer=csv.writer(f, delimiter='\t'); writer.writerow(['dataset','context','observed_units','independence_boundary','source']); writer.writerows(inventory)
audit = {'verified_fields':['GSM','count_column','mouse','age','age category','genotype','sex'],
         'libraries':17,'unique_mice':12,'shared_young_stage_mice':5,
         'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [soft,source,matrix]},
         'scope':'Metadata verification only; no model refitting; unresolved animal identities explicitly retained'}
(out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))
