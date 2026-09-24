"""ID-level reuse inventory from frozen exports, without statistical refitting."""
from pathlib import Path
import hashlib
import itertools
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/validation/cohort_overlap'
OUT.mkdir(exist_ok=True)
inputs=[]
def read(name, sep='\t'):
    p=ROOT/'data'/name
    inputs.append({'path':str(p.relative_to(ROOT)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    return pd.read_csv(p,sep=sep)
groups={}; namespaces={}; refs={}
def add(name, values, namespace, source):
    values=list(values)
    assert not pd.isna(values).any()
    assert len(values)==len(set(values)), f'Duplicate IDs in {name}'
    groups[name]=set(map(str,values)); namespaces[name]=namespace; refs[name]=source

p='validation/lmo2_subgroups_existing/TALL_clinical_annotation_full.csv'
d=read(p,',')
add('TARGET_all',d.Sample_ID,'TARGET patient',p)
add('TARGET_classified',d.loc[d.MolecularSubtype.ne('Unknown'),'Sample_ID'],'TARGET patient',p)
add('TARGET_LMO_associated',d.loc[d.MolecularSubtype.isin(['LMO1/2','LMO2_LYL1']),'Sample_ID'],'TARGET patient',p)
p='validation/lmo2_within_subtype/sample_inputs.tsv'; d=read(p)
q=read('validation/lmo2_within_subtype/within_subtype.tsv')
labels=q.loc[q.cohort.eq('TARGET-ALL-P2') & q.inference,'subtype']
add('TARGET_interaction',d.loc[d.cohort.eq('TARGET-ALL-P2') & d.subtype.isin(labels),'sample'],'TARGET patient',p)
p='validation/clinical_endpoint_audit/patient_protocol_mrd_manifest.tsv'; d=read(p)
add('Pharmacotype_all',d['sample'],'Pharmacotype patient',p)
add('Pharmacotype_paired_MRD',d.loc[d.paired_mrd,'sample'],'Pharmacotype patient',p)
add('Pharmacotype_persistence',d.loc[d.persistence_model_eligible,'sample'],'Pharmacotype patient',p)
p='SourceData_Fig7F_mtx_points.tsv'; d=read(p)
d=d[d.lineage.eq('T-ALL')].dropna(subset=['ZEB1','value','model_id'])
for source in ['GDSC1','GDSC2']:
    add(source+'_MTX',d.loc[d.source.eq(source),'model_id'],'DepMap ModelID',p)
p='SourceData_Fig7C_ZEB1_dependency.tsv'; d=read(p).dropna(subset=['ZEB1_GE'])
add('DepMap_ZEB1_CRISPR',d.ModelID,'DepMap ModelID',p)
p='validation/perturbation_units/GSE188225_verified_samples.tsv'; d=read(p)
for stage in ['DN2','DN3a']:
    part=d[d.stage.eq(stage) & d.age_category.eq('Young')]
    add('GSE188225_young_'+stage,part.mouse_id,'GSE188225 mouse',p)

membership=[{'analysis_set':name,'id_namespace':namespaces[name],'id':id,'source_table':refs[name]}
            for name,ids in groups.items() for id in sorted(ids)]
pd.DataFrame(membership).to_csv(OUT/'memberships.tsv',sep='\t',index=False)
rows=[]
for left,right in itertools.combinations(groups,2):
    if namespaces[left]!=namespaces[right]:
        continue # different ID namespaces do not establish independent subjects
    a,b=groups[left],groups[right]
    rows.append({'set_a':left,'set_b':right,'id_namespace':namespaces[left],
                 'n_a':len(a),'n_b':len(b),'shared_n':len(a&b),'only_a_n':len(a-b),
                 'only_b_n':len(b-a),'relationship':'identical' if a==b else 'nested' if a<=b or b<=a else 'partly_overlapping' if a&b else 'no_matching_ID_in_this_namespace',
                 'shared_ids':';'.join(sorted(a&b))})
pd.DataFrame(rows).to_csv(OUT/'overlap_edges.tsv',sep='\t',index=False)
(OUT/'provenance.json').write_text(json.dumps({'inputs':inputs,'scope':'Exact IDs in specified local exports only; unmatched namespaces are not evidence of distinct patients. Shared cohorts may have subset-specific measurements; no statistics rerun.'},indent=2),encoding='utf-8')
print(pd.DataFrame(rows).drop(columns=['shared_ids']).to_string(index=False))
