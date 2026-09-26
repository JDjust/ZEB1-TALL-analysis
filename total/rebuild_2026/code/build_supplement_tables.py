"""Assemble manuscript supplement tables strictly from frozen results."""
from pathlib import Path
import csv, shutil

base = Path(__file__).resolve().parents[1]
valid = base.parent / 'data' / 'validation'
out = base / 'data' / 'source_data_rebuilt' / 'supplementary_tables'
out.mkdir(parents=True, exist_ok=True)

def read(rel):
    with (valid / rel).open(encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write(name, rows, fields):
    with (out / name).open('w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=fields)
        w.writeheader(); w.writerows(rows)

models = [
    ('Induction failure','OR','polonen_round1/model_induction_failure.tsv',[0,1]),
    ('M2/M3 morphology','OR','polonen_round1/model_poor_morphology.tsv',[0,1]),
    ('MRD >=0.1%','OR','polonen_round1/model_mrd.tsv',[0,1]),
    ('MRD >=0.01%','OR','polonen_round1/model_mrd.tsv',[2,3]),
    ('Event-free survival','HR','polonen_round1/model_survival.tsv',[0,1]),
    ('Overall survival','HR','polonen_round1/model_survival.tsv',[2,3]),
]
clin=[]
for ep,measure,rel,indices in models:
    rows=read(rel)
    for adjustment,idx in zip(('Unadjusted','Subtype-adjusted'),indices):
        row=rows[idx]
        clin.append(dict(endpoint=ep,measure=measure,adjustment=adjustment,
            method=row.get('method','Cox model'),n=row['n'],events=row['events'],
            estimate=row.get('or',row.get('hr')),
            ci_low=row['ci_low'],ci_high=row['ci_high'],p=row['p'],
            source=rel,model=row['model']))
def bh(values):
    order=sorted(range(len(values)), key=lambda i: values[i])
    result=[1.0]*len(values)
    carry=1.0
    for rank,i in reversed(list(enumerate(order, start=1))):
        carry=min(carry, values[i]*len(values)/rank)
        result[i]=carry
    return result
for adjustment in ('Unadjusted','Subtype-adjusted'):
    positions=[i for i,r in enumerate(clin) if r['adjustment']==adjustment]
    q=bh([float(clin[i]['p']) for i in positions])
    for i,value in zip(positions,q):
        clin[i]['endpoint_bh_fdr']=f'{value:.12g}'
write('Table_S1_clinical_models.tsv',clin,list(clin[0]))

datasets=[
 ('GSE142522','normal thymus','sorted bulk stage libraries','6 stage summaries','normal stage shape; CD8 SP recovery'),
 ('GSE195812','normal thymus','FACS stage libraries pooled from six donors','8 stage summaries','normal stage shape and existing cell embedding'),
 ('GSE206710','normal thymus','three-donor scRNA stage aggregate','3 stage summaries','DN-to-DP shape; reference units'),
 ('Park/HTA','normal thymus','donor-stage aggregates','4 broad stages; 15 eligible paired donors','donor-level DN to CD4 SP support'),
 ('Pölönen 2024','diagnostic T-ALL','one patient/one RNA sample','1,309 patients; 17 molecular subtypes','primary balance, developmental residual and clinical models'),
 ('GSE146901','ETP/non-ETP T-ALL and healthy peripheral T','RNA sample','8 ETP + 10 non-ETP + 4 peripheral T','independent expression polarity; pooled loops negative'),
 ('Lim 2024','T-ALL malignant single-cell states','patient-level pseudobulk','54 Day0; 41 eligible ZBTB16 state pairs','within-patient state localization; Day28 sensitivity'),
 ('GSE162280','BCL11B-rearranged immature leukemia','case RNA library','12 rearranged cases','within-lesion partner convergence only'),
 ('HiChIP collection','T-ALL, CD34 and thymus','filtered contact library','4 non-ETP + 3 ETP + 2 CD34 + 1 usable thymus','orthogonal contact support'),
 ('scATAC collection','primary/PDX T-ALL','peak set','7 peak sets','negative independent chromatin comparison'),
]
write('Table_S2_dataset_units.tsv',
      [dict(zip(('dataset','modality','biological_unit','count','role'),r)) for r in datasets],
      ['dataset','modality','biological_unit','count','role'])

excluded=[
 ('MTX response','historical drug-response screens','outside the frozen developmental/subtype question; no independent biomarker claim retained'),
 ('DDR/NHEJ','exploratory pathway summaries','not used to infer a ZEB mechanism from the available validation evidence'),
 ('NOTCH dependency','perturbation/dependency exploration','no subtype-specific dependency conclusion carried into this paper'),
 ('Dormancy','exploratory state interpretation','not established by the frozen patient-level results'),
 ('GPL19197 TLX1 KD/OE','legacy perturbation-array interpretation','not accepted as a valid independent perturbation validation'),
 ('HOXA historical contrast','older bulk cohort exploration','not among the robust scoped TLX/TAL signals'),
]
write('Table_S3_deprioritized_lines.tsv',
      [dict(zip(('line','prior_material','reason_not_in_main_story'),r)) for r in excluded],
      ['line','prior_material','reason_not_in_main_story'])

old = {r['subtype']: r for r in read('zeb_developmental_residual/subtype_residual_summary.tsv')}
with (base/'data/source_data_rebuilt/revision2_statistics/subtype_effects_bootstrap.tsv').open(
        encoding='utf-8-sig', newline='') as fh:
    effects = list(csv.DictReader(fh, delimiter='\t'))
def subtype_key(s):
    return s.replace('<U+03B3><U+03B4>', 'γδ').replace('<U+03B1><U+03B2>', 'αβ')
old = {subtype_key(k): v for k, v in old.items()}
full=[]
for e in effects:
    o=old[subtype_key(e['subtype'])]
    full.append({**e, 'dev_median':o['dev_median'],
                 'balance_median':o['balance_median'],
                 'normal_reference_median_sensitivity':o['residual_normal_median'],
                 'n_outside_normal_dev_range':o['n_outside_normal_dev_range']})
write('Table_S4_full_subtype_residuals.tsv', full, list(full[0]))
print('Supplementary tables S1-S4 written')
