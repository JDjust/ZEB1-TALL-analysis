"""Audit exact source endpoint categories; never turn unknown events into censoring."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/validation/target_survival'
source = OUT / 'target_clinical_combined.tsv'
clin = pd.read_csv(source, sep='\t', keep_default_na=False)
clin['_k'] = clin['TARGET USI'].str.replace('-', '', regex=False).str.extract(r'([A-Z0-9]{5,6})$')[0]
old = pd.read_csv(ROOT / 'data/F7_target_surv_merged.tsv', sep='\t')
m = old.merge(clin, on='_k', how='left', validate='one_to_one', indicator=True)
assert m['_merge'].eq('both').all()
categories = {'None': 0, 'Censored': 0, 'Relapse': 1, 'Death': 1,
              'Progression': 1, 'Second Malignant Neoplasm': 1}
assert set(m['First Event']) <= set(categories) | {''}
assert set(m['Vital Status']) <= {'Alive', 'Dead', ''}
m['audited_efs_event'] = m['First Event'].map(categories)
m['audited_os_event'] = m['Vital Status'].map({'Alive': 0, 'Dead': 1})
mapping = {'os_time': 'Overall Survival Time in Days', 'efs_time': 'Event Free Survival Time in Days',
           'age_days': 'Age at Diagnosis in Days'}
checks = {}
for target, raw in mapping.items():
    checks[target] = bool(np.allclose(m[target], pd.to_numeric(m[raw], errors='coerce'), equal_nan=True))
for endpoint in ['os', 'efs']:
    checks[endpoint + '_event'] = bool(np.allclose(m[endpoint+'_event'], m['audited_'+endpoint+'_event'], equal_nan=True))
assert all(checks.values()), checks
workbook = OUT/'TARGET_ALL_ClinicalData_Phase_II_Validation_20230727.xlsx'
raw = pd.read_excel(workbook, keep_default_na=False)
raw.columns = raw.columns.str.replace(r'\s+', ' ', regex=True).str.strip()
joined = m.merge(raw, on='TARGET USI', how='left', validate='one_to_one', suffixes=('', '_xlsx'), indicator='workbook_match')
assert joined.workbook_match.eq('both').all()
for col in ['First Event','Vital Status']:
    checks[col+'_xlsx'] = bool(joined[col].eq(joined[col+'_xlsx']).all())
for col in mapping.values():
    checks[col+'_xlsx'] = bool(np.allclose(pd.to_numeric(joined[col],errors='coerce'),pd.to_numeric(joined[col+'_xlsx'],errors='coerce'),equal_nan=True))
assert all(checks.values()), checks
cols = ['_k','TARGET USI','source_file','First Event','Vital Status', 'ZEB1', 'age_days',
        'os_time','os_event','efs_time','efs_event','audited_os_event','audited_efs_event']
for endpoint in ['os','efs']:
    m[endpoint+'_included'] = m[[endpoint+'_time',endpoint+'_event','ZEB1','age_days']].notna().all(axis=1) & m[endpoint+'_time'].gt(0)
    cols += [endpoint+'_included']
m[cols].to_csv(OUT/'patient_endpoint_audit.tsv',sep='\t',index=False)
m.groupby(['First Event','Vital Status'],dropna=False).size().rename('n').reset_index().to_csv(OUT/'event_categories.tsv',sep='\t',index=False)
refit = pd.read_csv(OUT/'cox_refit.tsv',sep='\t')
prior = pd.read_csv(ROOT/'data/F7_target_cox.tsv',sep='\t')
comparison = refit[refit.ties.eq('breslow')].merge(prior,on=['endpoint','model','term'],suffixes=('_refit','_original'),validate='one_to_one')
for col in ['hr','ci_lo','ci_hi','p']:
    comparison[col+'_abs_difference'] = abs(comparison[col+'_refit']-comparison[col+'_original'])
comparison.to_csv(OUT/'original_vs_survival_package.tsv',sep='\t',index=False)
record = {'scope':'source exported clinical table to patient endpoint and R survival diagnostics',
          'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
          'endpoint_checks':checks,'matched_patients':len(m),
          'literal_None_no_event':int(m['First Event'].eq('None').sum()),
          'truly_blank_First_Event':int(m['First Event'].eq('').sum()),
          'max_HR_difference_breslow':float(comparison.hr_abs_difference.max()),
          'raw_clinical_workbook_rechecked':True,
          'workbook_sha256':hashlib.sha256(workbook.read_bytes()).hexdigest()}
(OUT/'audit_summary.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2))
