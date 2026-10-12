"""Unchanged fixed normal-reference modules; AIIMS35 FPKM is not count/CPM."""
from pathlib import Path
import json,datetime,hashlib,re,gzip,xml.etree.ElementTree as ET
import pandas as pd,numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
H=Path(__file__).resolve().parent; O=H/'independent_AIIMS35_maturation'; S=H/'public_source_qualification'/'recent_clinical2025_2026'
members_path=H/'multigene_maturity_proxy'/'FIXED200_MODULE_MEMBERS.tsv'; mapping_path=H/'multigene_maturity_proxy'/'GENE_MAPPING_LEDGER.tsv'; normal_path=H.parents[1]/'tasks/ZEB1_TARGET_DISCOVERY/inputs/NORMAL_DONOR_STATE_COUNTS.tsv.gz'
contract_path=O/'EXPLORATORY_FIXED_MODULE_CONTRACT.json'; contract=json.loads(contract_path.read_text()); assert 'FPKM' in contract['score'] and '10000' in contract['family']
source=json.loads((O/'SOURCE_ACQUISITION_AND_IDENTITY_QUALIFICATION.json').read_text()); p=O/'GSE216117_genes.fpkm_table.txt.gz'
assert hashlib.sha256(p.read_bytes()).hexdigest()==source['source_sha256']
assert hashlib.sha256(members_path.read_bytes()).hexdigest()==source['frozen_members_sha256']
assert hashlib.sha256(contract_path.read_bytes()).hexdigest()==source['contract_sha256']
with gzip.open(p,'rt') as f: d=pd.read_csv(f,sep='\t',index_col=0)
d.index=d.index.astype(str).str.replace(r'\.\d+$','',regex=True); assert d.index.is_unique and d.columns.is_unique
assert d.shape==(58219,35) and np.isfinite(d.to_numpy()).all() and (d>=0).all().all()
units=[]
for s in source['samples']:
    title=[x.split(' = ',1)[1] for x in s['fields'] if x.startswith('!Sample_title = ')]; assert len(title)==1
    if title[0] in d.columns:
        units.append({'gsm':s['gsm'],'expression_unit_id':title[0],'origin':'AIIMS_New_Delhi_GSE216117','unit_basis':'35 T-ALL discovery patients / 35 published RNA samples','clinical_id_mapping':'UNESTABLISHED_DI_1_IS_NOT_ASSUMED_D1','source_disease_state':';'.join(x.split(': ',1)[-1] for x in s['fields'] if 'disease state:' in x)})
u=pd.DataFrame(units); assert len(u)==35 and u.gsm.is_unique and u.expression_unit_id.is_unique and set(u.expression_unit_id)==set(d.columns)
u.to_csv(O/'ALL35_DIAGNOSTIC_EXPRESSION_UNIT_LEDGER.tsv',sep='\t',index=False)
root=ET.parse(S/'PMC10847906.xml').getroot(); passages=[]
for sec in root.iter('sec'):
    title=sec.find('title'); title=' '.join(title.itertext()) if title is not None else ''
    if title in ['Patients','Clinical characteristics of patients in the discovery cohort','Data Availability']:
        passages.append({'title':title,'text':' '.join(' '.join(sec.itertext()).split())})
assert any('discovery' in x['text'] and '35' in x['text'] for x in passages)
normal=set(pd.read_csv(normal_path,sep='\t',usecols=[0]).iloc[:,0].astype(str)); mp=pd.read_csv(mapping_path,sep='\t')
mp=mp[mp.eligible&mp.symbol.isin(normal)&mp.stable.isin(d.index)].copy(); assert mp.stable.is_unique and mp.symbol.is_unique
mp.to_csv(O/'EXACT_STABLE_MAPPING_LEDGER.tsv',sep='\t',index=False); e=d.loc[mp.stable].copy(); e.index=mp.symbol.to_numpy()
m=pd.read_csv(members_path,sep='\t'); assert len(m)==200 and m.symbol.nunique()==166 and all(m.groupby(['axis','direction']).size()==50)
cv=[]
for (axis,direction),g in m.groupby(['axis','direction']): cv.append({'axis':axis,'direction':direction,'n_original':50,'n_present':int(g.symbol.isin(e.index).sum()),'missing_symbols':';'.join(g.loc[~g.symbol.isin(e.index),'symbol'])})
cov=pd.DataFrame(cv); passed=len(e)>=20000 and all(cov.n_present>=45) and all(cov.groupby('axis').n_present.sum()>=90) and 'ZEB1' in e.index
ip=H/'public_source_qualification'/'new_independent_clinical_leads'/'PMC10924344.xml'; ir=ET.parse(ip).getroot(); india_aff=[' '.join(' '.join(x.itertext()).split()) for x in ir.iter('aff')]
q={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'QUALIFIED_FOR_FIXED_EXPLORATORY_AIIMS35_MATURATION_ASSOCIATION' if passed else 'COVERAGE_GATE_FAILED_NO_ASSOCIATION_FITTED','n_patients':35,'n_expression_gene_rows':len(d),'n_exact_shared_ranking_symbols':len(e),'module_coverage':cv,'source_representation':'Published hg38 Cufflinks FPKM. No count/CPM conversion.','original_patient_basis':passages,'cross_source_provenance':{'AIIMS35':'AIIMS New Delhi, different source institution from PGIMER India25 and AALL0434; no harmonized global identifiers to conclusively exclude any transfers','India25_affiliations':india_aff,'policy':'No pooling; separate cohorts. No individual clinical joins without exact mapping.'},'source_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [p,members_path,mapping_path,normal_path,contract_path,S/'PMC10847906.xml',ip,Path(__file__)]},'boundary':'External fixed normal-maturation expression association. No linked individual OS times/MRD; not an OS replication, causal perturbation or mature-window clinical interaction.'}
(O/'INPUT_QUALIFICATION.json').write_text(json.dumps(q,indent=2),encoding='utf-8'); print('AIIMS35_INPUT_GATE',q['status'],'ranking',len(e),'coverage',cv,flush=True)
if not passed: raise SystemExit(0)
scores=pd.DataFrame(index=u.expression_unit_id); scores['ZEB1_FPKM']=e.loc['ZEB1']; scores['ZEB1_log2_FPKM_plus1']=np.log2(e.loc['ZEB1']+1)
ranks=e.rank(axis=0,method='average')/len(e); rng=np.random.default_rng(20261010); results=[]
for axis,g in m.groupby('axis'):
    up=g.loc[g.direction.eq('up')&g.symbol.isin(e.index),'symbol']; down=g.loc[g.direction.eq('down')&g.symbol.isin(e.index),'symbol']
    y=ranks.loc[up].mean()-ranks.loc[down].mean(); scores[axis+'_raw']=y; x=scores.ZEB1_log2_FPKM_plus1.to_numpy(); y=y.reindex(scores.index).to_numpy(); rho,pv=stats.spearmanr(x,y); boot=[]
    for _ in range(10000):
        ix=rng.integers(0,35,35); v=stats.spearmanr(x[ix],y[ix]).statistic
        if np.isfinite(v): boot.append(v)
    ci=np.quantile(boot,[.025,.975]); results.append({'axis':axis,'n_patients':35,'spearman_rho':float(rho),'p_two_sided':float(pv),'bootstrap95CI_low':float(ci[0]),'bootstrap95CI_high':float(ci[1]),'n_valid_bootstrap':len(boot),'members_up_used':len(up),'members_down_used':len(down)})
for r,qq in zip(results,multipletests([r['p_two_sided'] for r in results],method='fdr_bh')[1]): r['BH2_q']=float(qq)
scores.rename_axis('expression_unit_id').to_csv(O/'ALL35_FIXED_PATIENT_SCORES.tsv',sep='\t'); pd.DataFrame(results).to_csv(O/'ALL2_ASSOCIATION_RESULTS.tsv',sep='\t',index=False)
axis_columns=[a+'_raw' for a in m.axis.unique()]; inter=stats.spearmanr(scores[axis_columns[0]],scores[axis_columns[1]])
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'EXPLORATORY_AIIMS35_FIXED_MATURATION_ASSOCIATIONS_ESTIMATED','results':results,'n_ranking_symbols':len(e),'maturation_axes_correlation':{'axes':axis_columns,'spearman_rho':float(inter.statistic),'p':float(inter.pvalue)},'input_qualification_sha256':hashlib.sha256((O/'INPUT_QUALIFICATION.json').read_bytes()).hexdigest(),'boundary':q['boundary'],'cross_modality_limit':'Gene-length normalized FPKM ranks differ from NOPHO count ranks and main-cohort TPM ranks; no pooled effects or absolute maturity thresholds.'}
(O/'RESULTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8'); print('AIIMS35_FIXED_RESULTS',json.dumps(out),flush=True)
