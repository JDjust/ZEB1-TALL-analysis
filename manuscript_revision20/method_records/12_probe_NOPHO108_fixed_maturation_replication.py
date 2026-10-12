"""Fixed-module independent expression bridge; never reconstruct missing OS/MRD."""
from pathlib import Path
import json,datetime,hashlib,urllib.request,gzip,io,re
import pandas as pd,numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
H=Path(__file__).resolve().parent;S=H/'public_source_qualification'/'recent_clinical2025_2026';O=H/'independent_NOPHO108_maturation';O.mkdir(exist_ok=True)
members_path=H/'multigene_maturity_proxy'/'FIXED200_MODULE_MEMBERS.tsv';mapping_path=H/'multigene_maturity_proxy'/'GENE_MAPPING_LEDGER.tsv';normal_path=H.parents[1]/'tasks/ZEB1_TARGET_DISCOVERY/inputs/NORMAL_DONOR_STATE_COUNTS.tsv.gz'
contract={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'role':'Independent normal-reference fixed-module maturation expression association. Not independent OS or response validation.','prior_exposure':'Original Blood2025 and clinTALL2026 papers, source headers and CIMP aggregate outcomes read; no ZEB1/maturity association coefficients viewed. Exploratory external probe.','patient_source':'108 diagnostic NOPHO ALL2008 T-ALL specimens in GSE272023; exclude3 normal reference samples. Original paper describes one sample per108 patient. No COG AALL0434 overlap inferred from origin; direct sample IDs compared to main1309 and held India/AIEOP cohorts.','members':'Exactly same200 normal donor-selected members; no reselection/substitutions. Exact eligible stable ENSG map.','gate':'Same90/100 per axis and45/50 per direction; ranking universe at least20000 shared original normal symbols. If fail, no coefficients.','score':'Within-sample average count ranks divided by number of common genes; mean frozen up minus mean frozen down. Ranks library-scale invariant but ties, gene length, transcript representation and sparsity remain limitations.','ZEB1':'log2(CPM+1), library denominator all published RNA gene counts, before common-symbol restriction. Not outcome-selected probe or dichotomized expression.','family':'Two two-sided Spearman associations (ZEB1 vs each module), BH2;10000 paired-patient bootstrap95%CI, seed20261010. Keep both directions. No CIMP subset selection, outcome fit or tuning.','interpretation':'Association with fixed normal-state similarity is not full developmental trajectory, human perturbation mechanism, window-specific clinical effect or OS replication.'}
(O/'EXPLORATORY_FIXED_MODULE_CONTRACT.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
headers=json.loads((S/'GSE272023_ALL_SAMPLE_HEADERS.json').read_text());led=[]
for s in headers['samples']:
    desc=[l.split(' = ',1)[1] for l in s['fields'] if l.startswith('!Sample_description')];pid=next((x for x in desc if re.fullmatch(r'TALL\d+',x)),None)
    if pid:led.append({'gsm':s['gsm'],'patient_id':pid,'origin':'NOPHO_ALL2008','diagnostic':any('taken at diagnosis' in l for l in s['fields'])})
units=pd.DataFrame(led);assert len(units)==108 and units.gsm.is_unique and units.patient_id.is_unique and units.diagnostic.all();units.to_csv(O/'ALL108_PATIENT_UNIT_LEDGER.tsv',sep='\t',index=False)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272023/suppl/GSE272023_MatrixRNAseqRawCounts.txt.gz';p=O/'GSE272023_MatrixRNAseqRawCounts.txt.gz'
if not p.exists():
    with urllib.request.urlopen(url,timeout=35) as f:data=f.read()
    gzip.decompress(data);p.write_bytes(data)
with gzip.open(p,'rt') as f:d=pd.read_csv(f,sep='\t',index_col=0)
print('NOPHO108_ORIGINAL_COUNT_MATRIX',d.shape,'genes',d.index[:4].tolist(),'columns',d.columns[:8].tolist(),flush=True)
d.index=d.index.astype(str).str.replace(r'\.\d+$','',regex=True);assert d.index.is_unique and d.columns.is_unique
assert set(units.patient_id)<=set(d.columns);patients=d.loc[:,units.patient_id].copy();assert np.isfinite(patients.to_numpy()).all() and (patients>=0).all().all() and np.allclose(patients.to_numpy(),np.round(patients.to_numpy()),rtol=0,atol=1e-8)
assert len(set(d.columns)-set(units.patient_id))==3
normal=set(pd.read_csv(normal_path,sep='\t',usecols=[0]).iloc[:,0].astype(str));mp=pd.read_csv(mapping_path,sep='\t');mp=mp[mp.eligible&mp.symbol.isin(normal)&mp.stable.isin(d.index)].copy();assert mp.stable.is_unique and mp.symbol.is_unique
mp.to_csv(O/'EXACT_STABLE_MAPPING_LEDGER.tsv',sep='\t',index=False);e=patients.loc[mp.stable].copy();e.index=mp.symbol.to_numpy();m=pd.read_csv(members_path,sep='\t');assert len(m)==200 and all(m.groupby(['axis','direction']).size()==50)
cv=[]
for (axis,direction),g in m.groupby(['axis','direction']):cv.append({'axis':axis,'direction':direction,'n_original':50,'n_present':int(g.symbol.isin(e.index).sum()),'missing_symbols':';'.join(g.loc[~g.symbol.isin(e.index),'symbol'])})
cov=pd.DataFrame(cv);passed=len(e)>=20000 and all(cov.n_present>=45) and all(cov.groupby('axis').n_present.sum()>=90) and 'ZEB1' in e.index
q={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'QUALIFIED_FOR_FIXED_EXPLORATORY_INDEPENDENT_MATURATION_ASSOCIATION' if passed else 'COVERAGE_GATE_FAILED_NO_ASSOCIATION_FITTED','n_patients':108,'n_normal_references_excluded':3,'n_expression_gene_rows':len(d),'n_exact_shared_ranking_symbols':len(e),'source_url':url,'module_coverage':cv,'source_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [p,members_path,mapping_path,normal_path,S/'GSE272023_ALL_SAMPLE_HEADERS.json',O/'EXPLORATORY_FIXED_MODULE_CONTRACT.json',Path(__file__)]},'boundary':contract['interpretation']}
(O/'INPUT_QUALIFICATION.json').write_text(json.dumps(q,indent=2),encoding='utf-8');print('NOPHO108_INPUT_GATE',q['status'],'n_ranking',len(e),'coverage',cv,flush=True)
if not passed:raise SystemExit(0)
lib=patients.sum(axis=0);assert (lib>0).all();scores=pd.DataFrame(index=units.patient_id);scores['library_total']=lib;scores['ZEB1_raw_count']=e.loc['ZEB1'];scores['ZEB1_log2_CPM_plus1']=np.log2(e.loc['ZEB1']/lib*1e6+1)
ranks=e.rank(axis=0,method='average')/len(e);rng=np.random.default_rng(20261010);results=[]
for axis,g in m.groupby('axis'):
    up=g.loc[g.direction.eq('up')&g.symbol.isin(e.index),'symbol'];down=g.loc[g.direction.eq('down')&g.symbol.isin(e.index),'symbol'];y=ranks.loc[up].mean()-ranks.loc[down].mean();scores[axis+'_raw']=y;x=scores.ZEB1_log2_CPM_plus1.to_numpy();y=y.reindex(scores.index).to_numpy();rho,pv=stats.spearmanr(x,y);boot=[]
    for _ in range(10000):
        ix=rng.integers(0,108,108);v=stats.spearmanr(x[ix],y[ix]).statistic
        if np.isfinite(v):boot.append(v)
    ci=np.quantile(boot,[.025,.975]);results.append({'axis':axis,'n_patients':108,'spearman_rho':float(rho),'p_two_sided':float(pv),'bootstrap95CI_low':float(ci[0]),'bootstrap95CI_high':float(ci[1]),'n_valid_bootstrap':len(boot),'members_up_used':len(up),'members_down_used':len(down)})
for r,qq in zip(results,multipletests([r['p_two_sided'] for r in results],method='fdr_bh')[1]):r['BH2_q']=float(qq)
scores.rename_axis('patient_id').to_csv(O/'ALL108_FIXED_PATIENT_SCORES.tsv',sep='\t');pd.DataFrame(results).to_csv(O/'ALL2_ASSOCIATION_RESULTS.tsv',sep='\t',index=False)
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'EXPLORATORY_INDEPENDENT_NOPHO108_FIXED_MATURATION_ASSOCIATIONS_ESTIMATED','results':results,'n_ranking_symbols':len(e),'input_qualification_sha256':hashlib.sha256((O/'INPUT_QUALIFICATION.json').read_bytes()).hexdigest(),'boundary':contract['interpretation']}
(O/'RESULTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print('NOPHO108_FIXED_MATURATION_ASSOCIATIONS',json.dumps(results),flush=True)
