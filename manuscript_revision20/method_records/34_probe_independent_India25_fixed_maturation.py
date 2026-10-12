"""Exploratory independent expression/state probe with frozen normal modules.
No OS reconstruction; no clinical response model; no patient-driven gene selection.
"""
from pathlib import Path
import json,datetime,hashlib,re,zipfile
import numpy as np,pandas as pd,openpyxl
from scipy import stats
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];S=HERE/'public_source_qualification'/'new_independent_clinical_leads';OUT=HERE/'independent_India25_maturation';OUT.mkdir(exist_ok=True)
xp=S/'PMC10924344__12885_2024_12063_MOESM1_ESM.xlsx';members_path=HERE/'multigene_maturity_proxy'/'FIXED200_MODULE_MEMBERS.tsv';mapping_path=HERE/'multigene_maturity_proxy'/'GENE_MAPPING_LEDGER.tsv';normal_path=ROOT/'tasks/ZEB1_TARGET_DISCOVERY/inputs/NORMAL_DONOR_STATE_COUNTS.tsv.gz'
contract={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'role':'Independent normal-reference maturation association, not clinical OS replication','prior_source_exposure':'Main article, Figure1 and survival summaries were read during source qualification; this is exploratory, not prospectively blinded confirmation. No ZEB1/state correlation coefficients have been viewed.','modules':'Exactly frozen200 normal donor-selected genes, no reselection; exclude ZEB1/ZEB2 and cell-cycle programmes as originally defined.','mapping':'Exact stable ENSG via eligible original mapping, remove ambiguous symbols/IDs; author symbols retained as provenance.','coverage_gate':'At least90 of100 members per axis including at least45 per direction; ranking universe at least20000 original shared symbols. Missing symbols remain missing; no substitute members.','score':'Within-sample average ranks / original-common intersected universe; mean frozen up genes minus mean frozen down genes. Raw scores; no cross-cohort z threshold.','test_family':'Two two-sided Spearman tests of ZEB1 expression vs state similarities, BH2; 10000 paired-patient bootstrap CI, seed20261010. Conditional/dichotomized/outcome models not part of this probe.','unit':'25 independent diagnostic patient expression columns; exclude3 named control columns before associations.','endpoint_boundary':'No OS times/MRD patient table; do not impute time from aggregate KM or median followup.','normalization_boundary':'Author supplementary DESeq2 expression values are floating point normalized counts; not raw counts/TPM. Rank representation mitigates scales but not platform/normalization bias.'}
amend_path=OUT/'OFFICIAL_ANNOTATION_AMENDMENT.json';amended=amend_path.exists();prefix='ANNOTATION_AMENDED_' if amended else ''
if amended:contract['mapping_amendment']=json.loads(amend_path.read_text())
(OUT/(prefix+'EXPLORATORY_PROBE_CONTRACT.json')).write_text(json.dumps(contract,indent=2),encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
book=openpyxl.load_workbook(xp,read_only=True,data_only=True);sh=book.active;rows=sh.iter_rows(values_only=True);header=list(next(rows));clinical=[str(x) for x in header[16:]];assert len(clinical)==25 and len(set(clinical))==25 and all(re.fullmatch(r'\d+',x) for x in clinical)
vals=[];ids=[];author_symbols=[]
for r in rows:
    stable=str(r[0]).split('.')[0]
    if not stable.startswith('ENSG'):continue
    try:v=np.asarray(r[16:],dtype=float)
    except (ValueError,TypeError):raise RuntimeError('Nonnumeric author expression entries for '+stable)
    assert len(v)==25 and np.isfinite(v).all() and (v>=0).all()
    ids.append(stable);author_symbols.append(str(r[1]));vals.append(v)
book.close();expr=pd.DataFrame(vals,index=ids,columns=clinical);assert expr.index.is_unique
mapping=pd.read_csv(mapping_path,sep='\t');mapping=mapping[mapping.eligible];assert mapping.stable.is_unique and mapping.symbol.is_unique
normal_symbols=set(pd.read_csv(normal_path,sep='\t',usecols=[0]).iloc[:,0].astype(str));approved=mapping[mapping.symbol.isin(normal_symbols)]
original_common_n=25288;common=approved[approved.stable.isin(expr.index)].copy();author_by_id=dict(zip(ids,author_symbols));common['identity_method']='EXACT_STABLE_ENSG'
if amended:
    official=json.loads((OUT/'MISSING_MEMBER_OFFICIAL_IDENTITY_QUALIFICATION.json').read_text());official_by_symbol={r['symbol']:r for r in official['sources']}
    for m in contract['mapping_amendment']['approved_same_gene_mappings']:
        o=official_by_symbol[m['symbol']];assert len(o['records'])==1 and o['records'][0]['hgnc_id']==m['HGNC'] and len(o['author_matches'])==1
        assert o['author_matches'][0]=={'stable':m['author_stable'],'author_symbol':m['author_previous_symbol']}
        assert m['author_stable'] not in set(common.stable) and m['symbol'] not in set(common.symbol) and m['author_stable'] in expr.index
        common=pd.concat([common,pd.DataFrame([{'stable':m['author_stable'],'symbol':m['symbol'],'eligible':True,'identity_method':'HGNC_UNIQUE_CURRENT_OR_PREVIOUS_SYMBOL'}])],ignore_index=True)
common['author_symbol']=[author_by_id[i] for i in common.stable];common['symbol_disagreement']=common.author_symbol!=common.symbol;common.to_csv(OUT/(prefix+'EXACT_ENSG_MAPPING_LEDGER.tsv'),sep='\t',index=False)
mapped=expr.loc[common.stable].copy();mapped.index=common.symbol.to_numpy();members=pd.read_csv(members_path,sep='\t');assert len(members)==200 and all(members.groupby(['axis','direction']).size()==50)
coverage=[]
for (axis,direction),g in members.groupby(['axis','direction']):
    present=[s for s in g.symbol if s in mapped.index];missing=[s for s in g.symbol if s not in mapped.index];coverage.append({'axis':axis,'direction':direction,'n_original':50,'n_present':len(present),'missing_symbols':';'.join(missing)})
cov=pd.DataFrame(coverage);cov.to_csv(OUT/(prefix+'FIXED_MODULE_INPUT_COVERAGE.tsv'),sep='\t',index=False)
passed=len(mapped)>=20000 and all(cov.n_present>=45) and all(cov.groupby('axis').n_present.sum()>=90)
q={'status':'QUALIFIED_FOR_EXPLORATORY_INDEPENDENT_MATURATION_ASSOCIATION' if passed else 'COVERAGE_GATE_FAILED_NO_ASSOCIATION_FITTED','n_patient_columns':25,'n_expression_records':len(expr),'n_exact_mapped_ranking_symbols':len(mapped),'reference_common_symbol_count':original_common_n,'mapped_symbol_disagreements_n':int(common.symbol_disagreement.sum()),'module_coverage':coverage,'sources':{str(p):sha(p) for p in [xp,members_path,mapping_path,normal_path,Path(__file__)]},'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'boundary':'Input qualification does not establish independent clinical OS or full maturation trajectory.'}
if amended:q['annotation_amendment_sha256']=sha(amend_path);q['official_identity_receipt_sha256']=sha(OUT/'MISSING_MEMBER_OFFICIAL_IDENTITY_QUALIFICATION.json')
(OUT/(prefix+'INPUT_QUALIFICATION.json')).write_text(json.dumps(q,indent=2),encoding='utf-8');print('INDIA25_INPUT_QUALIFICATION',q['status'],'ranking_symbols',len(mapped),'module_coverage',coverage,flush=True)
if not passed:raise SystemExit(0)
assert 'ZEB1' in mapped.index
ranks=mapped.rank(axis=0,method='average')/len(mapped);scores=pd.DataFrame(index=clinical);scores['ZEB1_author_normalized_counts']=mapped.loc['ZEB1'];scores['ZEB1_log2_counts_plus1']=np.log2(mapped.loc['ZEB1']+1);results=[];rng=np.random.default_rng(20261010)
for axis,g in members.groupby('axis'):
    up=[x for x in g.loc[g.direction.eq('up'),'symbol'] if x in mapped.index];down=[x for x in g.loc[g.direction.eq('down'),'symbol'] if x in mapped.index];sc=ranks.loc[up].mean()-ranks.loc[down].mean();scores[axis+'_raw']=sc
    x=scores.ZEB1_log2_counts_plus1.to_numpy();y=sc.reindex(clinical).to_numpy();rho,p=stats.spearmanr(x,y);boot=[]
    for _ in range(10000):
        ix=rng.integers(0,25,25);r=stats.spearmanr(x[ix],y[ix]).statistic
        if np.isfinite(r):boot.append(r)
    lo,hi=np.quantile(boot,[.025,.975]);results.append({'axis':axis,'n_patients':25,'spearman_rho':float(rho),'p_two_sided':float(p),'bootstrap95CI_low':float(lo),'bootstrap95CI_high':float(hi),'n_valid_bootstrap':len(boot),'members_up_used':len(up),'members_down_used':len(down)})
pvals=np.array([r['p_two_sided'] for r in results]);order=np.argsort(pvals);qs=np.minimum.accumulate((pvals[order]*len(pvals)/np.arange(1,len(pvals)+1))[::-1])[::-1].clip(max=1)
for i,z in zip(order,qs):results[i]['BH2_q']=float(z)
scores.rename_axis('patient_id').to_csv(OUT/'ALL25_PATIENT_FIXED_SCORES.tsv',sep='\t');pd.DataFrame(results).to_csv(OUT/'ALL2_ASSOCIATION_RESULTS.tsv',sep='\t',index=False)
r={'status':'EXPLORATORY_INDEPENDENT_EXPRESSION_ASSOCIATIONS_ESTIMATED','results':results,'n_reference_fixed_members':200,'n_ranking_symbols':len(mapped),'source_qualification_sha256':sha(OUT/(prefix+'INPUT_QUALIFICATION.json')),'contract_sha256':sha(OUT/(prefix+'EXPLORATORY_PROBE_CONTRACT.json')),'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'limitations':'Small independent cohort, raw rank-state scores, author normalized-count platform, no independent normal reference. No OS/MRD, no window-specific clinical effect, no mechanism. Preserve both axis outcomes without selecting significant model.'}
(OUT/'RESULTS.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print('INDIA25_FIXED_MATURATION_ASSOCIATIONS',json.dumps(results),flush=True)
