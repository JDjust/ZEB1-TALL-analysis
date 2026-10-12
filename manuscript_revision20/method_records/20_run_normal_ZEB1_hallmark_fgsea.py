from pathlib import Path
import json,hashlib,datetime,subprocess
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];OUT=HERE/'normal_human_ZEB1_perturbation'
assert json.loads((OUT/'VECTOR_READBACK.json').read_text())['passed']
contract=HERE/'NORMAL_HUMAN_ZEB1_PROGRAMME_CONTRACT.json'
gmt=ROOT/'.project-artifacts/work/20261006T010250Z-75920166/analysis/Q2_FUNCTIONAL/inputs/h.all.v2025.1.Hs.symbols.gmt'
sets={}
for line in gmt.read_text(encoding='utf-8').splitlines():
    fields=line.split('\t');sets[fields[0]]=fields[2:]
assert len(sets)==50 and all(k.startswith('HALLMARK_') for k in sets)
jobs=[];sources=[gmt,contract,OUT/'EXTRACTION_RECEIPT.json',OUT/'VECTOR_READBACK.json'];coverage=[]
for condition in ['Rest','Stim8hr','Stim48hr']:
    source=OUT/f'ZEB1_{condition}_ALL_GENE_DE.tsv.gz';sources.append(source)
    df=pd.read_csv(source,sep='\t')
    rank=df.loc[np.isfinite(df.zscore)&~df.gene_name.eq('ZEB1'),['gene_name','zscore']].sort_values(['zscore','gene_name'],ascending=[False,True])
    assert rank.gene_name.is_unique
    rank.to_csv(OUT/f'{condition}_FULL_WALD_RANK.tsv',sep='\t',index=False)
    measured=set(rank.gene_name)
    for name,members in sets.items():coverage.append({'condition':condition,'pathway':name,'original_size':len(set(members)),'measured_size':len(set(members)&measured)})
    csv=rank.to_csv(sep='\t',index=False)
    jobs.append('data <- read.delim(textConnection('+json.dumps(csv,ensure_ascii=False)+'),check.names=FALSE)\nrank <- setNames(data$zscore,data$gene_name)\nresults[["'+condition+'"]] <- analyze(rank)')
paths='pathways <- list('+',\n'.join(json.dumps(name)+'=c('+','.join(json.dumps(g) for g in members)+')' for name,members in sets.items())+')\n'
r_code=r'''suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(jsonlite))
analyze <- function(rank){
 warnings <- character();set.seed(20261010)
 result <- withCallingHandlers(fgseaMultilevel(pathways=pathways,stats=rank,minSize=15,maxSize=500,eps=0,sampleSize=101,nPermSimple=10000,nproc=1),warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
 result$leadingEdge <- vapply(result$leadingEdge,paste,collapse=";",character(1))
 list(table=as.data.frame(result),warnings=warnings,n_ranked_genes=length(rank))
}
results <- list()
'''+paths+'\n'.join(jobs)
r_code+='\ncat("RESULT_BEGIN\\n");cat(toJSON(list(environment=list(R=R.version.string,fgsea=as.character(packageVersion("fgsea"))),results=results),auto_unbox=TRUE,digits=16,na="null"));cat("\\nRESULT_END\\n")\n'
script=OUT/'normal_ZEB1_hallmark_fgsea.R';script.write_text(r_code,encoding='utf-8')
print('FGSEA_START','3conditions','50unchanged_Hallmarks','all10281_nonZEB1_genes','fgseaMultilevel_eps0_nPermSimple10000_seed20261010',flush=True)
p=subprocess.run([r'D:\Git\Git\usr\bin\ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15','hulu','/opt/bioinfo/envs/r-bio/bin/Rscript','-'],input=r_code,capture_output=True,text=True,encoding='utf-8',timeout=240)
(OUT/'FGSEA_REMOTE_LOG.txt').write_text(p.stdout+'\nSTDERR\n'+p.stderr,encoding='utf-8')
if p.returncode:raise RuntimeError(p.stderr[-3000:])
result=json.loads(p.stdout.split('RESULT_BEGIN\n',1)[1].split('\nRESULT_END',1)[0])
tables=[]
for condition,record in result['results'].items():
    table=pd.DataFrame(record['table']);table['condition']=condition
    table=table.rename(columns={'padj':'author_method_within_condition_padj'})
    tables.append(table)
combined=pd.DataFrame(coverage).merge(pd.concat(tables,ignore_index=True),on=['condition','pathway'],how='left',validate='one_to_one')
combined['status']=np.where(combined.pval.notna(),'ESTIMATED','UNESTIMATED_COVERAGE_OR_METHOD')
assert len(combined)==150
q=multipletests(combined.pval.fillna(1).to_numpy(),method='fdr_bh')[1]
combined['BH_q_full150']=np.where(combined.pval.notna(),q,np.nan)
combined.to_csv(OUT/'FULL150_HALLMARK_RESULTS.tsv',sep='\t',index=False)
result['full_family_n']=150;result['n_estimated']=int(combined.pval.notna().sum())
result['BH_family']='All50 original sets x3 conditions; missing estimates reserved with p=1 only for conservative family denominator, actual p/q left unavailable'
result['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources+[script,Path(__file__)]}
result['created_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
result['key_programmes']=combined[combined.pathway.isin(['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1','HALLMARK_MYC_TARGETS_V2'])].to_dict('records')
result['all_full_family_q_lt_005']=combined.loc[combined.BH_q_full150.lt(.05)].to_dict('records')
(OUT/'HALLMARK_RESULTS.json').write_text(json.dumps(result,indent=2,default=str),encoding='utf-8')
print('FIXED_MAINLINE_PROGRAMMES',json.dumps([{k:r[k] for k in ['condition','pathway','NES','pval','BH_q_full150','measured_size','status']} for r in result['key_programmes']]),flush=True)
print('FGSEA_COMPLETE','n_estimated',result['n_estimated'],'family150','warnings',{k:v['warnings'] for k,v in result['results'].items()},'n_q005',len(result['all_full_family_q_lt_005']),flush=True)
