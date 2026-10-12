from pathlib import Path
import subprocess,json,hashlib,datetime,pandas as pd
from statsmodels.stats.multitest import multipletests
HERE=Path(__file__).resolve().parent;OUT=HERE/'multigene_maturity_proxy';amendment=HERE/'MULTIGENE_MATURITY_JOINT_AMENDMENT.json';primary=OUT/'CLINICAL_RESULTS.json';source=OUT/'MULTIGENE_CLINICAL_SENSITIVITY.R'
code=source.read_text();assert code.count('+ SP_squared + strata(subtype_g)')==1
code=code.replace('+ SP_squared + strata(subtype_g)','+ SP_squared + dev_centered + dev_squared + strata(subtype_g)')
assert code.count('~ZEB1+cortical_vs_early+cortical_squared+SP_vs_cortical+SP_squared,')==1
code=code.replace('~ZEB1+cortical_vs_early+cortical_squared+SP_vs_cortical+SP_squared,','~ZEB1+cortical_vs_early+cortical_squared+SP_vs_cortical+SP_squared+dev_centered+dev_squared,')
script=OUT/'JOINT_MATURITY_SENSITIVITY.R';script.write_text(code,encoding='utf-8');print('JOINT_MATURITY_SENSITIVITY_START_UNCHANGED_LANDMARK_PATIENTS',flush=True)
p=subprocess.run([r'D:\Git\Git\usr\bin\ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15','hulu','/opt/bioinfo/envs/r-bio/bin/Rscript','-'],input=code,capture_output=True,text=True,encoding='utf-8',timeout=240)
(OUT/'JOINT_CLINICAL_REMOTE_LOG.txt').write_text(p.stdout+'\nSTDERR\n'+p.stderr,encoding='utf-8')
if p.returncode:raise RuntimeError(p.stderr[-2500:])
result=json.loads(p.stdout.split('RESULT_BEGIN\n',1)[1].split('\nRESULT_END',1)[0]);initial=json.loads(primary.read_text())
family=[(label,ep,obj['results'][ep]['ZEB1_block_LRT']['p']) for label,obj in [('primary',initial),('joint_sensitivity',result)] for ep in ['OS','EFS']]
q=multipletests([r[2] for r in family],method='fdr_bh')[1];result['BH4_FAMILY']=[{'model':label,'endpoint':ep,'p':pv,'BH_q_four':float(qq)} for (label,ep,pv),qq in zip(family,q)]
for ep,obj in result['results'].items():
    for name in ['static_reduced','static_full','extended_reduced','extended_constant','extended_full']:pd.DataFrame(obj[name]['coefficients']).to_csv(OUT/f'{ep}_JOINT_{name}_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(obj['static_PH']).to_csv(OUT/f'{ep}_JOINT_STATIC_PH.tsv',sep='\t',index=False)
    print('JOINT_MATURITY_RESULT',ep,'static_ZEB1',[x for x in obj['static_full']['coefficients'] if x['term']=='ZEB1'],'block',obj['ZEB1_block_LRT'],'warnings',obj['extended_full']['warnings'],'condition',obj['maturity_design_condition_number'],flush=True)
result.update(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [amendment,primary,source,script,Path(__file__)]},created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),interpretation='Exploratory post-primary joint-proxy sensitivity, no model selection; preserves original primary results. All6 maturity nuisance terms included, still measured-proxy and same-cohort limits.')
(OUT/'JOINT_CLINICAL_RESULTS.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print('JOINT_MATURITY_BH4',result['BH4_FAMILY'],flush=True)
