"""Exact bigWig window summaries and locus bins; no occupancy specificity inference."""
import argparse,json,re,hashlib
from pathlib import Path
import numpy as np,pandas as pd,pyBigWig
p=argparse.ArgumentParser();p.add_argument('--project',required=True);p.add_argument('--output',required=True);a=p.parse_args()
root=Path(a.project);out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
tss0=31319216-1  # Ensembl canonical transcript start, convert to zero-based.
wins={'legacy_12kb':(31316447,31328447),'canonical_TSS_pm1kb':(tss0-1000,tss0+1000),
      'canonical_TSS_pm2kb':(tss0-2000,tss0+2000),'legacy_gene_body':(31318447,31529814)}
rows=[];bins=[];manifest=[]
files=sorted((root/'module6/processed/gse154675_bw').glob('*.bigwig'))
for f in files:
    m=re.search(r'(ARR|DU528|HSB2|CCRFCEM)_a(LMO2|TAL1|LDB1|GATA2)',f.name)
    if not m:continue
    bw=pyBigWig.open(str(f));chrom='chr10' if 'chr10' in bw.chroms() else '10'
    length=bw.chroms(chrom)
    # Length is an assembly consistency check, not a substitute for GEO processing metadata.
    if length!=133797422:raise ValueError(f'{f.name}: chr10 length {length} incompatible with GRCh38')
    manifest.append({'file':f.name,'bytes':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns,'chr10_length':length})
    for window,(start,end) in wins.items():
        v=np.array(bw.values(chrom,start,end));finite=np.isfinite(v)
        rows.append(dict(file=f.name,line=m[1],antibody=m[2],window=window,start0=start,end0=end,
                         covered_fraction=float(finite.mean()),mean_covered=float(v[finite].mean()) if finite.any() else np.nan,
                         mean_missing_as_zero=float(np.nansum(v)/(end-start)),
                         mean_zoom_default=bw.stats(chrom,start,end,type='mean')[0]))
    for start in range(tss0-10000,tss0+10000,200):
        v=np.array(bw.values(chrom,start,start+200));finite=np.isfinite(v)
        bins.append(dict(file=f.name,line=m[1],antibody=m[2],start0=start,end0=start+200,
                         relative_mid_bp=start+100-tss0,mean_missing_as_zero=float(np.nansum(v)/200),covered_fraction=float(finite.mean())))
    bw.close()
assert len(manifest)==16
pd.DataFrame(rows).to_csv(out/'window_summaries.tsv',sep='\t',index=False)
pd.DataFrame(bins).to_csv(out/'locus_200bp_bins.tsv',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(out/'bigwig_manifest.tsv',sep='\t',index=False)
old=pd.read_csv(root/'module6/tables/M6_r2_GSE154675_occupancy.tsv',sep='\t')
compare=pd.DataFrame(rows).query('window == "legacy_12kb"').merge(old[['file','hg38_prom']],on='file',validate='one_to_one')
compare['default_abs_diff']=abs(compare.mean_zoom_default-compare.hg38_prom)
compare.to_csv(out/'legacy_comparison.tsv',sep='\t',index=False)
record={'canonical_TSS_1based':tss0+1,'coordinate_system':'0-based half-open for bigWig',
        'windows':wins,'bigwig_count':len(manifest),'pyBigWig_version':pyBigWig.__version__,
        'legacy_default_max_abs_diff':float(compare.default_abs_diff.max()),
        'limitations':'No matched input normalization; overlapping ZEB1-AS1 region; metadata manifest is not whole-file hashing; no causal or locus-exclusive binding claim.'}
(out/'run.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2),flush=True)
