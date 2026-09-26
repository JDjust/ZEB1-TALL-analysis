"""Make a local, reviewable release candidate; this script publishes nothing."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import csv
import hashlib

base=Path(__file__).resolve().parents[1]
project=base.parents[1]
out=base/'data/manifests/reviewer_release_candidate.zip'
sources=base/'data/manifests/panel_source_resolution.tsv'
with sources.open(encoding='utf-8-sig',newline='') as fh:
    rows=list(csv.DictReader(fh,delimiter='\t'))
files=set()
for r in rows:
    files.add(project/Path(r['resolved_source']))
for name in ('model_induction_failure.tsv','model_poor_morphology.tsv',
             'model_mrd.tsv','model_survival.tsv'):
    files.add(project/'total/data/validation/polonen_round1'/name)
files.update((base/'code').glob('*.R'))
files.update((base/'code').glob('*.py'))
for obsolete in ('build_f3_ggplot.R','build_s1_table_pdf.py',
                 'build_s4_table_pdf.py'):
    files.discard(base/'code'/obsolete)
files.update((base/'data/source_data_rebuilt').rglob('*.tsv'))
files.update((base/'data/source_data_rebuilt/revision2_statistics').glob('*.txt'))
files.update((base/'figures/main').glob('*.pdf'))
files.update((base/'figures/supplementary').glob('*.pdf'))
files.update((base/'figures/main').glob('*_panel_manifest.tsv'))
files.update((base/'figures/supplementary').glob('*_panel_manifest.tsv'))
files.update((base/'manuscript').glob('*.md'))
files.update((base/'manuscript').glob('*.pdf'))
files.update((base/'data/manifests').glob('panel_source_resolution.tsv'))
files.add(base/'README.md')
missing=[p for p in files if not p.is_file()]
if missing:
    raise FileNotFoundError(missing)
note='''# Reviewer release candidate (local; not published)

This archive preserves the ZEB1 directory layout required by the R plotting
scripts. It contains frozen source tables used by figure panels, the Revision 2
patient-bootstrap input/output, manuscript and supplement PDFs/Markdown, all
14 figure PDFs, source-data tables and active production scripts. The original
raw public study archives and author-only metadata are not redistributed here.

Run from `ZEB1/total/rebuild_2026` using the order in README.md. Figure 3
primary inference is within the 1,309-patient diagnostic cohort. Its normal
thymus projection is a separately standardized sensitivity only. The bundle
needs author review and an approved public/reviewer-accessible deposit before
the manuscript availability statement can be finalized.
'''
entries=[]
with ZipFile(out,'w',compression=ZIP_DEFLATED,compresslevel=7) as z:
    z.writestr('REVIEWER_RELEASE_README.md',note)
    for p in sorted(files,key=lambda x:str(x).lower()):
        arc='ZEB1/'+p.relative_to(project).as_posix()
        z.write(p,arcname=arc)
        entries.append((arc,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()))
manifest=base/'data/manifests/reviewer_release_candidate_contents.tsv'
with manifest.open('w',encoding='utf-8',newline='') as fh:
    w=csv.writer(fh,delimiter='\t')
    w.writerow(['archive_path','bytes','sha256'])
    w.writerows(entries)
print(f'{out} | files={len(entries)} | bytes={out.stat().st_size}')
