"""Resolve every figure-manifest upstream source to an existing frozen input file."""
from pathlib import Path
import csv
import hashlib

base=Path(__file__).resolve().parents[1]
project=base.parents[1]
validation=project/'total/data/validation'
figures=list((base/'figures/main').glob('*_panel_manifest.tsv'))
figures+=list((base/'figures/supplementary').glob('*_panel_manifest.tsv'))
out=base/'data/manifests/panel_source_resolution.tsv'
rows=[]
errors=[]
for manifest in sorted(figures):
    with manifest.open(encoding='utf-8-sig',newline='') as fh:
        panels=list(csv.DictReader(fh,delimiter='\t'))
    for panel in panels:
        for token in panel['source'].split(' + '):
            token=token.strip()
            candidates=[project/token,validation/token]
            matches=[x for x in candidates if x.is_file()]
            if not matches and '/' not in token and '\\' not in token:
                matches=list(validation.rglob(token))
            if len(matches)!=1:
                errors.append(f'{manifest.name} {panel["panel"]}: {token} -> {len(matches)} matches')
                continue
            path=matches[0].resolve()
            with path.open('rb') as fh:
                digest=hashlib.file_digest(fh,'sha256').hexdigest()
            rows.append(dict(figure=panel['figure'],panel=panel['panel'],
                             manifest=str(manifest.relative_to(base)),raw_source=token,
                             resolved_source=str(path.relative_to(project)),bytes=path.stat().st_size,
                             sha256=digest))
out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',encoding='utf-8',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ['figure','panel','manifest','raw_source','resolved_source','bytes','sha256'],delimiter='\t')
    w.writeheader();w.writerows(rows)
print(f'Resolved {len(rows)} upstream panel-source references; ambiguous or missing: {len(errors)}')
for error in errors: print(error)
if errors: raise SystemExit(1)
