"""Validate every final-delivery file against its recorded SHA-256 digest."""
from pathlib import Path
import csv
import hashlib

base=Path(__file__).resolve().parents[1]
manifest=base/'data/manifests/final_delivery_manifest.tsv'
with manifest.open(encoding='utf-8-sig',newline='') as fh:
    rows=list(csv.DictReader(fh,delimiter='\t'))
errors=[]
for row in rows:
    path=Path(row['path'])
    if not path.is_absolute(): path=base/path
    if not path.is_file():
        errors.append(f'Missing: {path}')
        continue
    if path.stat().st_size!=int(row['bytes']):
        errors.append(f'Size mismatch: {path}')
        continue
    with path.open('rb') as fh:
        digest=hashlib.file_digest(fh,'sha256').hexdigest()
    if digest!=row['sha256']:
        errors.append(f'Hash mismatch: {path}')
print(f'Checked {len(rows)} delivered files; mismatches: {len(errors)}')
for error in errors: print(error)
if errors: raise SystemExit(1)
