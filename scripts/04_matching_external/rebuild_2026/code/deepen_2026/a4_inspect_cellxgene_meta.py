import json
from pathlib import Path

p = Path(r"C:\Users\Ls180\.cursor\projects\d-bioinformation\agent-tools\c6fd6d64-53e7-411e-ab24-0aa33cb0bf8f.txt")
js = json.loads(p.read_text(encoding="utf-8"))
want = ("T cell subset", "merged thymus visium pediatric", "merged thymus visium fetal", "thymus scRNA-seq atlas")
for d in js.get("datasets", []):
    title = d.get("title") or d.get("name") or ""
    if not any(w in title for w in want):
        continue
    print("=" * 60, flush=True)
    print(title, flush=True)
    print("id", d.get("dataset_id") or d.get("id"), "cells", d.get("cell_count"), flush=True)
    for key in ("obs", "schema", "annotations", "extra_layers"):
        if key in d:
            print(key, str(d[key])[:400], flush=True)
    keys = sorted(d.keys())
    print("keys", keys, flush=True)
    for k in keys:
        if k in {"title", "name", "assets", "dataset_id", "id", "cell_count", "citation", "description", "mean_genes_per_cell"}:
            continue
        val = d[k]
        s = str(val)
        if len(s) < 500:
            print(k, s, flush=True)
        elif "cell_type" in k.lower() or "obs" in k.lower() or "annot" in k.lower():
            print(k, s[:400], flush=True)
