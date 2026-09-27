import json
from pathlib import Path

p = Path(r"C:\Users\Ls180\.cursor\projects\d-bioinformation\agent-tools\c6fd6d64-53e7-411e-ab24-0aa33cb0bf8f.txt")
out = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a4_yayon")
out.mkdir(parents=True, exist_ok=True)
js = json.loads(p.read_text(encoding="utf-8"))
rows = []
for d in js.get("datasets", []):
    title = d.get("title") or d.get("name")
    cells = d.get("cell_count")
    assets = d.get("assets") or []
    h5 = next((a for a in assets if str(a.get("filetype", "")).lower() in {"h5ad", "rds"} or str(a.get("url", "")).endswith(".h5ad")), None)
    url = (h5 or {}).get("url") if h5 else None
    size = (h5 or {}).get("filesize") if h5 else None
    rows.append((title, cells, size, url, d.get("dataset_id") or d.get("id")))
    print(f"{cells}\t{size}\t{title}", flush=True)
    if url:
        print("  ", url[:180], flush=True)
(out / "cellxgene_datasets.tsv").write_text(
    "title\tcells\tsize\tid\turl\n"
    + "\n".join(f"{t}\t{c}\t{s}\t{i}\t{u}" for t, c, s, u, i in rows),
    encoding="utf-8",
)
print("n_datasets", len(rows), flush=True)
