from pathlib import Path
from urllib.request import Request, urlopen

out = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a7_lineage")
out.mkdir(parents=True, exist_ok=True)
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253355/suppl/"
print("GET", url, flush=True)
req = Request(url, headers={"User-Agent": "a7-list"})
with urlopen(req, timeout=15) as r:
    html = r.read().decode("utf-8", "replace")
(out / "GSE253355_suppl.html").write_text(html, encoding="utf-8")
for line in html.splitlines():
    if any(k in line for k in ("rds", "Seurat", "h5ad", "MSC", "Atlas")):
        print(line[:300], flush=True)
