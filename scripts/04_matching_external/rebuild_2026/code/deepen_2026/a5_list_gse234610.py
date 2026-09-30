from pathlib import Path
from urllib.request import Request, urlopen

out = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a5_gse234610")
out.mkdir(parents=True, exist_ok=True)
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234610/suppl/"
print("GET", url, flush=True)
req = Request(url, headers={"User-Agent": "a5-list"})
with urlopen(req, timeout=15) as r:
    html = r.read().decode("utf-8", "replace")
(out / "GSE234610_suppl.html").write_text(html, encoding="utf-8")
print("bytes", len(html), flush=True)
for line in html.splitlines():
    if any(k in line.lower() for k in ("h5", "rds", "h5ad", "meta", "annot", "seurat", "csv", "tsv", "gz")):
        print(line[:240], flush=True)
