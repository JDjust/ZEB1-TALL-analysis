"""Download Yayon merged paediatric Visium h5ad from CELLxGENE (not NCBI)."""
from pathlib import Path
from urllib.request import Request, urlopen

out = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
out.mkdir(parents=True, exist_ok=True)
dest = out / "merged_thymus_visium_pediatric.h5ad"
url = "https://datasets.cellxgene.cziscience.com/9532bbf2-da9f-4cc1-bf52-df9a3dfea3d1.h5ad"
print("GET", url, flush=True)
req = Request(url, headers={"User-Agent": "a4-yayon"})
with urlopen(req, timeout=30) as r:
    # stream to disk; do not load into pandas
    n = 0
    with dest.open("wb") as f:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
            n += len(chunk)
            if n % (20 * 1024 * 1024) < 1024 * 1024:
                print("bytes", n, flush=True)
print("wrote", dest, dest.stat().st_size, flush=True)
