"""Download Yayon T-cell atlas h5ad from CELLxGENE for A4-A stage labels."""
from pathlib import Path
from urllib.request import Request, urlopen

out = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a4_yayon")
dest = out / "thymus_scrna_tcell_subset.h5ad"
url = "https://datasets.cellxgene.cziscience.com/9e93e320-fecf-4f29-88b0-f052c945a335.h5ad"
print("GET", url, flush=True)
req = Request(url, headers={"User-Agent": "a4-yayon"})
with urlopen(req, timeout=30) as r:
    n = 0
    with dest.open("wb") as f:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
            n += len(chunk)
            if n % (100 * 1024 * 1024) < 1024 * 1024:
                print("bytes", n, flush=True)
print("wrote", dest, dest.stat().st_size, flush=True)
