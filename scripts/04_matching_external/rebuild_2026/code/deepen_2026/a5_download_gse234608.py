from pathlib import Path
from urllib.parse import quote
import urllib.request

out = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a5_gse234608")
out.mkdir(parents=True, exist_ok=True)
dest = out / "GSE234608_raw_counts.csv.gz"
urls = [
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234608&format=file&file=" + quote("GSE234608_raw_counts.csv.gz"),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234608/suppl/GSE234608_raw_counts.csv.gz",
]
if dest.exists() and dest.stat().st_size > 1_000_000:
    print("exists", dest, dest.stat().st_size)
else:
    last = None
    for url in urls:
        print("GET", url, flush=True)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 a5-zeb1"})
            with urllib.request.urlopen(req, timeout=120) as r:
                dest.write_bytes(r.read())
            print("wrote", dest, dest.stat().st_size, flush=True)
            if dest.stat().st_size > 1_000_000:
                break
        except Exception as e:
            last = e
            print("fail", e, flush=True)
    else:
        raise SystemExit(last)
print("done")
