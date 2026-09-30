from pathlib import Path
import gzip
import urllib.request

out = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a5_gse234608")
out.mkdir(parents=True, exist_ok=True)
urls = [
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234608/matrix/GSE234608_series_matrix.txt.gz",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234608&targ=self&form=text&view=full",
]
for url in urls:
    dest = out / Path(url.split("?")[0]).name
    if "acc.cgi" in url:
        dest = out / "GSE234608_geo_soft.txt"
    print("GET", url, flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 a5-zeb1"})
        with urllib.request.urlopen(req, timeout=60) as r:
            dest.write_bytes(r.read())
        print("wrote", dest, dest.stat().st_size, flush=True)
        if dest.suffix == ".gz" or dest.name.endswith(".gz"):
            with gzip.open(dest, "rt", encoding="utf-8", errors="replace") as f:
                for i, line in enumerate(f):
                    if i > 80:
                        break
                    if line.startswith("!") or "ETP" in line or "MPAL" in line or "T-ALL" in line or "characteristic" in line.lower():
                        print(line[:300].rstrip())
        else:
            text = dest.read_text(encoding="utf-8", errors="replace")
            for line in text.splitlines()[:120]:
                if "Sample_title" in line or "characteristic" in line.lower() or "ETP" in line or "MPAL" in line:
                    print(line[:300])
    except Exception as e:
        print("FAIL", e, flush=True)
