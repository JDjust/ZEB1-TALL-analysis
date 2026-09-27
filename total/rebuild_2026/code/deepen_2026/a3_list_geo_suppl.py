"""List GEO series supplement files. Do not download the 23 GB SuperSeries RAW tar."""
from pathlib import Path
import re
import urllib.request

OUT = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a3_gse165209")
OUT.mkdir(parents=True, exist_ok=True)

for acc in ["GSE173432", "GSE165016", "GSE165207"]:
    nnn = "GSE" + acc[3:6] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/suppl/"
    print("====", acc, url)
    try:
        html = urllib.request.urlopen(url, timeout=40).read().decode("utf-8", "replace")
        links = re.findall(r'href="([^"]+)"', html)
        files = [x for x in links if not x.startswith("?") and x not in ("/", "../")]
        for x in files:
            print(" ", x)
        (OUT / f"{acc}_suppl_listing.txt").write_text("\n".join(files), encoding="utf-8")
    except Exception as e:
        print("FAIL", acc, e)
