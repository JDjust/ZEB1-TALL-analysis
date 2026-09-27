from pathlib import Path
import re
import urllib.request

out = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a4_yayon")
out.mkdir(parents=True, exist_ok=True)
for acc in ["GSE271304", "GSE247917", "GSE247918"]:
    nnn = "GSE" + acc[3:6] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/suppl/"
    print("====", acc, url, flush=True)
    try:
        html = urllib.request.urlopen(url, timeout=40).read().decode("utf-8", "replace")
        files = [x for x in re.findall(r'href="([^"]+)"', html) if not x.startswith("?")]
        print("\n".join("  " + x for x in files), flush=True)
        (out / f"{acc}_suppl.txt").write_text("\n".join(files), encoding="utf-8")
    except Exception as e:
        print("FAIL", e, flush=True)
    try:
        fl = urllib.request.urlopen(url + "filelist.txt", timeout=40).read().decode("utf-8", "replace")
        (out / f"{acc}_filelist.txt").write_text(fl, encoding="utf-8")
        print(fl[:2000], flush=True)
    except Exception as e:
        print("no filelist", e, flush=True)
