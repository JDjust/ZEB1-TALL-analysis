from pathlib import Path
import re
import urllib.request

names = [
    "GSM5024536", "GSM5024537", "GSM5024539", "GSM5024540", "GSM5024541",
    "GSM5265324", "GSM5265325", "GSM5265326",
    "GSM5028224", "GSM5028225", "GSM5028226", "GSM5028227", "GSM5028228",
    "GSM5028229", "GSM5028230", "GSM5028231", "GSM5028232",
]
out = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a3_gse165209\gsm_suppl_listing.txt")
lines = []
for gsm in names:
    nnn = gsm[:7] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{nnn}/{gsm}/suppl/"
    print("====", gsm)
    try:
        html = urllib.request.urlopen(url, timeout=30).read().decode("utf-8", "replace")
        files = [x for x in re.findall(r'href="([^"]+)"', html) if x.endswith((".gz", ".bedpe", ".narrowPeak", ".txt"))]
        print(" ", files)
        lines.append(gsm + "\t" + "\t".join(files))
    except Exception as e:
        print(" FAIL", e)
        lines.append(gsm + "\tFAIL\t" + str(e))
out.write_text("\n".join(lines), encoding="utf-8")
