"""Download only processed ChIP peaks and HiChIP BEDPE. Never the SuperSeries RAW tar."""
from pathlib import Path
import re
import urllib.request

ROOT = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a3_gse165209")
CHIP = ROOT / "gse165016"
HIC = ROOT / "gse165207"
CHIP.mkdir(parents=True, exist_ok=True)
HIC.mkdir(parents=True, exist_ok=True)


def get(url, dest: Path):
    print("GET", url)
    urllib.request.urlretrieve(url, dest)
    print(" wrote", dest, dest.stat().st_size)


def list_suppl(acc: str) -> list[str]:
    nnn = "GSE" + acc[3:6] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/suppl/filelist.txt"
    dest = ROOT / f"{acc}_filelist.txt"
    get(url, dest)
    return dest.read_text(encoding="utf-8", errors="replace").splitlines()


def gsm_processed(acc: str) -> list[tuple[str, str, str]]:
    """Return (gsm, filename, url) for sample-level processed files."""
    # GEO SOFT family is large; use the series matrix / FTP samples listing instead.
    nnn = "GSM" + acc[3:6] + "nnn"  # unused placeholder
    return []


def main():
    for acc, outdir, keep_re in [
        ("GSE165016", CHIP, re.compile(r"(narrowPeak|broadPeak|bed\.gz|peaks)", re.I)),
        ("GSE165207", HIC, re.compile(r"(bedpe|loops|hiccups|FitHiChIP|mango)", re.I)),
    ]:
        rows = list_suppl(acc)
        print("====", acc, "filelist rows", len(rows))
        for line in rows[:8]:
            print(" ", line)
        # filelist.txt format: name<TAB>size
        wanted = []
        for line in rows:
            name = line.split("\t")[0].split()[0] if line.strip() else ""
            if keep_re.search(name) and not name.lower().endswith(".hic") and "bigwig" not in name.lower() and not name.lower().endswith(".bw"):
                wanted.append(name)
        print("wanted", acc, wanted)
        nnn = "GSE" + acc[3:6] + "nnn"
        for name in wanted:
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/suppl/{name}"
            dest = outdir / Path(name).name
            if dest.exists() and dest.stat().st_size > 0:
                print("exists", dest)
                continue
            try:
                get(url, dest)
            except Exception as e:
                print("FAIL", name, e)


if __name__ == "__main__":
    main()
