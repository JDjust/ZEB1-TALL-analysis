"""Download only narrowPeak and BEDPE from GSM sample FTP. Skip RAW/HIC/BW."""
from pathlib import Path
import urllib.request

ROOT = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a3_gse165209")
CHIP = ROOT / "gse165016"
HIC = ROOT / "gse165207"
CHIP.mkdir(parents=True, exist_ok=True)
HIC.mkdir(parents=True, exist_ok=True)

CHIP_FILES = [
    "GSM5024536_DND41_H3K27ac_peaks.narrowPeak.gz",
    "GSM5024537_DND41_BCL11B_peaks.narrowPeak.gz",
    "GSM5024539_SJTALL005006_H3K27ac_peaks.narrowPeak.gz",
    "GSM5024540_SJTALL005006_BCL11B_rep1_peaks.narrowPeak.gz",
    "GSM5024541_SJTALL005006_BCL11B_rep2_peaks.narrowPeak.gz",
    "GSM5265324_BCL11B_E15804_peaks.narrowPeak.gz",
    "GSM5265325_BCL11B_E22952_peaks.narrowPeak.gz",
    "GSM5265326_BCL11B_E14259_peaks.narrowPeak.gz",
]
BEDPE_FILES = [
    "GSM5028224_H3K27ac_HiChIP_loops_SJALL068279.bedpe.gz",
    "GSM5028225_H3K27ac_HiChIP_loops_SJAUL068292.bedpe.gz",
    "GSM5028226_H3K27ac_HiChIP_loops_SJMPAL011914.bedpe.gz",
    "GSM5028227_H3K27ac_HiChIP_loops_SJTALL005006.bedpe.gz",
    "GSM5028228_H3K27ac_HiChIP_loops_SJMPAL011911.bedpe.gz",
    "GSM5028229_H3K27ac_HiChIP_loops_DND41_rep1.bedpe.gz",
    "GSM5028230_H3K27ac_HiChIP_loops_DND41_rep2.bedpe.gz",
    "GSM5028231_H3K27ac_HiChIP_loops_Jurkat.bedpe.gz",
    "GSM5028232_H3K27ac_HiChIP_loops_CD34_5M.bedpe.gz",
]


def gsm_from_name(name):
    return name.split("_", 1)[0]


def fetch(name, dest_dir):
    dest = dest_dir / name
    if dest.exists() and dest.stat().st_size > 1000:
        print("exists", dest.name, dest.stat().st_size, flush=True)
        return
    gsm = gsm_from_name(name)
    nnn = gsm[:7] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{nnn}/{gsm}/suppl/{name}"
    print("GET", url, flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "a3-zeb1/1.0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        dest.write_bytes(r.read())
    print(" wrote", dest.name, dest.stat().st_size, flush=True)


def main():
    for name in CHIP_FILES:
        fetch(name, CHIP)
    for name in BEDPE_FILES:
        fetch(name, HIC)
    print("ALL_DONE", flush=True)


if __name__ == "__main__":
    main()
