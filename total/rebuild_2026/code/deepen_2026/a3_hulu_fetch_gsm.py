#!/usr/bin/env python3
from pathlib import Path
import urllib.request

base = Path("/data-b/liangfuhua/zeb1_a3")
chip = base / "gse165016"
hic = base / "gse165207"
chip.mkdir(parents=True, exist_ok=True)
hic.mkdir(parents=True, exist_ok=True)

chip_files = [
    "GSM5024536_DND41_H3K27ac_peaks.narrowPeak.gz",
    "GSM5024537_DND41_BCL11B_peaks.narrowPeak.gz",
    "GSM5024539_SJTALL005006_H3K27ac_peaks.narrowPeak.gz",
    "GSM5024540_SJTALL005006_BCL11B_rep1_peaks.narrowPeak.gz",
    "GSM5024541_SJTALL005006_BCL11B_rep2_peaks.narrowPeak.gz",
    "GSM5265324_BCL11B_E15804_peaks.narrowPeak.gz",
    "GSM5265325_BCL11B_E22952_peaks.narrowPeak.gz",
    "GSM5265326_BCL11B_E14259_peaks.narrowPeak.gz",
]
bedpe = [
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


def fetch(name, dest_dir):
    dest = dest_dir / name
    if dest.exists() and dest.stat().st_size > 1000:
        print("exists", dest.name, dest.stat().st_size, flush=True)
        return
    gsm = name.split("_", 1)[0]
    nnn = gsm[:7] + "nnn"
    urls = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{nnn}/{gsm}/suppl/{name}",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gsm}&format=file&file={name}",
    ]
    last = None
    for url in urls:
        print("GET", url, flush=True)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "a3-zeb1/1.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                dest.write_bytes(r.read())
            if dest.stat().st_size > 1000:
                print(" wrote", dest.name, dest.stat().st_size, flush=True)
                return
        except Exception as e:
            last = e
            print(" fail", e, flush=True)
    raise RuntimeError(f"could not fetch {name}: {last}")


for n in chip_files:
    fetch(n, chip)
for n in bedpe:
    fetch(n, hic)
print("ALL_DONE", flush=True)
