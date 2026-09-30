from pathlib import Path
from urllib.parse import quote
import urllib.request

ROOT = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a3_gse165209")
jobs = {
    ROOT / "gse165016": [
        "GSM5024536_DND41_H3K27ac_peaks.narrowPeak.gz",
        "GSM5024537_DND41_BCL11B_peaks.narrowPeak.gz",
        "GSM5024539_SJTALL005006_H3K27ac_peaks.narrowPeak.gz",
        "GSM5024540_SJTALL005006_BCL11B_rep1_peaks.narrowPeak.gz",
        "GSM5024541_SJTALL005006_BCL11B_rep2_peaks.narrowPeak.gz",
        "GSM5265324_BCL11B_E15804_peaks.narrowPeak.gz",
        "GSM5265325_BCL11B_E22952_peaks.narrowPeak.gz",
        "GSM5265326_BCL11B_E14259_peaks.narrowPeak.gz",
    ],
    ROOT / "gse165207": [
        "GSM5028224_H3K27ac_HiChIP_loops_SJALL068279.bedpe.gz",
        "GSM5028225_H3K27ac_HiChIP_loops_SJAUL068292.bedpe.gz",
        "GSM5028226_H3K27ac_HiChIP_loops_SJMPAL011914.bedpe.gz",
        "GSM5028227_H3K27ac_HiChIP_loops_SJTALL005006.bedpe.gz",
        "GSM5028228_H3K27ac_HiChIP_loops_SJMPAL011911.bedpe.gz",
        "GSM5028229_H3K27ac_HiChIP_loops_DND41_rep1.bedpe.gz",
        "GSM5028230_H3K27ac_HiChIP_loops_DND41_rep2.bedpe.gz",
        "GSM5028231_H3K27ac_HiChIP_loops_Jurkat.bedpe.gz",
        "GSM5028232_H3K27ac_HiChIP_loops_CD34_5M.bedpe.gz",
    ],
}


def main():
    for dest_dir, names in jobs.items():
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            dest = dest_dir / name
            if dest.exists() and dest.stat().st_size > 1000:
                print("exists", name, dest.stat().st_size, flush=True)
                continue
            gsm = name.split("_", 1)[0]
            url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gsm}&format=file&file={quote(name)}"
            print("GET", url, flush=True)
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 a3-zeb1"})
            with urllib.request.urlopen(req, timeout=90) as r:
                dest.write_bytes(r.read())
            print(" wrote", name, dest.stat().st_size, flush=True)
    print("ALL_DONE", flush=True)


if __name__ == "__main__":
    main()
