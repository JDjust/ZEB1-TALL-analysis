from pathlib import Path
import gzip
import urllib.request

ROOT = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a3_gse165209")
jobs = [
    (ROOT / "gse165016" / "GSM5265325_BCL11B_E22952_peaks.narrowPeak.gz", 572747),
    (ROOT / "gse165016" / "GSM5265326_BCL11B_E14259_peaks.narrowPeak.gz", 670468),
    (ROOT / "gse165207" / "GSM5028224_H3K27ac_HiChIP_loops_SJALL068279.bedpe.gz", 117852),
    (ROOT / "gse165207" / "GSM5028225_H3K27ac_HiChIP_loops_SJAUL068292.bedpe.gz", 204305),
    (ROOT / "gse165207" / "GSM5028226_H3K27ac_HiChIP_loops_SJMPAL011914.bedpe.gz", 91312),
    (ROOT / "gse165207" / "GSM5028227_H3K27ac_HiChIP_loops_SJTALL005006.bedpe.gz", 176000),
    (ROOT / "gse165207" / "GSM5028228_H3K27ac_HiChIP_loops_SJMPAL011911.bedpe.gz", 122315),
    (ROOT / "gse165207" / "GSM5028229_H3K27ac_HiChIP_loops_DND41_rep1.bedpe.gz", 343254),
    (ROOT / "gse165207" / "GSM5028230_H3K27ac_HiChIP_loops_DND41_rep2.bedpe.gz", 298472),
    (ROOT / "gse165207" / "GSM5028231_H3K27ac_HiChIP_loops_Jurkat.bedpe.gz", 99099),
    (ROOT / "gse165207" / "GSM5028232_H3K27ac_HiChIP_loops_CD34_5M.bedpe.gz", 314534),
]


def ok(path, expect):
    if not path.exists() or path.stat().st_size != expect:
        return False
    try:
        with gzip.open(path, "rb") as f:
            f.read(20)
        return True
    except Exception:
        return False


def main():
    for dest, expect in jobs:
        if ok(dest, expect):
            print("ok", dest.name, flush=True)
            continue
        if dest.exists():
            dest.unlink()
        gsm = dest.name.split("_", 1)[0]
        nnn = gsm[:7] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{nnn}/{gsm}/suppl/{dest.name}"
        print("GET", url, flush=True)
        req = urllib.request.Request(url, headers={"User-Agent": "a3-zeb1/1.0"})
        with urllib.request.urlopen(req, timeout=180) as r:
            dest.write_bytes(r.read())
        print(" wrote", dest.name, dest.stat().st_size, "expect", expect, "ok", ok(dest, expect), flush=True)
    print("REMAINING_DONE", flush=True)


if __name__ == "__main__":
    main()
