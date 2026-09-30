"""GSE162280: is the ZEB2-dominant state limited to ZEB2-BCL11B fusion?

Twelve BCL11B-rearranged leukemias. No ordinary T-ALL controls are in this
series, so this is a within-BCL11B-R comparison. Partner class comes from the
author karyotype (GEO) matched to Di Giacomo et al., Blood 2021, Table 1 and
the molecular sections: every t(2;14)(q22;3;q32) produced a ZEB2-BCL11B
fusion; t(6;14)/6q25 and t(7;14) moved enhancers next to BCL11B and did not
produce a fusion. Gene-level ZEB2 counts in the fusion cases can include
fusion-transcript reads, so those cases are not used as proof of ZEB2-high.
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path

import pandas as pd

RAW = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\data\validation\gse162280\raw")
OUT = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\data\validation\gse162280")

# Author Table 1 number, UPN, GEO karyotype, paper partner, phenotype.
MANIFEST = [
    ("1", "455", "Case1", "46,XY,t(2;14)(q22;q32)", "ZEB2-BCL11B fusion", "yes", "AML"),
    ("2", "456", "Case2", "46,XY,t(2;14)(q22;q32)", "ZEB2-BCL11B fusion", "yes", "T/myeloid MPAL"),
    ("3", "457", "Case3", "46,XY,t(2;14)(q22;q32)", "ZEB2-BCL11B fusion", "yes", "AML M0"),
    ("4", "141", "Case4", "46,XX,t(2;14)(q22;q32),del(5)(q13q31)", "ZEB2-BCL11B fusion", "yes", "AML M0"),
    ("7", "459", "Case7", "46,XY.ish t(2;14)(q22;q32)", "ZEB2-BCL11B fusion", "yes", "ETP-ALL"),
    ("8", "143", "Case8", "46,XY,t(6;14)(q25;q32)", "6q25.3 ARID1B enhancer, no fusion", "no", "AML M1"),
    ("9", "145", "Case9", "46,XY,t(6;14)(q25;q32),del(11)(p13p15)", "6q25.3 ARID1B enhancer, no fusion", "no", "not labeled on the Table 1 fragment"),
    ("10", "460", "Case10", "46,XY,add(19)(q13).ish ins(6;14)(q25;q32)", "6q25.3 ARID1B enhancer, no fusion", "no", "ETP-ALL"),
    ("11", "260", "Case11", "46,XY.ish t(6;14)(q25;q32)", "6q25.3 ARID1B enhancer, no fusion", "no", "ETP-ALL"),
    ("12", "462", "Case12", "47,XX,+4/46,XX.ish t(6;14)(q25;q32)", "6q25.3 ARID1B enhancer, no fusion", "no", "ETP-ALL"),
    ("17", "149", "Case17", "46,XY,t(7;14)(q22;q32)", "7q21.2 CDK6 enhancer, no fusion", "no", "ETP-ALL"),
    ("18", "147", "Case18", "46,XY,t(7;14)(q22;q32)", "7q21.2 CDK6 enhancer, no fusion", "no", "T/myeloid MPAL"),
]

# HTSeq files are keyed by gene symbol, not Ensembl id.
GENES = {"ZEB1": "ZEB1", "ZEB2": "ZEB2", "BCL11B": "BCL11B"}


def read_counts(path: Path) -> dict[str, float]:
    out = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("__") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0].split(".")[0]
            # Case 7 has one line whose newline was lost, gluing a count to the next id.
            token = parts[1]
            digits = []
            for ch in token:
                if ch.isdigit() or ch == ".":
                    digits.append(ch)
                else:
                    break
            if not digits:
                continue
            out[gene] = float("".join(digits))
    return out


def main() -> None:
    rows = []
    for table_no, upn, case, karyotype, partner, zeb2_partner, phenotype in MANIFEST:
        hits = list(RAW.glob(f"*_{case}_HTSeqCounts.txt.gz"))
        if len(hits) != 1:
            raise SystemExit(f"{case}: found {hits}")
        counts = read_counts(hits[0])
        lib = sum(counts.values())
        rec = {
            "table_no": table_no,
            "upn": upn,
            "case": case,
            "karyotype": karyotype,
            "partner": partner,
            "zeb2_partner": zeb2_partner,
            "phenotype": phenotype,
            "library_counts": lib,
        }
        for symbol, ens in GENES.items():
            if ens not in counts:
                raise SystemExit(f"{symbol} {ens} missing in {case}")
            rec[symbol] = counts[ens]
            rec[f"{symbol}_cpm"] = 1e6 * counts[ens] / lib
            rec[f"log2_{symbol}_cpm"] = math.log2(rec[f"{symbol}_cpm"] + 1)
        rec["log2_ZEB1_over_ZEB2"] = math.log2((rec["ZEB1"] + 1) / (rec["ZEB2"] + 1))
        rec["zeb2_count_exceeds_zeb1"] = int(rec["ZEB2"] > rec["ZEB1"])
        rows.append(rec)
    df = pd.DataFrame(rows)
    # z within the 12, for a balance number. Fusion-case ZEB2 is contaminated.
    for symbol in ("ZEB1", "ZEB2"):
        x = df[f"log2_{symbol}_cpm"]
        df[f"z_{symbol}"] = (x - x.mean()) / x.std(ddof=1)
    df["balance"] = df["z_ZEB1"] - df["z_ZEB2"]
    df.to_csv(OUT / "case_manifest_and_expression.tsv", sep="\t", index=False)

    def summarize(sub: pd.DataFrame, label: str) -> dict:
        return {
            "group": label,
            "n": int(len(sub)),
            "n_ZEB2_counts_above_ZEB1": int(sub["zeb2_count_exceeds_zeb1"].sum()),
            "median_log2_ZEB1_cpm": float(sub["log2_ZEB1_cpm"].median()),
            "median_log2_ZEB2_cpm": float(sub["log2_ZEB2_cpm"].median()),
            "median_log2_ZEB1_over_ZEB2": float(sub["log2_ZEB1_over_ZEB2"].median()),
            "median_balance": float(sub["balance"].median()),
            "min_log2_ratio": float(sub["log2_ZEB1_over_ZEB2"].min()),
            "max_log2_ratio": float(sub["log2_ZEB1_over_ZEB2"].max()),
        }

    summary = pd.DataFrame([
        summarize(df, "all 12"),
        summarize(df[df.zeb2_partner == "yes"], "ZEB2-BCL11B fusion"),
        summarize(df[df.zeb2_partner == "no"], "non-ZEB2 partner"),
        summarize(df[df.partner.str.startswith("6q25")], "6q25 ARID1B enhancer"),
        summarize(df[df.partner.str.startswith("7q21")], "7q21 CDK6 enhancer"),
    ])
    summary.to_csv(OUT / "partner_group_summary.tsv", sep="\t", index=False)
    print(df[["case", "zeb2_partner", "partner", "ZEB1", "ZEB2", "log2_ZEB1_over_ZEB2", "balance"]].to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
