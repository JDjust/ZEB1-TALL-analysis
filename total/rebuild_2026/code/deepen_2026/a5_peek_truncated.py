import zlib
from pathlib import Path

p = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a5_gse234608\GSE234608_raw_counts.csv.gz")
data = p.read_bytes()
print("gz_bytes", len(data), flush=True)
d = zlib.decompressobj(16 + zlib.MAX_WBITS)
out = d.decompress(data)
print("partial_bytes", len(out), flush=True)
text = out.decode("utf-8", "replace")
lines = text.splitlines()
print("n_lines", len(lines), flush=True)
print("header", lines[0][:240] if lines else None, flush=True)
print("last", lines[-1][:80] if lines else None, flush=True)
for g in ("ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"):
    hit = any(l.startswith(g + ",") or l.startswith('"' + g + '"') for l in lines)
    print(g, hit, flush=True)
