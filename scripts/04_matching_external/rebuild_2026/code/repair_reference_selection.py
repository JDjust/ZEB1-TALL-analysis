"""Keep original-cohort references for the historical meta-analysis (editorial fix)."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[3]
path = root / "total/rebuild_2026/manuscript/ZEB1_ZEB2_TALL_manuscript_polished.md"
archive = root / "_archive/old_versions/rebuild_2026_manuscript_outputs/ZEB1_ZEB2_TALL_manuscript_prepolish.md"
text = path.read_text(encoding="utf-8")
head, refs_text = text.split("## References\n\n", 1)
current = {int(m.group(1)): m.group(2) for m in
           re.finditer(r"(?m)^(\d+)\. (.+)$", refs_text)}
original = {int(m.group(1)): m.group(2) for m in
            re.finditer(r"(?m)^(\d+)\. (.+)$", archive.read_text(encoding="utf-8").split("## References", 1)[1])}
old_drop = {10, 21, 23, 25, 26, 27, 29, 31, 46, 60}
old_kept = [n for n in range(1, 61) if n not in old_drop]
assert len(current) == len(old_kept) == 50
new_drop = {10, 25, 26, 29, 31, 46, 53, 55, 57, 60}
new_kept = [n for n in range(1, 61) if n not in new_drop]
old_from_current = {i: old for i, old in enumerate(old_kept, 1)}
new_from_old = {old: i for i, old in enumerate(new_kept, 1)}

def convert(match):
    body = match.group(1)
    if not re.fullmatch(r"[\d\s,;–-]+", body):
        return match.group(0)
    ids = []
    for part in re.split(r"[,;]", body):
        nums = [int(x) for x in re.findall(r"\d+", part)]
        ids.extend(range(nums[0], nums[1] + 1) if len(nums) == 2 else nums)
    originals = {old_from_current[n] for n in ids}
    values = sorted(new_from_old[n] for n in originals if n in new_from_old)
    assert values
    return "[" + ",".join(map(str, values)) + "]"

head = re.sub(r"\[([^\]]+)\]", convert, head)
# The first editorial pass lost source-cohort citations in this one sentence.
head = head.replace(
    "lower ZEB2 in TAL-class disease (Fig. 5F; Supplementary Fig. S7) [20].",
    "lower ZEB2 in TAL-class disease (Fig. 5F; Supplementary Fig. S7) [21,23,27].",
)
head = head.replace(
    "clinical relevance of ZEB balance is mediated through molecular subtype composition",
    "clinical relevance of ZEB balance is largely explained by molecular subtype composition",
)
refs = {}
for old in new_kept:
    if old in old_kept:
        refs[old] = current[old_kept.index(old) + 1]
    else:
        refs[old] = original[old]
text = head + "## References\n\n" + "\n\n".join(
    f"{new_from_old[old]}. {refs[old]}" for old in new_kept) + "\n"
path.write_text(text, encoding="utf-8")
print("Repaired historical source citations and retained 50 references.")
