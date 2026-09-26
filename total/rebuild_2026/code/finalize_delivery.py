"""Copy reviewed exports into stable delivery paths and create integrity manifests."""
from pathlib import Path
import csv
import hashlib
import re
import shutil
from PIL import Image, ImageOps, ImageDraw, ImageFont

BASE = Path(__file__).resolve().parents[1]
EDITION = BASE / "edition_20260925/figures"
FIG = BASE / "figures"
MAN = BASE / "manuscript"
CONVENIENCE = BASE.parent / "manuscript"
DELIVERED = []

def copy_file(src, dst):
    if not src.is_file():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    DELIVERED.append(dst)

for group, prefix in (("main", "F"), ("supplementary", "S")):
    for n in range(1, 8):
        stem = f"{prefix}{n}"
        for ext in ("pdf", "png", "svg", "tiff"):
            copy_file(EDITION/group/f"{stem}.{ext}", FIG/group/f"{stem}.{ext}")
        for suffix in ("panel_manifest.tsv", "R_sessionInfo.txt"):
            copy_file(EDITION/group/f"{stem}_{suffix}", FIG/group/f"{stem}_{suffix}")

convenience = [
    "ZEB1_ZEB2_TALL_manuscript.md", "ZEB1_ZEB2_TALL_manuscript.docx",
    "ZEB1_ZEB2_TALL_manuscript.pdf", "Supplementary_Information.md",
    "Supplementary_Information.pdf", "ZEB1_ZEB2_TALL_Supplementary.docx",
    "Supplementary_Table_S1.pdf", "Supplementary_Table_S4.pdf",
    "Figure_Legends.md", "Supplementary_Figure_Legends_consolidated.md",
    "References.md",
]
for name in convenience:
    copy_file(MAN/name, CONVENIENCE/name)
for name in ("Table_S1_clinical_models.tsv", "Table_S2_dataset_units.tsv",
             "Table_S3_deprioritized_lines.tsv", "Table_S4_full_subtype_residuals.tsv"):
    copy_file(BASE/"data/source_data_rebuilt/supplementary_tables"/name,
              CONVENIENCE/name)
copy_file(BASE/"data/source_data_rebuilt/revision2_statistics/clinical_age_sex_wbc_sensitivity.tsv",
          CONVENIENCE/"Table_S1b_clinical_covariate_sensitivity.tsv")

review = FIG / "exports_for_review"
review.mkdir(parents=True, exist_ok=True)
for group, prefix in (("main", "F"), ("supplementary", "S")):
    thumbs=[]
    for n in range(1,8):
        with Image.open(FIG/group/f"{prefix}{n}.png") as im:
            rgb=im.convert("RGB")
            thumb=ImageOps.contain(rgb,(1500,1080),Image.Resampling.LANCZOS)
            sheet=Image.new("RGB",(1550,1135),"white")
            draw=ImageDraw.Draw(sheet)
            draw.text((22,16),f"{prefix}{n}",fill="#152b39",font=ImageFont.truetype("arial.ttf",28))
            sheet.paste(thumb,((1550-thumb.width)//2,50))
            thumbs.append(sheet)
    canvas=Image.new("RGB",(1550,1135*7),"#f5f8fa")
    for idx,thumb in enumerate(thumbs):
        canvas.paste(thumb,(0,1135*idx))
    output=review/f"{group}_contact_sheet.png"
    canvas.save(output,optimize=True)
    DELIVERED.append(output)

manifest_dir = BASE / "data/manifests"
manifest_dir.mkdir(parents=True,exist_ok=True)
files = sorted(set(DELIVERED + list((BASE/"data/source_data_rebuilt").rglob("*.tsv"))
                   + list((BASE/"code").glob("*.R")) + list((BASE/"code").glob("*.py"))
                   + [BASE/"README.md", BASE/"DELIVERY_CHECKLIST.md",
                      BASE/"logs/FINAL_QC.md", BASE/"logs/REVIEWER_AUDIT_FINAL.md",
                      BASE/"logs/VISUAL_QC_FINAL.md",
                      BASE/"data/manifests/panel_source_resolution.tsv",
                      BASE/"data/manifests/reviewer_release_candidate.zip",
                      BASE/"data/manifests/reviewer_release_candidate_contents.tsv"]
                   + [MAN/name for name in convenience]))
manifest=manifest_dir/"final_delivery_manifest.tsv"
with manifest.open("w",encoding="utf-8",newline="") as fh:
    writer=csv.writer(fh,delimiter="\t")
    writer.writerow(["path","bytes","sha256"])
    for path in files:
        with path.open("rb") as stream:
            digest=hashlib.file_digest(stream,"sha256").hexdigest()
        writer.writerow([str(path.relative_to(BASE)) if path.is_relative_to(BASE) else str(path),
                         path.stat().st_size,digest])

main_text=(MAN/"ZEB1_ZEB2_TALL_manuscript.md").read_text(encoding="utf-8")
refs=set()
for citation in re.findall(r"\[([\d,–\- ]+)\]",main_text.split("## References")[0]):
    for part in citation.split(","):
        part=part.strip()
        match=re.fullmatch(r"(\d+)\s*[–-]\s*(\d+)",part)
        if match:
            refs.update(range(int(match.group(1)),int(match.group(2))+1))
        elif part.isdigit():
            refs.add(int(part))
refs=sorted(refs)
if refs != list(range(1,61)):
    raise AssertionError(f"Reference citations incomplete: {refs}")
for group,prefix in (("main","F"),("supplementary","S")):
    for n in range(1,8):
        stem=f"{prefix}{n}"
        for ext in ("pdf","png","svg","tiff"):
            if (FIG/group/f"{stem}.{ext}").stat().st_size < 10000:
                raise AssertionError(f"Export unexpectedly small: {stem}.{ext}")
print(f"Copied {len(DELIVERED)} reviewed assets; manifest {len(files)} files; all 60 numbered references cited")
