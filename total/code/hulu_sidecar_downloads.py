#!/usr/bin/env python3
"""Sidecar: fetch PRISM despite figshare API 403, and Park Human T-cell h5ad."""
from __future__ import annotations

import json
import subprocess
import urllib.request
from pathlib import Path

RAW = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F")
PRISM = RAW / "PRISM"
PARK = RAW / "Park_HTA"
PRISM.mkdir(parents=True, exist_ok=True)
PARK.mkdir(parents=True, exist_ok=True)

UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
WGET = [
    "wget", "-c", "--tries=25", "--timeout=60", "--waitretry=8",
    "--retry-connrefused", "--no-verbose", f"--user-agent={UA}",
]


def log(m):
    print(m, flush=True)


def http_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as resp:
        return json.loads(resp.read().decode())


def wget(url, dest):
    dest = Path(dest)
    cmd = WGET + ["-O", str(dest), url]
    log("WGET " + url)
    rc = subprocess.run(cmd).returncode
    log(f"  rc={rc} size={dest.stat().st_size if dest.exists() else 0}")
    return rc == 0 and dest.exists() and dest.stat().st_size > 1000


def figshare_files(art):
    urls = [
        f"https://api.figshare.com/v2/articles/{art}/files",
        f"https://api.figshare.com/v2/articles/{art}",
        f"https://plus.figshare.com/api/articles/{art}",
    ]
    for u in urls:
        try:
            info = http_json(u)
            log(f"OK {u} type={type(info).__name__}")
            if isinstance(info, list):
                return info
            if isinstance(info, dict) and "files" in info:
                return info["files"]
        except Exception as e:
            log(f"fail {u} {e}")
    return []


def main():
    wanted = ("Matrix", "Compound", "LFC", "Repurposing", "PRISM", "prism")
    for art in (25917643, 25917655, 11383929, 24667905):
        files = figshare_files(art)
        log(f"article {art} n={len(files)}")
        for f in files:
            name = f.get("name") or ""
            log(f"  {name} {f.get('size')} {f.get('id')}")
            if not any(k.lower() in name.lower() for k in wanted + ("methotrexate", "dose-response", "compound")):
                continue
            if not any(k in name for k in ("Matrix", "Compound", "dose-response", "LFC", "Repurposing")):
                continue
            dest = PRISM / name
            url = f.get("download_url") or (f"https://figshare.com/ndownloader/files/{f['id']}" if f.get("id") else None)
            alt = f"https://ndownloader.figshare.com/files/{f['id']}" if f.get("id") else None
            if dest.exists() and f.get("size") and dest.stat().st_size >= 0.95 * f["size"]:
                log("  exists " + name)
                continue
            ok = False
            for u in (url, alt, f"https://plus.figshare.com/ndownloader/files/{f.get('id')}"):
                if u and wget(u, dest):
                    ok = True
                    break
            log(f"  got {name} {ok}")
        # whole-article zip fallback
        z = PRISM / f"article_{art}.zip"
        if not files:
            for u in (
                f"https://figshare.com/ndownloader/articles/{art}",
                f"https://plus.figshare.com/ndownloader/articles/{art}",
            ):
                if wget(u, z):
                    break

    # Park Human T cells (not TEC)
    try:
        info = http_json("https://api.cellxgene.cziscience.com/curation/v1/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910")
        for ds in info.get("datasets", []):
            title = ds.get("title") or ""
            log("dataset " + title)
            if "T cell" not in title and "Figure 1" not in title:
                continue
            for a in ds.get("assets") or []:
                url = a.get("url")
                if not url:
                    continue
                safe = "Human_T_cells.h5ad" if "T cell" in title else "Park_Figure1.h5ad"
                dest = PARK / safe
                if dest.exists() and dest.stat().st_size > 10**7:
                    log("exists " + safe)
                    continue
                wget(url, dest)
    except Exception as e:
        log("cellxgene " + str(e))
    log("SIDECAR_DONE")


if __name__ == "__main__":
    main()
