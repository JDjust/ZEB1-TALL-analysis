#!/usr/bin/env python3
"""Download Park/HTA, PRISM, GDSC on hulu; run decoupler + TARGET Cox.

Writes compact tables for Figure 4F / 6F / 7F (and a Park table for S2).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import traceback
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
OUT = ROOT / "total_fig467F"
RAW = DATA / "external_fig467F"
PYDEPS = OUT / "pydeps"
LOG = OUT / "logs"
for d in (OUT / "tables", RAW / "PRISM", RAW / "GDSC", RAW / "Park_HTA", PYDEPS, LOG):
    d.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(PYDEPS))
PY = sys.executable
WGET = [
    "wget", "-c", "--tries=25", "--timeout=60", "--waitretry=8",
    "--retry-connrefused", "--no-verbose",
]
UA = {"User-Agent": "ZEB1-TALL-analysis/1.0"}


def log(msg):
    print(msg, flush=True)


def write_tsv(df, name):
    path = OUT / "tables" / name
    df.to_csv(path, sep="\t", index=False)
    log(f"wrote {path} nrow={len(df)}")
    return path


def wget_to(url, dest, timeout_s=0):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = WGET + ["-O", str(dest), url]
    log("WGET " + url)
    p = subprocess.run(cmd, timeout=timeout_s or None)
    ok = p.returncode == 0 and dest.exists() and dest.stat().st_size > 1000
    log(f"  rc={p.returncode} size={dest.stat().st_size if dest.exists() else 0} ok={ok}")
    return ok


def start_wget(url, dest):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    lg = dest.with_suffix(dest.suffix + ".wget.log")
    cmd = WGET + ["-O", str(dest), url]
    log("WGET-BG " + url)
    fh = open(lg, "ab")
    return subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT)


def http_json(url, timeout=90):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def pip_install(*pkgs):
    cmd = [PY, "-m", "pip", "install", "--target", str(PYDEPS), "--upgrade", *pkgs]
    log("PIP " + " ".join(pkgs))
    p = subprocess.run(cmd)
    log(f"  pip rc={p.returncode}")
    return p.returncode == 0


def fdr_bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    mask = np.isfinite(p)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranked = np.empty(n, dtype=float)
    prev = 1.0
    for i, k in enumerate(order[::-1]):
        rank = n - i
        val = min(prev, pv[k] * n / rank)
        ranked[k] = val
        prev = val
    out[mask] = ranked
    return out


def spearman_ci(x, y, n_boot=1500, seed=1):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 4:
        return dict(n=n, rho=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        r, _ = stats.spearmanr(x[i], y[i])
        if np.isfinite(r):
            boots.append(r)
    lo, hi = np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    return dict(n=n, rho=float(rho), p=float(p), ci_lo=float(lo), ci_hi=float(hi))


def coxph_fit(time, event, X, maxiter=40):
    """Breslow Cox PH. X is n x p (no intercept). Returns beta, se, n, n_event."""
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    m = np.isfinite(t) & np.isfinite(e) & np.all(np.isfinite(X), axis=1) & (t > 0)
    t, e, X = t[m], e[m], X[m]
    n, p = X.shape
    n_event = int(e.sum())
    if n < 20 or n_event < 8:
        return dict(n=int(n), n_event=n_event, beta=np.full(p, np.nan),
                    se=np.full(p, np.nan))
    order = np.argsort(t, kind="mergesort")
    t, e, X = t[order], e[order], X[order]
    beta = np.zeros(p)
    hess = np.eye(p)
    for _ in range(maxiter):
        eta = X @ beta
        eta = np.clip(eta, -20, 20)
        risk = np.exp(eta)
        R = np.cumsum(risk[::-1])[::-1]
        SX = np.cumsum((X * risk[:, None])[::-1], axis=0)[::-1]
        grad = np.zeros(p)
        H = np.zeros((p, p))
        i = 0
        while i < n:
            j = i
            while j < n and t[j] == t[i]:
                j += 1
            d = e[i:j].sum()
            if d > 0 and R[i] > 0:
                mean = SX[i] / R[i]
                # variance of covariates on the risk set (outer product approximation)
                wX = X[i:] * risk[i:, None]
                second = (wX.T @ X[i:]) / R[i]
                var = second - np.outer(mean, mean)
                grad += X[i:j][e[i:j] == 1].sum(axis=0) - d * mean
                H += d * var
            i = j
        try:
            step = np.linalg.solve(H + 1e-8 * np.eye(p), grad)
        except np.linalg.LinAlgError:
            break
        beta = beta + step
        hess = H
        if np.max(np.abs(step)) < 1e-8:
            break
    try:
        cov = np.linalg.inv(hess + 1e-8 * np.eye(p))
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)
    return dict(n=int(n), n_event=n_event, beta=beta, se=se)


def cox_row(endpoint, model, time, event, names, cols, covariates):
    fit = coxph_fit(time, event, cols)
    rows = []
    for i, name in enumerate(names):
        b = fit["beta"][i]
        se = fit["se"][i]
        hr = np.exp(b) if np.isfinite(b) else np.nan
        z = b / se if (np.isfinite(b) and np.isfinite(se) and se > 0) else np.nan
        p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
        rows.append({
            "endpoint": endpoint, "model": model, "term": name,
            "n": fit["n"], "n_event": fit["n_event"],
            "hr": float(hr) if np.isfinite(hr) else np.nan,
            "ci_lo": float(np.exp(b - 1.96 * se)) if np.isfinite(hr) and np.isfinite(se) else np.nan,
            "ci_hi": float(np.exp(b + 1.96 * se)) if np.isfinite(hr) and np.isfinite(se) else np.nan,
            "p": p, "covariates": covariates,
        })
    return rows


def zscore(x):
    x = np.asarray(x, float)
    s = np.nanstd(x)
    if not np.isfinite(s) or s == 0:
        return x * np.nan
    return (x - np.nanmean(x)) / s


def pick_col(cols, needles):
    cols = list(cols)
    low = [str(c).replace("\n", " ").replace("\r", " ").strip() for c in cols]
    for nd in needles:
        ndl = nd.lower()
        for c, lc in zip(cols, low):
            if lc.lower() == ndl:
                return c
    for nd in needles:
        ndl = nd.lower()
        for c, lc in zip(cols, low):
            if ndl in lc.lower():
                return c
    return None


def gene_token(col):
    return str(col).split(" ")[0].split("(")[0]


# ---------------------------------------------------------------------------
# downloads
# ---------------------------------------------------------------------------
def download_prism():
    dest_dir = RAW / "PRISM"
    wanted = [
        "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
        "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv",
        "Repurposing_Public_24Q2_Extended_Primary_Quality_Control.csv",
    ]
    arts = [25917643, 25917655, 24667905]
    got = {}
    for art in arts:
        try:
            files = http_json(f"https://api.figshare.com/v2/articles/{art}/files")
        except Exception as e:
            log(f"figshare {art} fail {e}")
            continue
        log(f"figshare {art} n_files={len(files)}")
        for f in files:
            name = f.get("name") or ""
            log(f"  file {name} size={f.get('size')}")
            if name not in wanted and not any(k in name for k in ("Prism", "PRISM", "Repurposing", "methotrexate", "Methotrexate")):
                if "Compound" not in name and "Matrix" not in name and "LFC" not in name:
                    continue
            dest = dest_dir / name
            url = f.get("download_url")
            if dest.exists() and f.get("size") and dest.stat().st_size == f["size"]:
                log(f"  skip exists {name}")
                got[name] = dest
                continue
            if url and wget_to(url, dest):
                got[name] = dest
    # 19Q4 secondary IC50 as extra
    try:
        files = http_json("https://api.figshare.com/v2/articles/11383929/files")
        for f in files:
            name = f.get("name") or ""
            if "dose-response-curve-parameters" in name or "secondary-screen-dose" in name:
                dest = dest_dir / name
                if (not dest.exists() or dest.stat().st_size < 1000) and f.get("download_url"):
                    wget_to(f["download_url"], dest)
                if dest.exists():
                    got[name] = dest
    except Exception as e:
        log(f"prism 19q4 fail {e}")
    return got


def download_gdsc():
    dest_dir = RAW / "GDSC"
    urls = {
        "GDSC2_fitted_dose_response_24Jul22.csv":
            "https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_24Jul22.csv",
        "GDSC1_fitted_dose_response_24Jul22.csv":
            "https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_24Jul22.csv",
        "Cell_Lines_Details.xlsx":
            "https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx",
        "GDSC2_fitted_dose_response_27Oct23.xlsx":
            "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx",
    }
    got = {}
    for name, url in urls.items():
        dest = dest_dir / name
        if dest.exists() and dest.stat().st_size > 10000:
            log(f"gdsc exists {name} {dest.stat().st_size}")
            got[name] = dest
            continue
        if wget_to(url, dest):
            got[name] = dest
    return got


def start_park_download():
    dest_dir = RAW / "Park_HTA"
    procs = []
    # CELLxGENE collections
    coll_ids = [
        "de13e3e2-23b6-40ed-a413-e9e12d7d3910",
        "de13e19e-6e0c-4f43-b9aa-30daab96ca6c",
    ]
    h5ads = []
    for cid in coll_ids:
        try:
            info = http_json(f"https://api.cellxgene.cziscience.com/curation/v1/collections/{cid}")
        except Exception as e:
            log(f"cellxgene {cid} {e}")
            continue
        log(f"cellxgene collection {cid} keys={list(info)[:12]}")
        for ds in info.get("datasets", []) or []:
            title = ds.get("title") or ds.get("name") or ""
            assets = ds.get("assets") or ds.get("dataset_assets") or []
            log(f"  dataset {title} n_assets={len(assets)}")
            for a in assets:
                ftype = (a.get("filetype") or a.get("file_type") or "").lower()
                url = a.get("url") or a.get("s3_uri") or ""
                if "h5ad" in ftype or str(url).endswith(".h5ad"):
                    h5ads.append((title, url, a.get("filesize") or a.get("filesize_raw") or 0))
    h5ads.sort(key=lambda x: x[2] or 10**18)
    for title, url, sz in h5ads[:3]:
        if not url:
            continue
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in title)[:80] + ".h5ad"
        dest = dest_dir / safe
        if dest.exists() and dest.stat().st_size > 10**6:
            log(f"park h5ad exists {dest}")
            continue
        log(f"park cellxgene {title} size~{sz}")
        procs.append(start_wget(url, dest))
        break
    zdest = dest_dir / "thymus_annotated_matrix_files.zip"
    if not (zdest.exists() and zdest.stat().st_size > 10**8):
        procs.append(start_wget(
            "https://zenodo.org/records/5500511/files/thymus_annotated_matrix_files.zip?download=1",
            zdest,
        ))
    else:
        log("park zip exists")
    mdest = dest_dir / "sample_metadata_fix.xlsx"
    if not mdest.exists():
        procs.append(start_wget(
            "https://zenodo.org/records/5500511/files/sample_metadata_fix.xlsx?download=1",
            mdest,
        ))
    return procs


def wait_procs(procs, tag, max_s=10800):
    t0 = time.time()
    for p in procs:
        remain = max(1, int(max_s - (time.time() - t0)))
        try:
            rc = p.wait(timeout=remain)
            log(f"{tag} proc pid={p.pid} rc={rc}")
        except subprocess.TimeoutExpired:
            log(f"{tag} timeout, killing {p.pid}")
            p.kill()


# ---------------------------------------------------------------------------
# TARGET survival
# ---------------------------------------------------------------------------
def analyze_target_cox():
    clin_p = ROOT / "module1/tables/target_clinical_combined.tsv"
    samp_p = ROOT / "module4/processed/target_samples.tsv"
    tpm_p = ROOT / "module4/processed/target_log2tpm_filtered.tsv.gz"
    clin = pd.read_csv(clin_p, sep="\t")
    samp = pd.read_csv(samp_p, sep="\t")
    log("clinical cols: " + " | ".join(map(str, clin.columns)))
    log("samples cols: " + " | ".join(map(str, samp.columns)))
    log(f"clinical {clin.shape} samples {samp.shape}")

    usi_c = pick_col(clin.columns, ["TARGET USI", "usi", "case_submitter_id", "patient"])
    os_t = pick_col(clin.columns, [
        "Overall Survival Time in Days", "overall survival time", "OS.time", "OS time",
    ])
    os_e = pick_col(clin.columns, ["Vital Status", "OS", "overall survival indicator", "vital"])
    efs_t = pick_col(clin.columns, [
        "Event Free Survival Time in Days", "event free survival time", "EFS.time", "EFS time",
    ])
    efs_e = pick_col(clin.columns, ["First Event", "EFS", "event free"])
    age_c = pick_col(clin.columns, ["Age at Diagnosis in Days", "age at diagnosis", "Age"])
    wbc_c = pick_col(clin.columns, ["WBC at Diagnosis", "WBC"])
    log(f"picked usi={usi_c} os_t={os_t} os_e={os_e} efs_t={efs_t} efs_e={efs_e} age={age_c}")

    samp = samp.copy()
    sid = pick_col(samp.columns, ["sample", "barcode", "submitter", "USI", "patient", "id"])
    if sid is None:
        sid = samp.columns[0]
    samp["_usi6"] = samp[sid].astype(str).str.replace("-", "", regex=False).str.extract(r"([A-Z]{5,6})$")[0]
    if usi_c:
        clin["_usi6"] = clin[usi_c].astype(str).str.replace("-", "", regex=False).str.extract(r"([A-Z]{5,6})$")[0]
    else:
        clin["_usi6"] = clin.iloc[:, 0].astype(str).str.extract(r"([A-Z]{5,6})$")[0]

    # ZEB1
    zeb = None
    for c in samp.columns:
        if gene_token(c).upper() == "ZEB1":
            zeb = samp[c]
            break
    if zeb is None and tpm_p.exists():
        tpm = pd.read_csv(tpm_p, sep="\t", index_col=0)
        log(f"tpm {tpm.shape}")
        g = "ZEB1" if "ZEB1" in tpm.index else None
        if g is None:
            hits = [i for i in tpm.index.astype(str) if str(i).startswith("ZEB1")]
            g = hits[0] if hits else None
        if g:
            zmap = tpm.loc[g]
            key = pick_col(samp.columns, ["sample_id", "sample", "barcode", "id"])
            zeb = samp[key].astype(str).map(zmap)
    if zeb is None:
        raise RuntimeError("no ZEB1 in TARGET samples/tpm")
    samp["ZEB1"] = pd.to_numeric(zeb, errors="coerce")

    m = samp.merge(clin, on="_usi6", how="left", suffixes=("", "_clin"))
    log(f"merged {m.shape} ZEB1 non-na={m['ZEB1'].notna().sum()}")

    def vital_event(s):
        x = s.astype(str).str.lower()
        ev = np.where(x.str.contains("dead|deceased|died|1|true"), 1,
                      np.where(x.str.contains("alive|0|false|censored"), 0, np.nan))
        return pd.to_numeric(ev, errors="coerce")

    def efs_event(s):
        x = s.astype(str).str.lower()
        # Censored / none = 0; relapse, death, induction failure = 1
        ev = np.where(x.str.contains("censor|none|no event|alive"), 0,
                      np.where(x.str.contains("relapse|death|dead|fail|progress|event"), 1, np.nan))
        # if already 0/1
        num = pd.to_numeric(s, errors="coerce")
        ev = np.where(num.isin([0, 1]), num, ev)
        return pd.to_numeric(ev, errors="coerce")

    if os_t:
        m["os_time"] = pd.to_numeric(m[os_t], errors="coerce")
    else:
        m["os_time"] = np.nan
    m["os_event"] = vital_event(m[os_e]) if os_e else np.nan
    if efs_t:
        m["efs_time"] = pd.to_numeric(m[efs_t], errors="coerce")
    else:
        m["efs_time"] = np.nan
    m["efs_event"] = efs_event(m[efs_e]) if efs_e else np.nan
    if age_c:
        m["age_days"] = pd.to_numeric(m[age_c], errors="coerce")
    else:
        m["age_days"] = np.nan

    sub_c = pick_col(samp.columns, ["level1", "subtype", "Level1", "group"])
    if sub_c:
        m["subtype"] = m[sub_c].astype(str)
    else:
        m["subtype"] = "NA"

    keep_cols = ["_usi6", sid, "ZEB1", "os_time", "os_event", "efs_time", "efs_event",
                 "age_days", "subtype"]
    slim = m[keep_cols].drop_duplicates("_usi6")
    write_tsv(slim, "F7_target_surv_merged.tsv")

    z = zscore(slim["ZEB1"])
    agez = zscore(slim["age_days"])
    dummies = pd.get_dummies(slim["subtype"].fillna("NA"), drop_first=True).astype(float)
    rows = []
    for ep, tt, ee in (("OS", "os_time", "os_event"), ("EFS", "efs_time", "efs_event")):
        rows += cox_row(ep, "ZEB1_unadj", slim[tt], slim[ee], ["ZEB1_z"], z, "ZEB1 (per SD)")
        X2 = np.column_stack([z, agez])
        rows += cox_row(ep, "ZEB1_age", slim[tt], slim[ee], ["ZEB1_z", "age_z"], X2,
                        "ZEB1 (per SD) + age")
        if dummies.shape[1] >= 1:
            X3 = np.column_stack([z, dummies.to_numpy()])
            names = ["ZEB1_z"] + [f"sub_{c}" for c in dummies.columns]
            rows += cox_row(ep, "ZEB1_subtype", slim[tt], slim[ee], names, X3,
                            "ZEB1 (per SD) + subtype")
    out = pd.DataFrame(rows)
    write_tsv(out, "F7_target_cox.tsv")
    return out


# ---------------------------------------------------------------------------
# decoupler on GSE110635
# ---------------------------------------------------------------------------
def analyze_decoupler():
    pip_install("decoupler")
    import decoupler as dc
    log("decoupler " + str(getattr(dc, "__version__", "?")))
    counts = pd.read_csv(
        ROOT / "module9/processed/gse110635_symbol_counts.tsv.gz", sep="\t", index_col=0,
    )
    cd = pd.read_csv(ROOT / "module9/processed/gse110635_symbol_coldata.tsv", sep="\t")
    log(f"counts {counts.shape} coldata {cd.shape} groups={cd['group'].value_counts().to_dict() if 'group' in cd.columns else cd.head()}")
    # samples as columns
    lib = counts.sum(axis=0).replace(0, np.nan)
    logcpm = np.log1p(counts.divide(lib, axis=1) * 1e6)
    mat = logcpm.T
    mat.index = mat.index.astype(str)
    cd = cd.copy()
    key = pick_col(cd.columns, ["gsm", "sample", "id"])
    cd["_id"] = cd[key].astype(str)
    cd = cd.set_index("_id")
    common = [i for i in mat.index if i in cd.index]
    if len(common) < 4:
        # try title / file stem
        log("id overlap low, trying column names vs gsm")
        common = [c for c in mat.index if c in set(cd.index)]
    mat = mat.loc[common]
    cd = cd.loc[common]
    net = None
    for getter in (
        lambda: dc.get_collectri(organism="human", split_complexes=False),
        lambda: dc.op.collectri(organism="human"),
        lambda: dc.translate_net(dc.get_collectri()),
    ):
        try:
            net = getter()
            log(f"collectri {getattr(net, 'shape', type(net))}")
            break
        except Exception as e:
            log(f"collectri getter fail {e}")
    if net is None:
        raise RuntimeError("could not load CollecTRI")
    # normalize net columns
    if isinstance(net, pd.DataFrame):
        rename = {}
        for c in net.columns:
            cl = c.lower()
            if cl in ("source", "tf", "tf.symbol"):
                rename[c] = "source"
            elif cl in ("target", "target_genesymbol", "gene"):
                rename[c] = "target"
            elif cl in ("weight", "mor", "likelihood"):
                rename[c] = "weight"
        net = net.rename(columns=rename)
        if "weight" not in net.columns:
            net["weight"] = 1.0
    acts = None
    for runner in (
        lambda: dc.run_ulm(mat=mat, net=net, source="source", target="target",
                           weight="weight", min_n=5),
        lambda: dc.mt.ulm(mat, net),
    ):
        try:
            out = runner()
            if isinstance(out, tuple):
                acts = out[0]
            else:
                acts = out
            log(f"ulm acts {getattr(acts, 'shape', type(acts))}")
            break
        except Exception as e:
            log(f"ulm runner fail {e}")
            traceback.print_exc()
    if acts is None:
        raise RuntimeError("ULM failed")
    if not isinstance(acts, pd.DataFrame):
        acts = pd.DataFrame(acts)
    # samples x TFs
    if acts.shape[0] != mat.shape[0] and acts.shape[1] == mat.shape[0]:
        acts = acts.T
    acts.index = acts.index.astype(str)
    grp = cd["group"].astype(str) if "group" in cd.columns else cd.iloc[:, 0].astype(str)
    ctrl = grp.str.lower().str.contains("control|scrambled|ntc|sint|neg")
    kd = grp.str.lower().str.contains("sitlx|kd|knock")
    if int(kd.sum()) == 0:
        kd = ~ctrl
    log(f"ctrl={int(ctrl.sum())} kd={int(kd.sum())} groups={grp.value_counts().to_dict()}")
    rows = []
    ctrl_ids = [i for i in cd.index[ctrl] if i in acts.index]
    kd_ids = [i for i in cd.index[kd] if i in acts.index]
    for tf in acts.columns:
        a = pd.to_numeric(acts.reindex(ctrl_ids)[tf], errors="coerce").dropna()
        b = pd.to_numeric(acts.reindex(kd_ids)[tf], errors="coerce").dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        delta = float(b.mean() - a.mean())
        t, p = stats.ttest_ind(b, a, equal_var=False)
        rows.append({
            "tf": tf, "mean_ctrl": float(a.mean()), "mean_kd": float(b.mean()),
            "delta_kd_minus_ctrl": delta, "t": float(t), "p": float(p),
            "n_ctrl": int(len(a)), "n_kd": int(len(b)),
        })
    out = pd.DataFrame(rows)
    out["fdr"] = fdr_bh(out["p"])
    out = out.sort_values("p")
    write_tsv(out, "F6_decoupler_ulm_all.tsv")
    focus = ["TLX1", "ZEB1", "GATA3", "TCF7", "BCL11B", "MYC", "MEF2C",
             "RUNX1", "TAL1", "LMO2", "LYL1", "LEF1", "NOTCH1", "MYB", "HOXA9"]
    sub = out[out["tf"].isin(focus)].copy()
    write_tsv(sub, "F6_decoupler_ulm_focus.tsv")
    return out


# ---------------------------------------------------------------------------
# MTX PRISM / GDSC
# ---------------------------------------------------------------------------
def load_depmap_tall():
    model = pd.read_csv(DATA / "DepMap_26Q1/Model.csv")
    model["lineage"] = np.where(
        model["OncotreePrimaryDisease"].eq("T-Lymphoblastic Leukemia/Lymphoma"), "T-ALL",
        np.where(model["OncotreePrimaryDisease"].eq("B-Lymphoblastic Leukemia/Lymphoma"), "B-ALL",
                 np.where(model["OncotreeLineage"].eq("Lymphoid"), "other_lymphoid",
                          np.where(model["OncotreeLineage"].eq("Myeloid"), "myeloid", "other"))),
    )
    expr_p = DATA / "DepMap_26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    # only need ZEB1 + ModelID
    header = pd.read_csv(expr_p, nrows=0)
    cols = list(header.columns)
    zeb_col = None
    for c in cols:
        if gene_token(c) == "ZEB1":
            zeb_col = c
            break
    idcol = cols[0]
    use = [idcol, zeb_col] if zeb_col else [idcol]
    expr = pd.read_csv(expr_p, usecols=use)
    expr = expr.rename(columns={idcol: "ModelID", zeb_col: "ZEB1"} if zeb_col else {idcol: "ModelID"})
    m = model.merge(expr, on="ModelID", how="left")
    log("depmap lineage " + str(m.lineage.value_counts().to_dict()))
    log(f"T-ALL n={int(m.lineage.eq('T-ALL').sum())} with ZEB1={int(m.loc[m.lineage.eq('T-ALL'), 'ZEB1'].notna().sum())}")
    return m


def analyze_mtx(dep):
    points = []
    stats_rows = []

    def add_points(src, metric, df, line_col, val_col, id_col=None):
        d = df.copy()
        d["_key"] = d[line_col].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
        dep2 = dep.copy()
        dep2["_key"] = dep2["StrippedCellLineName"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
        dep2["_key2"] = dep2["CellLineName"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
        merged = d.merge(dep2[["ModelID", "CellLineName", "lineage", "ZEB1", "_key"]], on="_key", how="left")
        miss = merged["ModelID"].isna()
        if miss.any():
            extra = d.merge(dep2[["ModelID", "CellLineName", "lineage", "ZEB1", "_key2"]].rename(columns={"_key2": "_key"}),
                            on="_key", how="left")
            merged.loc[miss, ["ModelID", "CellLineName", "lineage", "ZEB1"]] = extra.loc[miss, ["ModelID", "CellLineName", "lineage", "ZEB1"]].values
        if id_col and id_col in d.columns:
            by_id = d.merge(dep[["ModelID", "CellLineName", "lineage", "ZEB1"]],
                            left_on=id_col, right_on="ModelID", how="left")
            still = merged["ZEB1"].isna()
            if still.any() and "ZEB1" in by_id.columns:
                merged.loc[still, ["ModelID", "CellLineName", "lineage", "ZEB1"]] = by_id.loc[still, ["ModelID", "CellLineName", "lineage", "ZEB1"]].values
        for _, r in merged.iterrows():
            points.append({
                "source": src, "metric": metric, "drug": "Methotrexate",
                "cell_line": r.get("CellLineName") if pd.notna(r.get("CellLineName")) else r[line_col],
                "model_id": r.get("ModelID"), "lineage": r.get("lineage"),
                "ZEB1": r.get("ZEB1"), "value": r[val_col],
            })
        mp = pd.DataFrame(points)
        mp = mp[mp["source"].eq(src)]
        for lin_lab, sub in (("T-ALL", mp[mp.lineage.eq("T-ALL")]),
                             ("ALL_TB", mp[mp.lineage.isin(["T-ALL", "B-ALL"])]),
                             ("all_matched", mp.dropna(subset=["ZEB1"]))):
            sp = spearman_ci(sub["ZEB1"], sub["value"])
            stats_rows.append({
                "source": src, "metric": metric, "subset": lin_lab,
                **sp,
            })

    # GDSC
    for fname, src in (
        ("GDSC2_fitted_dose_response_24Jul22.csv", "GDSC2"),
        ("GDSC1_fitted_dose_response_24Jul22.csv", "GDSC1"),
    ):
        p = RAW / "GDSC" / fname
        xlsx = RAW / "GDSC" / "GDSC2_fitted_dose_response_27Oct23.xlsx"
        df = None
        if p.exists():
            df = pd.read_csv(p, low_memory=False)
        elif src == "GDSC2" and xlsx.exists():
            df = pd.read_excel(xlsx)
        if df is None:
            log(f"missing {src}")
            continue
        dcol = pick_col(df.columns, ["DRUG_NAME", "drug_name"])
        mtx = df[df[dcol].astype(str).str.contains("methotrexate", case=False, na=False)].copy()
        log(f"{src} mtx rows={len(mtx)} cols={list(df.columns)[:12]}")
        if mtx.empty:
            continue
        line = pick_col(mtx.columns, ["CELL_LINE_NAME", "cell_line_name"])
        val = pick_col(mtx.columns, ["LN_IC50", "ln_ic50", "IC50"])
        mtx = mtx.groupby(line, as_index=False)[val].mean()
        add_points(src, "LN_IC50", mtx, line, val)
        write_tsv(mtx.head(5000), f"F7_{src}_mtx_raw.tsv")

    # PRISM primary matrix
    prism_files = list((RAW / "PRISM").glob("*"))
    log("prism files " + str([f.name for f in prism_files]))
    matrix = None
    compounds = None
    for f in prism_files:
        n = f.name.lower()
        if "matrix" in n and f.suffix == ".csv":
            matrix = f
        if "compound" in n and f.suffix == ".csv":
            compounds = f
    if matrix is not None:
        log(f"reading prism {matrix}")
        # compound list to find methotrexate row id
        mtx_ids = ["methotrexate"]
        if compounds is not None:
            cpd = pd.read_csv(compounds, low_memory=False)
            log("prism compound cols " + str(list(cpd.columns)))
            namec = pick_col(cpd.columns, ["name", "drug_name", "compound", "IDs", "broad_id"])
            idc = pick_col(cpd.columns, ["broad_id", "column_name", "IDs", "id"])
            hit = cpd[cpd.apply(lambda r: r.astype(str).str.contains("methotrexate", case=False).any(), axis=1)]
            log(f"prism mtx compounds {len(hit)}")
            write_tsv(hit, "F7_PRISM_mtx_compounds.tsv")
            mtx_ids = hit.astype(str).stack().unique().tolist()
        # read matrix header only first
        hdr = pd.read_csv(matrix, nrows=0)
        log(f"prism matrix ncol={len(hdr.columns)} first={list(hdr.columns)[:8]}")
        idx_col = hdr.columns[0]
        # scan for methotrexate in first column via chunks
        chunks = []
        for ch in pd.read_csv(matrix, chunksize=400):
            key = ch.iloc[:, 0].astype(str)
            mask = key.str.contains("methotrexate", case=False, na=False)
            if not mask.any() and mtx_ids:
                mask = key.isin(set(map(str, mtx_ids))) | key.str.lower().isin({str(x).lower() for x in mtx_ids})
            if mask.any():
                chunks.append(ch.loc[mask])
        if chunks:
            mtxm = pd.concat(chunks, ignore_index=True)
            log(f"prism mtx matrix rows={len(mtxm)}")
            write_tsv(mtxm.iloc[:, :5], "F7_PRISM_mtx_matrix_head.tsv")
            # wide: row=compound, cols=ACH ids
            id_like = [c for c in mtxm.columns if str(c).startswith("ACH-")]
            if not id_like:
                # maybe rows are cell lines
                id_like = [c for c in mtxm.columns if str(c).startswith("BRD-")]
            if id_like:
                row = mtxm.iloc[0]
                long = pd.DataFrame({"ModelID": id_like, "value": pd.to_numeric(row[id_like], errors="coerce")})
                long = long.merge(dep[["ModelID", "CellLineName", "lineage", "ZEB1"]], on="ModelID", how="left")
                for _, r in long.iterrows():
                    points.append({
                        "source": "PRISM24Q2", "metric": "LFC_2.5uM", "drug": "Methotrexate",
                        "cell_line": r["CellLineName"], "model_id": r["ModelID"],
                        "lineage": r["lineage"], "ZEB1": r["ZEB1"], "value": r["value"],
                    })
                mp = pd.DataFrame([p for p in points if p["source"] == "PRISM24Q2"])
                for lin_lab, sub in (("T-ALL", mp[mp.lineage.eq("T-ALL")]),
                                     ("ALL_TB", mp[mp.lineage.isin(["T-ALL", "B-ALL"])]),
                                     ("all_matched", mp.dropna(subset=["ZEB1"]))):
                    sp = spearman_ci(sub["ZEB1"], sub["value"])
                    stats_rows.append({"source": "PRISM24Q2", "metric": "LFC_2.5uM", "subset": lin_lab, **sp})
            else:
                log("prism matrix layout not ACH-wide; trying long format")
                # long format with cell line and LFC columns
                line = pick_col(mtxm.columns, ["ccle_name", "cell_line", "depmap_id", "ModelID"])
                val = pick_col(mtxm.columns, ["LFC", "lfc", "ic50", "auc", "log2"])
                if line and val:
                    add_points("PRISM24Q2", "LFC", mtxm, line, val,
                               id_col=pick_col(mtxm.columns, ["depmap_id", "ModelID", "ACH"]))

    # 19Q4 dose-response parameters
    for f in prism_files:
        if "dose-response-curve-parameters" in f.name or "secondary-screen-dose" in f.name:
            dr = pd.read_csv(f, low_memory=False)
            log(f"prism DR {f.name} {dr.shape} {list(dr.columns)[:15]}")
            namec = pick_col(dr.columns, ["name", "drug_name", "compound"])
            mtx = dr[dr[namec].astype(str).str.contains("methotrexate", case=False, na=False)] if namec else dr
            if mtx.empty:
                continue
            line = pick_col(mtx.columns, ["ccle_name", "cell_line", "depmap_id"])
            val = pick_col(mtx.columns, ["ic50", "auc", "ec50"])
            if line and val:
                add_points("PRISM19Q4_DR", str(val), mtx, line, val)

    pts = pd.DataFrame(points)
    st = pd.DataFrame(stats_rows)
    if len(pts):
        write_tsv(pts, "F7_mtx_points.tsv")
    if len(st):
        write_tsv(st, "F7_mtx_stats.tsv")
    return pts, st


# ---------------------------------------------------------------------------
# Park / HTA
# ---------------------------------------------------------------------------
def analyze_park():
    dest_dir = RAW / "Park_HTA"
    genes = ["ZEB1", "ZEB2", "LMO2", "TCF7", "GATA3", "BCL11B", "CD3E", "CD1A", "DNTT", "CD34"]
    h5ads = list(dest_dir.glob("*.h5ad"))
    zips = list(dest_dir.glob("*.zip"))
    log(f"park h5ads={ [p.name for p in h5ads] } zips={[p.name for p in zips]}")
    adata_path = None
    if zips:
        z = zips[0]
        if z.stat().st_size > 10**7:
            import zipfile
            with zipfile.ZipFile(z) as zf:
                names = zf.namelist()
                log("zip n=" + str(len(names)))
                for n in names[:40]:
                    log("  " + n)
                csvs = [n for n in names if n.lower().endswith(".csv") and "meta" in n.lower()]
                h5 = [n for n in names if n.lower().endswith(".h5ad")]
                # prefer annotated T / fig1
                pref = [n for n in h5 if any(k in n.lower() for k in ("fig1", "tcell", "t_cell", "human", "anno"))]
                pick = (pref or h5)[:1]
                extract_dir = dest_dir / "unzipped"
                extract_dir.mkdir(exist_ok=True)
                for n in pick + csvs[:3]:
                    outp = extract_dir / Path(n).name
                    if not outp.exists():
                        log("extract " + n)
                        zf.extract(n, extract_dir)
                        # flatten if nested
                    nested = extract_dir / n
                    if nested.exists() and nested.suffix == ".h5ad":
                        adata_path = nested
                if adata_path is None:
                    found = list(extract_dir.rglob("*.h5ad"))
                    if found:
                        # smallest annotated-looking
                        found.sort(key=lambda p: p.stat().st_size)
                        adata_path = found[0]
    if adata_path is None and h5ads:
        h5ads.sort(key=lambda p: p.stat().st_size)
        adata_path = h5ads[0]
    if adata_path is None:
        log("no park h5ad; writing GSE206710 server-side table as note")
        g = pd.read_csv(ROOT / "module2/tables/M2_r2_GSE206710_donor_fraction.tsv", sep="\t")
        g["atlas"] = "GSE206710_server"
        write_tsv(g, "F4_park_fallback_GSE206710.tsv")
        return g
    log(f"reading park {adata_path} size={adata_path.stat().st_size}")
    import anndata as ad
    adata = ad.read_h5ad(adata_path, backed="r")
    log(f"adata n_obs={adata.n_obs} n_vars={adata.n_vars} obs={list(adata.obs.columns)[:30]}")
    # gene names
    var_names = adata.var_names.astype(str)
    if "gene_symbols" in adata.var.columns:
        symbols = adata.var["gene_symbols"].astype(str)
    elif "symbol" in adata.var.columns:
        symbols = adata.var["symbol"].astype(str)
    else:
        symbols = var_names
    want_idx = []
    want_gene = []
    sym_list = list(symbols)
    for g in genes:
        if g in var_names:
            want_idx.append(list(var_names).index(g) if False else int(np.where(var_names == g)[0][0]))
            want_gene.append(g)
        elif g in set(sym_list):
            want_idx.append(sym_list.index(g))
            want_gene.append(g)
    log("found genes " + str(want_gene))
    sub = adata[:, want_idx].to_memory()
    X = sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)
    expr = pd.DataFrame(X, columns=want_gene, index=sub.obs_names.astype(str))
    obs = sub.obs.copy()
    ctype = pick_col(obs.columns, [
        "Anno_level_fig1", "cell_type", "celltype", "annotation", "anno",
        "CellType", "annot", "Fig1",
    ])
    donor = pick_col(obs.columns, ["donor", "Donor", "sample", "Sample", "donor_id", "individual"])
    age = pick_col(obs.columns, ["Age", "age", "development_stage"])
    log(f"ctype={ctype} donor={donor} age={age}")
    expr["celltype"] = obs[ctype].astype(str).values if ctype else "unknown"
    expr["donor"] = obs[donor].astype(str).values if donor else "unknown"
    if age:
        expr["age"] = obs[age].astype(str).values
    agg = expr.groupby(["donor", "celltype"], as_index=False)[want_gene].mean()
    ncell = expr.groupby(["donor", "celltype"]).size().rename("n_cells").reset_index()
    agg = agg.merge(ncell, on=["donor", "celltype"])
    write_tsv(agg, "F4_park_donor_celltype.tsv")
    # also celltype means (still report n_donors)
    ct = expr.groupby("celltype", as_index=False).agg(
        n_cells=("ZEB1", "size") if "ZEB1" in expr.columns else (want_gene[0], "size"),
        n_donors=("donor", "nunique"),
        **{g: (g, "mean") for g in want_gene if g in expr.columns},
    )
    write_tsv(ct, "F4_park_celltype_means.tsv")
    return agg


def main():
    summary = {"started": time.strftime("%Y-%m-%d %H:%M:%S")}
    park_procs = []
    try:
        park_procs = start_park_download()
    except Exception:
        traceback.print_exc()
        summary["park_start"] = "fail"
    try:
        summary["prism_files"] = [str(p) for p in download_prism()]
    except Exception as e:
        summary["prism"] = str(e)
        traceback.print_exc()
    try:
        summary["gdsc_files"] = [str(p) for p in download_gdsc()]
    except Exception as e:
        summary["gdsc"] = str(e)
        traceback.print_exc()
    try:
        cox = analyze_target_cox()
        summary["cox_rows"] = int(len(cox))
    except Exception as e:
        summary["cox"] = str(e)
        traceback.print_exc()
    try:
        dc = analyze_decoupler()
        summary["decoupler_tfs"] = int(len(dc))
    except Exception as e:
        summary["decoupler"] = str(e)
        traceback.print_exc()
    try:
        dep = load_depmap_tall()
        pts, st = analyze_mtx(dep)
        summary["mtx_points"] = int(len(pts))
        summary["mtx_stats"] = int(len(st))
    except Exception as e:
        summary["mtx"] = str(e)
        traceback.print_exc()
    try:
        wait_procs(park_procs, "park", max_s=10800)
        analyze_park()
        summary["park"] = "ok"
    except Exception as e:
        summary["park"] = str(e)
        traceback.print_exc()
    summary["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    (OUT / "tables" / "RUN_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log("SUMMARY " + json.dumps(summary))
    log("DONE")


if __name__ == "__main__":
    main()
