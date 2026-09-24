"""Rebuild all publication artifacts from the current exported analysis tables.

This does not claim to rerun raw sequencing analyses. Each child has a timeout,
log, heartbeat and persistent status; failed children stop the build.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    "build_dataset_registry.py",
    "build_reference_catalog.py",
    "fig4_crosscohort_source.py",
    "tlx_meta_sensitivity.py",
    "scrna_patient_qc_sensitivity.py",
    "matched_gene_set_control.py",
    "audit_pharmacotype_endpoints.py",
    "mrd_protocol_sensitivity.py",
    "summarize_harmony_sensitivity.py",
    "thymus_reference_validation.py",
    "run_target_survival_diagnostics.py",
    "gse144035_paired_audit.py",
    "fig1_disease_specificity.py", "fig2_development.py", "fig3_subtype.py",
    "fig4_transcriptome.py", "fig5_singlecell.py", "fig6_mechanism.py",
    "fig7_therapy.py", "figS_all.py", "figS13_probe_annotation.py",
    "figS14_target_qc.py", "figS15_three_gene_patient.py",
    "build_claims_registry.py",
    "build_manuscript.py", "build_supplementary.py", "figure_overview.py",
]
RENDER_STEPS = STEPS[STEPS.index('fig1_disease_specificity.py'):]


def write_json(path, obj):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    temp.replace(path)


def hashes(paths):
    records = []
    for p in sorted(set(paths)):
        digest = hashlib.sha256()
        with p.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        records.append({"path": os.path.relpath(p, ROOT).replace("\\", "/"),
                        "bytes": p.stat().st_size, "sha256": digest.hexdigest()})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-timeout", type=int, default=180)
    parser.add_argument("--refresh-derived", action="store_true",
                        help="Explicitly rerun legacy derived analyses before rendering; default reuses frozen results")
    parser.add_argument("--list-steps", action="store_true",
                        help="Print the selected steps without executing or writing files")
    args = parser.parse_args()
    steps = STEPS if args.refresh_derived else RENDER_STEPS
    if args.list_steps:
        print('\n'.join(steps))
        return
    if args.step_timeout < 30:
        parser.error("--step-timeout must be at least 30 seconds")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run = ROOT / "reproducibility" / stamp
    run.mkdir(parents=True, exist_ok=False)
    packages = ["numpy", "pandas", "matplotlib", "scipy", "python-docx", "Pillow"]
    environment = {"python": sys.version, "executable": sys.executable,
                   "platform": platform.platform(),
                   "packages": {p: importlib.metadata.version(p) for p in packages}}
    write_json(run / "environment.json", environment)
    inputs = list((ROOT / "data").rglob("*.tsv")) + list((ROOT / "data").rglob("*.json"))
    inputs += list((ROOT / "data").rglob("*.csv"))
    inputs += list((ROOT.parent / "modules").glob("*/tables/*"))
    inputs += list((ROOT / "code").glob("*.py"))
    inputs += list((ROOT / "code").glob("*.R"))
    inputs += list((ROOT / "data/validation/target_survival").glob("*.xlsx"))
    write_json(run / "input_export_hashes.json", hashes(p for p in inputs if p.is_file()))
    state = {"scope": "exported-tables-to-figures-and-documents", "started_utc": stamp,
             "refresh_derived": args.refresh_derived,
             "status": "running", "steps": [], "step_timeout_seconds": args.step_timeout}
    status_path = run / "status.json"
    write_json(status_path, state)
    try:
        for script in steps:
            record = {"script": script, "status": "running"}
            state["steps"].append(record)
            write_json(status_path, state)
            print(f"START {script}", flush=True)
            started = time.monotonic()
            env = dict(os.environ, PYTHONIOENCODING="utf-8", MPLBACKEND="Agg")
            with (run / f"{Path(script).stem}.log").open("w", encoding="utf-8") as log:
                child = subprocess.Popen([sys.executable, str(ROOT / "code/run_traced_step.py"),
                                          str(ROOT / "code" / script),
                                          str(run / f"{Path(script).stem}.file_access.json")],
                                         cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, env=env)
                while True:
                    try:
                        code = child.wait(timeout=min(20, args.step_timeout))
                        break
                    except subprocess.TimeoutExpired:
                        elapsed = time.monotonic() - started
                        print(f"RUNNING {script}: {elapsed:.0f}s; log={log.name}", flush=True)
                        if elapsed >= args.step_timeout:
                            child.kill()
                            child.wait()
                            raise TimeoutError(f"{script} exceeded {args.step_timeout}s")
            record.update(status="passed" if code == 0 else "failed",
                          exit_code=code, elapsed_seconds=round(time.monotonic() - started, 2))
            write_json(status_path, state)
            if code != 0:
                raise RuntimeError(f"{script} failed; inspect its log in {run}")
            print(f"PASS {script}", flush=True)
        from docx import Document
        main_doc = ROOT / "manuscript/ZEB1_TALL_manuscript.docx"
        supp_doc = ROOT / "manuscript/ZEB1_TALL_Supplementary.docx"
        counts = {"main": len(Document(main_doc).inline_shapes),
                  "supplementary": len(Document(supp_doc).inline_shapes)}
        if counts != {"main": 7, "supplementary": 15}:
            raise RuntimeError(f"Unexpected embedded figure counts: {counts}")
        outputs = list((ROOT / "figures").rglob("*.pdf"))
        outputs += list((ROOT / "figures/png").glob("*.png"))
        outputs += list((ROOT / "figures/Supplementary/png").glob("*.png"))
        outputs += [main_doc, supp_doc]
        write_json(run / "output_hashes.json", hashes(outputs))
        state.update(status="passed", embedded_figures=counts)
    except BaseException as exc:
        state.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        if state["steps"] and state["steps"][-1]["status"] == "running":
            state["steps"][-1]["status"] = "failed"
        raise
    finally:
        write_json(status_path, state)
        print(f"BUILD {state['status']}: {status_path}", flush=True)


if __name__ == "__main__":
    main()
