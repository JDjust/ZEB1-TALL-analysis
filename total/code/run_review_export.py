"""Bounded Word export; terminates only the recorded dedicated Word instance."""
from pathlib import Path
import subprocess, time, sys
import psutil
root = Path(__file__).resolve().parents[1]
name = sys.argv[1] if len(sys.argv)>1 else 'ZEB1_TALL_manuscript'
assert name in ['ZEB1_TALL_manuscript', 'ZEB1_TALL_Supplementary']
pidfile = root/'manuscript/review_pdf'/f'{name}.wordpid'
if pidfile.exists(): pidfile.unlink()
started = time.time()
layout_only = '--layout-only' in sys.argv
command = ['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(root/'code/export_review_pdf.ps1'),'-DocumentName',name]
if layout_only: command.append('-LayoutOnly')
timeout = 60 if layout_only else 150
p = subprocess.Popen(command)
try:
    code = p.wait(timeout=timeout)
except subprocess.TimeoutExpired:
    print(f'Word operation exceeded {timeout} seconds; terminating recorded automation instance.', flush=True)
    if pidfile.exists():
        try:
            word = psutil.Process(int(pidfile.read_text(encoding='utf-8-sig').strip()))
            if word.name().lower()=='winword.exe' and word.create_time() >= started-2:
                word.terminate()
        except (psutil.NoSuchProcess, ValueError): pass
    p.terminate()
    p.wait(timeout=10)
    sys.exit(124)
sys.exit(code)
