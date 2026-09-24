"""Run the R survival diagnostics, then audit source fields and prior coefficients."""
from pathlib import Path
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
rscript = os.environ.get('RSCRIPT') or shutil.which('Rscript')
if not rscript and Path('D:/R/R-4.6.0/bin/Rscript.exe').is_file():
    rscript = 'D:/R/R-4.6.0/bin/Rscript.exe'
if not rscript:
    raise RuntimeError('Set RSCRIPT to an Rscript executable with the survival package installed.')
subprocess.run([rscript,str(ROOT/'code/target_survival_diagnostics.R'),str(ROOT)],check=True,timeout=90)
subprocess.run([sys.executable,str(ROOT/'code/audit_target_survival.py')],check=True,timeout=90)
