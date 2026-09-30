"""Run one build script while recording project file opens via Python audit events."""
from pathlib import Path
import json
import os
import runpy
import sys

script = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve()
root = script.parents[1]
project = root.parent
seen = set()


def audit(event, args):
    if event != "open" or not isinstance(args[0], (str, bytes, os.PathLike)):
        return
    try:
        path = Path(os.fsdecode(args[0])).resolve()
        relative = path.relative_to(project)
    except (ValueError, OSError):
        return
    if "reproducibility" in relative.parts or "__pycache__" in relative.parts:
        return
    mode = args[1]
    flags = args[2]
    write = any(c in mode for c in "wax+") if isinstance(mode, str) else bool(flags & (os.O_WRONLY | os.O_RDWR))
    seen.add((os.path.relpath(path, root).replace("\\", "/"), "write" if write else "read"))


sys.addaudithook(audit)
sys.path.insert(0, str(script.parent))
sys.argv = [str(script)]
try:
    runpy.run_path(str(script), run_name="__main__")
finally:
    output.write_text(json.dumps({"script": script.name,
                                 "scope": "observed Python project file-open events; not a raw-data lineage proof",
                                 "files": [{"path": p, "access": m} for p, m in sorted(seen)]},
                                indent=2), encoding="utf-8")
