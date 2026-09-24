"""Bundle byte-identical module tables used by figure scripts, without recomputation."""
from pathlib import Path
import ast, hashlib, json, shutil
from datetime import datetime, timezone
root = Path(__file__).resolve().parents[1]
references = {}
for script in (root/'code').glob('fig*.py'):
    for node in ast.walk(ast.parse(script.read_text(encoding='utf-8-sig'))):
        if not isinstance(node,ast.Call): continue
        is_tbl = isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name) and node.func.value.id=='S' and node.func.attr=='tbl'
        is_helper = isinstance(node.func,ast.Name) and node.func.id=='_p'
        if (is_tbl or is_helper) and len(node.args)==2 and all(isinstance(a,ast.Constant) and isinstance(a.value,str) for a in node.args):
            key=tuple(a.value for a in node.args)
            references.setdefault(key,[]).append(f'{script.name}:{node.lineno}')
rows=[]
for (module,name), uses in sorted(references.items()):
    source=root.parent/'modules'/module/'tables'/name
    target=root/'data/module_tables'/module/name
    assert source.is_file(),source
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,target)
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert digest(source)==digest(target)
    rows.append({'module':module,'filename':name,'source':str(source),'bundled':str(target.relative_to(root)),
                 'bytes':target.stat().st_size,'sha256':digest(target),'figure_calls':uses})
manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':'Literal S.tbl and _p calls in figure scripts; not raw-data or whole-pipeline packaging','files':rows}
(root/'data/module_tables/manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(f'Bundled {len(rows)} unchanged tables; {sum(x["bytes"] for x in rows):,} bytes.')
