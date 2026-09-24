"""Bounded document/output check; does not recompute any scientific analysis."""
from pathlib import Path
from docx import Document
import hashlib, json
root = Path(__file__).resolve().parents[1]
run = sorted((root / 'reproducibility').glob('*_editorial'))[-1]
records = {}
for name, count in [('ZEB1_TALL_manuscript.docx', 7), ('ZEB1_TALL_Supplementary.docx', 15)]:
    path = root / 'manuscript' / name
    doc = Document(path)
    assert len(doc.inline_shapes) == count
    text = '\n'.join(p.text for p in doc.paragraphs)
    assert 'Context-dependent ZEB1 expression and treatment-response associations' in text
    assert '[Missing ' not in text
    assert 'There is no independent clinical validation was available' not in text
    records[name] = {'embedded_figures': count, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
for path in [root/'figures/Figure2_development.pdf', root/'figures/Figure6_mechanism.pdf', root/'figures/main_overview.png', root/'figures/supplementary_overview.png']:
    assert path.stat().st_size > 1000
    records[str(path.relative_to(root))] = {'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
(run/'output_verification.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
(run/'README.md').write_text('''# Editorial polish

Existing analyses retained. No data download, server job or full package rebuild.
Rendered Figure 2 and Figure 6 only; regenerated main/supplementary documents and 7/15-figure overview sheets.
Manuscript: concise shared title, 246-word abstract, biological-context and MTX/MRD discussion, neutral evidence headings, corrected Figure 7 legend grammar.
Layout: shorter panel titles, headings kept with following text, figure paragraphs kept with captions, supplementary image width reduced to fit the text area.
Verification: DOCX opened programmatically; 7 and 15 embedded images, shared title and missing-image checks passed. Main/supplementary overview sheets and Figure 6 inspected visually. Word page-by-page rendering was not performed.
Rollback: the original four edited source files and both previous DOCX files are preserved in this directory. Restore source files and rerender affected outputs to undo.
''', encoding='utf-8')
print(json.dumps(records, indent=2))
