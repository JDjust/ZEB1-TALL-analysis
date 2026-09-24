"""Import only author-supplied front matter and declarations; retain provenance."""
from pathlib import Path
import hashlib
import json
from docx import Document

ROOT = Path(__file__).resolve().parents[1]
source = Path('D:/ZEB1/ZEB1-修改终版.docx')
paragraphs = [p.text for p in Document(source).paragraphs]
fields = {'authors': 2, 'affiliations': 3, 'correspondence': 4,
          'contributions': 8, 'acknowledgments': 82, 'funding': 88,
          'competing_interests': 90, 'ethics': 92}
record = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
          'scope': 'Author-supplied metadata only; scientific content not imported',
          'paragraph_indices_zero_based': fields,
          'verbatim': {key: paragraphs[index] for key, index in fields.items()},
          'remaining_confirmation': [
              'Meaning of # marks beside Qian Zhou and Qing yun Ni; no legend in source',
              'Qing yun Ni versus Qingyun Ni spelling differs within source',
              'Final approval applies to the newly revised manuscript',
              'Source ethics statement applicability to expanded datasets and access terms',
              'Public repository/access instructions and AI-use disclosure']}
out = ROOT / 'data/author_metadata.json'
out.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
print(f'Imported {len(fields)} metadata fields; source SHA256={record["source_sha256"]}')
