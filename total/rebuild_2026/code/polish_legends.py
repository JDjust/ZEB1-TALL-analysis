"""Keep journal legends focused on panel data; provenance remains in manifests."""
from pathlib import Path
import re

base=Path(__file__).resolve().parents[1]
man=base/'manuscript'
legends=[]
for i in range(1,8):
    path=man/f'Figure_{i}_legend.md'
    text=path.read_text(encoding='utf-8').strip()
    text=re.sub(r'\s*Source data: `data/source_data_rebuilt/[^`]+/`[.;]\s*(?:Plot code|assembly code): `code/[^`]+`\.', '', text)
    text=text.replace('previously frozen ','').replace('Frozen ','').replace('frozen ','')
    path.write_text(text+'\n',encoding='utf-8')
    legends.append(text)
main_legend='## Figure legends\n\n'+'\n\n'.join(legends)+'\n\n'
(man/'Figure_Legends.md').write_text('# Figure legends\n\n'+'\n\n'.join(legends)+'\n',encoding='utf-8')
main=man/'ZEB1_ZEB2_TALL_manuscript.md'
text=main.read_text(encoding='utf-8')
start=text.index('## Figure legends\n')
end=text.index('## References\n',start)
main.write_text(text[:start]+main_legend+text[end:],encoding='utf-8')

supp=man/'Supplementary_Information.md'
text=supp.read_text(encoding='utf-8')
text=re.sub(r' Panel source files and original R scripts are enumerated in `data/source_data_rebuilt/Supplement_S1/S1_panel_manifest.tsv`\.', '', text)
text=re.sub(r' Exact source paths are in `data/source_data_rebuilt/Supplement_S2/S2_panel_manifest.tsv`\.', '', text)
text=re.sub(r' Source paths are in `data/source_data_rebuilt/Supplement_S[3-7]/S[3-7]_panel_manifest.tsv`\.', '', text)
start=text.index('## Supplementary Figure S1.')
end=text.index('## Supplementary Table S1.',start)
section=text[start:end].replace('previously frozen ','').replace('Frozen ','').replace('frozen ','')
text=text[:start]+section+text[end:]
supp.write_text(text,encoding='utf-8')
end=text.index('## Supplementary Table S1.',start)
(man/'Supplementary_Figure_Legends_consolidated.md').write_text('# Supplementary figure legends\n\n'+text[start:end].strip()+'\n',encoding='utf-8')
print('Main and supplementary legends revised; panel provenance remains in manifests')
