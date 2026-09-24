"""Link source-paper MRD definitions to the 116-patient analysis table."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/validation'
OUT = DATA / 'clinical_endpoint_audit'


def main():
    OUT.mkdir(exist_ok=True)
    path = ROOT / 'data/reference_records/PMC9873558_fulltext.xml'
    root = ET.parse(path).getroot()
    ids = {e.get('pub-id-type'): e.text for e in root.findall('.//article-meta/article-id')}
    assert ids['doi'] == '10.1038/s41591-022-02112-7' and ids['pmid'] == '36604538'
    names = ['Genomic profiling, MRD and determination of genetic ancestry',
             'Pharmacotyping of primary ALL cells', 'Statistical analyses',
             'Patient and clinical treatment protocols']
    records = []
    for sec in root.findall('.//sec'):
        title_node = sec.find('title')
        title = ''.join(title_node.itertext()) if title_node is not None else ''
        if title in names:
            for i, p in enumerate(sec.findall('p'), 1):
                records.append(dict(section=title, paragraph=i, text=' '.join(p.itertext()),
                                    doi=ids['doi'], pmid=ids['pmid'], pmcid=ids['pmcid']))
    assert set(r['section'] for r in records) == set(names)
    pd.DataFrame(records).to_csv(OUT/'source_method_passages.tsv', sep='\t', index=False)
    d = pd.read_csv(DATA/'pharmacotype_tall_subtype.tsv', sep='\t')
    core = pd.read_csv(ROOT.parent/'modules/module7/tables/M7_r2_core_scores.tsv', sep='\t')
    assert len(d)==116 and not d['Patient ID'].duplicated().any() and not d['sample'].duplicated().any()
    joined = core.merge(d[['sample','Patient ID','Protocol']], on='sample', validate='one_to_one')
    assert len(joined)==116
    early, late = 'Day 15 MRD (%)', 'Day 42 or 46 MRD (%)'
    for col in [early, late]:
        source = d.set_index('sample')[col].reindex(joined['sample']).reset_index(drop=True)
        assert np.allclose(source, joined[col], atol=1e-10, rtol=1e-10, equal_nan=True)
    joined['paired_mrd'] = joined[early].notna() & joined[late].notna()
    joined['mid_positive'] = joined[early].ge(.01).where(joined[early].notna())
    joined['end_positive'] = joined[late].ge(.01).where(joined[late].notna())
    joined['persistence_model_eligible'] = joined.paired_mrd & joined[early].ge(.01)
    summaries = []
    for protocol, g in joined.groupby('Protocol', dropna=False):
        summaries.append(dict(protocol=protocol, n=len(g), early_available=g[early].notna().sum(),
            late_available=g[late].notna().sum(), paired=g.paired_mrd.sum(),
            model_eligible=g.persistence_model_eligible.sum(),
            late_positive_in_model=(g.persistence_model_eligible & g[late].ge(.01)).sum(),
            early_zero=g[early].eq(0).sum(), late_zero=g[late].eq(0).sum(),
            early_positive_below_limit=((g[early]>0)&(g[early]<.01)).sum(),
            late_positive_below_limit=((g[late]>0)&(g[late]<.01)).sum()))
    joined.to_csv(OUT/'patient_protocol_mrd_manifest.tsv', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(OUT/'protocol_endpoint_counts.tsv', sep='\t', index=False)
    (OUT/'source_verification.json').write_text(json.dumps(dict(
        source_doi=ids['doi'], source_pmid=ids['pmid'], source_pmcid=ids['pmcid'],
        fulltext_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        threshold_percent=.01, threshold_status='confirmed source study detection/negativity boundary',
        mid_induction_timing='day 19 Total XV; day 15 Total XVI/XVII',
        end_induction_timing='days 42-46',
        source_below_limit_log_rule='assign half detection limit (0.005 percent) before log transform',
        limitations=['source table zero values are below threshold, not proven complete biological clearance',
                     'individual reasons for missing MRD/drug assays not documented in current table',
                     'current 66-patient model does not adjust protocol, age or WBC; no external validation',
                     'MRD was not a prespecified endpoint of Total XVII according to source']
    ),indent=2))
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == '__main__':
    main()
