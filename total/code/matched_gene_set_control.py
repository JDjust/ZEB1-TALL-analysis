"""Competitive three-gene controls matched without using ZEB1 correlations.

Uses processed within-cohort log expression. This exploratory competitive
comparison tests expression-program specificity, not biological causality.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/validation'
OUT = DATA / 'matched_gene_sets'
GENES = ['GATA3', 'TCF7', 'BCL11B']
SEED = 20260924
N = 10000
POOL = 100


def main():
    OUT.mkdir(exist_ok=True)
    summaries, nulls, pools, metrics, patient_scores = [], [], [], [], []
    sources = []
    for cohort, name, expected_n in [
        ('TARGET', 'target_log2tpm_filtered.tsv.gz', 265),
        ('Pharmacotype', 'pharmacotype_log2fpkm_filtered.tsv.gz', 116)]:
        path = DATA / name
        frame = pd.read_csv(path, sep='\t', index_col=0)
        assert not frame.index.duplicated().any() and frame.shape[1] == expected_n
        x = frame.to_numpy(dtype=float)
        assert np.isfinite(x).all() and (x >= 0).all()
        avg, detection, sd = x.mean(axis=1), (x > 0).mean(axis=1), x.std(axis=1, ddof=1)
        valid = sd > 0
        frame = frame.loc[valid]
        x, avg, detection, sd = x[valid], avg[valid], detection[valid], sd[valid]
        z = (x - avg[:, None]) / sd[:, None]
        stats = pd.DataFrame(dict(gene=frame.index, mean_log_expression=avg,
                                  detection_fraction=detection, sd_log_expression=sd))
        stats['cohort'] = cohort
        metrics.append(stats)
        zi = frame.index.get_loc('ZEB1')
        focus = [frame.index.get_loc(g) for g in GENES]
        score = z[focus].mean(axis=0)
        obs = float(spearmanr(x[zi], score).statistic)
        zr = rankdata(x[zi]); zr = (zr-zr.mean()) / np.linalg.norm(zr-zr.mean())
        patient_scores.append(pd.DataFrame(dict(cohort=cohort, patient=frame.columns,
                                                ZEB1=x[zi], three_gene_score=score)))
        for scheme, columns in [('mean_detection', [avg,detection]),
                                ('mean_detection_sd', [avg,detection,sd])]:
            # Percentile ranks place covariates on comparable scales; tied
            # detection fractions receive their common midrank.
            features = np.column_stack([rankdata(v)/len(v) for v in columns])
            eligible = ~frame.index.isin(GENES + ['ZEB1'])
            candidates = []
            for gene, fi in zip(GENES, focus):
                distance = np.linalg.norm(features-features[fi], axis=1)
                distance[~eligible] = np.inf
                ix = np.argsort(distance, kind='stable')[:POOL]
                assert len(ix) == POOL and np.isfinite(distance[ix]).all()
                candidates.append(ix)
                for k in ix:
                    pools.append(dict(cohort=cohort, scheme=scheme, target_gene=gene,
                                      matched_gene=frame.index[k], percentile_distance=distance[k],
                                      mean_log_expression=avg[k], detection_fraction=detection[k],
                                      sd_log_expression=sd[k], target_mean=avg[fi],
                                      target_detection=detection[fi], target_sd=sd[fi]))
            offset = int(hashlib.sha256((cohort+scheme).encode()).hexdigest()[:8], 16)
            rng = np.random.default_rng(SEED+offset)
            draws = []
            while len(draws) < N:
                draw = [int(rng.choice(c)) for c in candidates]
                if len(set(draw)) == 3:
                    draws.append(draw)
            draws = np.array(draws)
            rs = []
            for start in range(0, N, 500):
                scores = z[draws[start:start+500]].mean(axis=1)
                ranked = rankdata(scores, axis=1)
                ranked -= ranked.mean(axis=1, keepdims=True)
                ranked /= np.linalg.norm(ranked, axis=1, keepdims=True)
                rs.extend(ranked @ zr)
            rs = np.asarray(rs)
            empirical = float((1+np.sum(rs >= obs-1e-12))/(N+1))
            summaries.append(dict(cohort=cohort, matching=scheme, n_patients=expected_n,
                n_background_genes=len(frame), n_sets=N, candidate_pool_per_gene=POOL,
                observed_rho=obs, null_median=np.median(rs), null_q025=np.quantile(rs,.025),
                null_q975=np.quantile(rs,.975), competitive_upper_tail_p=empirical,
                interpretation='exploratory competitive positive-tail comparison; not patient-null P'))
            nulls.append(pd.DataFrame(dict(cohort=cohort, matching=scheme,
                set_index=np.arange(1,N+1), rho=rs,
                gene_1=frame.index.to_numpy()[draws[:,0]],
                gene_2=frame.index.to_numpy()[draws[:,1]],
                gene_3=frame.index.to_numpy()[draws[:,2]])))
        sources.append(dict(file=name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    pd.DataFrame(summaries).to_csv(OUT/'competitive_summary.tsv', sep='\t', index=False)
    pd.concat(nulls).to_csv(OUT/'random_sets.tsv.gz', sep='\t', index=False)
    pd.DataFrame(pools).to_csv(OUT/'matching_pools.tsv', sep='\t', index=False)
    pd.concat(metrics).to_csv(OUT/'gene_matching_features.tsv', sep='\t', index=False)
    pd.concat(patient_scores).to_csv(OUT/'observed_patient_scores.tsv', sep='\t', index=False)
    (OUT/'run.json').write_text(json.dumps(dict(seed=SEED, n_draws=N, pool_size=POOL,
        genes=GENES, sources=sources, matching='nearest percentile-rank feature distance; no ZEB1 correlation used',
        duplicates='no repeated gene within a set; sets drawn with replacement',
        scope='post hoc competitive specificity control on previously selected genes; not independent validation'), indent=2))
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == '__main__':
    main()
