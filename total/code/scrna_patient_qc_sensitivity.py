"""Fixed-score patient sensitivity across cell-count thresholds and depth covariates."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/validation/scrna_patient_qc'
THRESHOLDS = [30, 100, 250, 500, 1000]
SEED = 20260924


def permute_p(x, y, rng):
    x, y = rankdata(x), rankdata(y)
    x = (x-x.mean()) / np.linalg.norm(x-x.mean())
    y = (y-y.mean()) / np.linalg.norm(y-y.mean())
    null = np.vstack([rng.permutation(y) for _ in range(10000)]) @ x
    return (1 + np.sum(abs(null) >= abs(x @ y)-1e-12)) / 10001


def main():
    qc = pd.read_csv(OUT / 'patient_cell_depth_detection.tsv', sep='\t')
    scores = pd.read_csv(ROOT / 'data/validation/three_gene_patient/patient_three_gene_score.tsv', sep='\t')
    d = scores.merge(qc, on=['cohort','patient'], validate='one_to_one', suffixes=('', '_qc'))
    assert len(d) == 25 and (d.n_cells == d.n_malignant_cells).all()
    assert (d.library_umi == d.summed_umi).all() and (d.ZEB1_count == d.zeb1_umi).all()
    rng = np.random.default_rng(SEED)
    rows, diagnostics = [], []
    for cohort, g in d.groupby('cohort'):
        for minimum in THRESHOLDS:
            h = g[g.n_malignant_cells >= minimum]
            row = dict(cohort=cohort, min_cells=minimum, n_patients=len(h),
                       excluded_patients=','.join(g.loc[g.n_malignant_cells < minimum, 'patient']),
                       score_scaling='frozen full-cohort gene means and SDs', status='evaluated')
            if len(h) < 6:
                row['status'] = 'not estimated: fewer than six patients'
            else:
                subset_seed = int(hashlib.sha256((cohort + '|'.join(sorted(h.patient))).encode()).hexdigest()[:8], 16)
                rng = np.random.default_rng(SEED + subset_seed)
                x, y = h.ZEB1.to_numpy(), h.three_gene_score.to_numpy()
                estimates = []
                for _ in range(2000):
                    ix = rng.integers(0, len(h), len(h))
                    if len(np.unique(x[ix])) >= 3 and len(np.unique(y[ix])) >= 3:
                        estimates.append(spearmanr(x[ix], y[ix]).statistic)
                lo, hi = np.quantile(estimates, [.025, .975])
                row.update(rho=spearmanr(x, y).statistic, permutation_p=permute_p(x, y, rng),
                           bootstrap_ci95_lo=lo, bootstrap_ci95_hi=hi)
            rows.append(row)
        for cov in ['median_umi','median_detected_genes','zeb1_detection_fraction','n_malignant_cells']:
            x, y, c = rankdata(g.ZEB1), rankdata(g.three_gene_score), rankdata(g[cov])
            design = np.column_stack([np.ones(len(c)), c])
            xr = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
            yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
            diagnostics.append(dict(cohort=cohort, covariate=cov, n_patients=len(g),
                zeb1_covariate_rho=spearmanr(g.ZEB1, g[cov]).statistic,
                score_covariate_rho=spearmanr(g.three_gene_score, g[cov]).statistic,
                score_zeb1_partial_rank_rho=np.corrcoef(xr, yr)[0,1],
                interpretation='descriptive sensitivity; no causal adjustment claim or selected model'))
    d.to_csv(OUT / 'patient_qc_scores.tsv', sep='\t', index=False)
    pd.DataFrame(rows).to_csv(OUT / 'cell_threshold_sensitivity.tsv', sep='\t', index=False)
    pd.DataFrame(diagnostics).to_csv(OUT / 'depth_detection_sensitivity.tsv', sep='\t', index=False)
    log_path = OUT / '01_prepare_module5.log'
    log = log_path.read_text(encoding='utf-8', errors='replace')
    required = ['scrublet doublets 121','omicverse harmony OK','leiden igraph OK',
                'Harmony stopped before convergence after 10 iterations']
    assert all(s in log for s in required)
    (OUT / 'analysis_provenance.json').write_text(json.dumps(dict(seed=SEED,
        thresholds=THRESHOLDS, threshold_definition='sensitivity grid defined before correlation recomputation',
        n_permutations=10000, n_bootstrap=2000, tests='exploratory nominal P; no threshold selected',
        upstream_log_sha256=hashlib.sha256(log_path.read_bytes()).hexdigest(),
        actual_pipeline_evidence=required, hulu_job=1730,
        interpretation='Harmony nonconvergence affects label confidence; pseudobulk uses unintegrated raw counts'), indent=2))
    print(pd.DataFrame(rows).to_string(index=False))
    print(pd.DataFrame(diagnostics).to_string(index=False))


if __name__ == '__main__':
    main()
