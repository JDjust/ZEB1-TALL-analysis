"""Audit saved diagnosis-cell counts; no new malignant labels or cell exclusions."""
import argparse
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from gse248287_malignant_pseudobulk import load_sample


def summarize(cohort, patient, counts, genes, selected, n_all, label):
    # Canonical orientation is cells by genes.
    x = counts[selected, :].tocsr()
    depth = np.asarray(x.sum(axis=1)).ravel()
    detected = np.asarray((x > 0).sum(axis=1)).ravel()
    zeb = np.asarray(x[:, pd.Index(genes).get_loc('ZEB1')].toarray()).ravel()
    return dict(cohort=cohort, patient=patient, n_diagnosis_cells=n_all,
                n_malignant_cells=x.shape[0], malignant_label_source=label,
                median_umi=float(np.median(depth)), median_detected_genes=float(np.median(detected)),
                zeb1_detected_cells=int((zeb > 0).sum()), zeb1_detection_fraction=float((zeb > 0).mean()),
                summed_umi=int(depth.sum()), zeb1_umi=int(zeb.sum()),
                qc_scope='retained processed diagnosis cells; not a raw-to-QC exclusion count')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    base = args.project / 'module5/processed'
    a = ad.read_h5ad(base / 'GSE227122_annotated.h5ad')
    rows = []
    dx = a.obs.timepoint.eq('Dx')
    for patient in sorted(a.obs.loc[dx, 'patient'].unique()):
        all_patient = dx & a.obs.patient.eq(patient)
        selected = (all_patient & a.obs.malignant.eq('malignant')).to_numpy()
        if selected.any():
            rows.append(summarize('GSE227122', str(patient), a.layers['counts'], a.var_names,
                                  selected, int(all_patient.sum()), 'project heuristic'))
    schema = {'GSE227122_obs_columns': list(a.obs.columns),
              'GSE227122_uns_keys': list(a.uns.keys()),
              'GSE227122_obsm_keys': list(a.obsm.keys()),
              'GSE227122_layers': list(a.layers.keys())}
    del a
    for folder in sorted((base / 'gse248287_dx').glob('P*_Dx')):
        meta, features, counts = load_sample(folder)
        rows.append(summarize('GSE248287', folder.name.split('_')[0], counts.T, features,
                              meta.Celltypes_all.eq('Malignant').to_numpy(), len(meta), 'author malignant'))
        schema['GSE248287_metadata_columns'] = list(meta.columns)
    args.output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output / 'patient_cell_depth_detection.tsv', sep='\t', index=False)
    (args.output / 'saved_object_schema.json').write_text(json.dumps(schema, indent=2))
    print(pd.DataFrame(rows).groupby('cohort').agg(n=('patient','size'),
          cells_min=('n_malignant_cells','min'), cells_max=('n_malignant_cells','max'),
          zeb1_detection_min=('zeb1_detection_fraction','min'),
          zeb1_detection_max=('zeb1_detection_fraction','max')).to_string(), flush=True)


if __name__ == '__main__':
    main()
