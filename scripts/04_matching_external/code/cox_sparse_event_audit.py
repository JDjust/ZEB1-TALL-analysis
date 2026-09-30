#!/usr/bin/env python3
"""Count TARGET survival events by subtype and flag Cox separation symptoms."""
from pathlib import Path

import numpy as np
import pandas as pd
import _style as S

surv = pd.read_csv(Path(S.DATA_DIR) / "F7_target_surv_merged.tsv", sep="\t")
cox = pd.read_csv(Path(S.DATA_DIR) / "F7_target_cox.tsv", sep="\t")
events = (surv.groupby("subtype", dropna=False)
          .agg(n_samples=("ZEB1", "size"), os_events=("os_event", "sum"),
               efs_events=("efs_event", "sum"))
          .reset_index())
events["zero_os_events"] = events.os_events.eq(0)
events["zero_efs_events"] = events.efs_events.eq(0)
events.to_csv(Path(S.DATA_DIR) / "SourceData_FigS12_subtype_events.tsv", sep="\t", index=False)
coef = cox[cox.model.eq("ZEB1_subtype")].copy()
coef["nonfinite_or_zero_ci"] = (~np.isfinite(coef.ci_lo) | ~np.isfinite(coef.ci_hi) |
                                 coef.ci_lo.eq(0) | coef.ci_hi.eq(0))
coef.to_csv(Path(S.DATA_DIR) / "SourceData_FigS12_cox_separation_flags.tsv", sep="\t", index=False)
print(events.to_string(index=False))
print("nonfinite/zero subtype-model intervals:", int(coef.nonfinite_or_zero_ci.sum()))
