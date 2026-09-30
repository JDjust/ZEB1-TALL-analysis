"""A6-2 quadrant: normal Δ vs leukemia β_R. Descriptive, not a P-value figure."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a6_yayon")
g = pd.read_csv(P / "a6_2_gene_effects.tsv", sep="\t")
col = {"concordant": "#2C7BB6", "discordant": "#D7191C", "neutral": "#BDBDBD"}
fig, ax = plt.subplots(figsize=(5.2, 5.0), dpi=200)
ax.axhline(0, color="#666666", lw=0.6)
ax.axvline(0, color="#666666", lw=0.6)
for cls, sub in g.groupby("class"):
    ax.scatter(
        sub["delta_normal"], sub["beta_R"],
        s=22, c=col[cls], alpha=0.85, linewidths=0, label=f"{cls} n={len(sub)}", zorder=3,
    )
ax.set_xlabel("Normal Δ (DP − early)")
ax.set_ylabel("Leukemia residual β")
ax.legend(frameon=False, fontsize=8, loc="upper right")
ax.set_title("A6-2  ρ = −0.50", fontsize=10)
fig.tight_layout()
fig.savefig(P / "a6_2_quadrant.png", bbox_inches="tight")
print("wrote", P / "a6_2_quadrant.png", flush=True)
