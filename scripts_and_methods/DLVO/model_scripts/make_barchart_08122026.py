import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = "/sessions/youthful-gifted-dijkstra/mnt/dairy_protein_USDA/scripts_and_methods/DLVO/outputs/real_data_08122026"

cs = pd.read_csv(f"{OUT}/condition_summary.csv")
tg = pd.read_csv(f"{OUT}/target_condition_predictions.csv")

rows = []
for _, r in cs.iterrows():
    label = f"pH {r['pH']:.1f}, {r['ionic_strength_mM']:.0f} mM"
    if r["radius_outlier"]:
        label += " *"
    rows.append({"label": label, "V_max_kT": r["V_max_kT"], "group": "measured grid"})
for _, r in tg.iterrows():
    rows.append({"label": r["condition"], "V_max_kT": r["V_max_kT"], "group": "target food condition"})

df = pd.DataFrame(rows).sort_values("V_max_kT").reset_index(drop=True)

XMIN, XMAX = -4, 16
df["plot_val"] = df["V_max_kT"].clip(lower=XMIN + 0.3)

fig, ax = plt.subplots(figsize=(8, 7.5))
colors = ["#c0392b" if g == "measured grid" else "#2980b9" for g in df["group"]]
bars = ax.barh(df["label"], df["plot_val"], color=colors)

for i, r in df.iterrows():
    if r["V_max_kT"] < XMIN:
        ax.annotate(f"{r['V_max_kT']:.1f} kT", xy=(XMIN + 0.3, i), xytext=(XMIN + 0.4, i),
                    va="center", ha="left", fontsize=7, color="white", fontweight="bold")

ax.axvline(15, color="green", linestyle="--", linewidth=1.2)
ax.axvline(10, color="orange", linestyle="--", linewidth=1.2)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlim(XMIN, XMAX)
ax.set_xlabel(r"Energy barrier $V_{max}$ ($k_BT$)  (bars truncated at -4 kT; true value labeled)")
ax.set_title("Zein-caseinate NP DLVO stability -- 08/12/2026 (single replicate, preliminary)\nred = measured pH x ionic-strength grid, blue = target food conditions", fontsize=10)

legend_elems = [
    Patch(facecolor="#c0392b", label="measured grid condition"),
    Patch(facecolor="#2980b9", label="target food condition"),
    plt.Line2D([0], [0], color="green", linestyle="--", label="stable threshold (15 kT)"),
    plt.Line2D([0], [0], color="orange", linestyle="--", label="marginal threshold (10 kT)"),
]
ax.legend(handles=legend_elems, fontsize=8, loc="upper right")
ax.text(0.99, 0.01, "* pH4/100mM: radius outlier (large aggregates near CCC)", transform=ax.transAxes,
        fontsize=7, ha="right", va="bottom", color="gray")

fig.tight_layout()
fig.savefig(f"{OUT}/condition_barrier_barchart.png", dpi=150)
print("saved")
