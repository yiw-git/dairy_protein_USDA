"""Replot the DLVO stability map in the ACS2026 figure style.

Reads the existing sweep output (no re-running of the DLVO model) and
re-renders stability_map.png with the deck's typography: bold dark headline
title, transparent background, gold/green kT contours, larger fonts.
"""
import json
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- ACS2026 deck palette -------------------------------------------------
RED="#E21833"; DARK="#282828"; GREY="#636363"; GOLD="#E9A200"; GREEN="#2E7D32"
BLUE="#1F6FB2"

OUT = "../../../scripts_and_methods/DLVO/outputs/real_data_08102026"  # relative to ACS2026/figures
DEST = "dlvo_stability_map.png"
MARGINAL_KT, STABLE_KT = 10.0, 15.0

sweep = pd.read_csv(f"{OUT}/stability_sweep_grid.csv")
tg    = pd.read_csv(f"{OUT}/target_condition_predictions.csv")
cs    = pd.read_csv(f"{OUT}/condition_summary.csv")
info  = json.load(open(f"{OUT}/zeta_model_info.json"))

pivot = sweep.pivot(index="ionic_strength_mM", columns="pH", values="V_max_kT")
X, Y, Z = pivot.columns.values, pivot.index.values, pivot.values

# half-slide panel: same 400 dpi / transparent convention as the deck figures
fig, ax = plt.subplots(figsize=(4.6, 3.15), dpi=400)

vmax = max(30.0, np.nanpercentile(Z, 99))
im = ax.pcolormesh(X, Y, Z, shading="auto", cmap="RdYlGn", vmin=-10, vmax=vmax)

ct = ax.contour(X, Y, Z, levels=[MARGINAL_KT, STABLE_KT],
                colors=[GOLD, GREEN], linewidths=1.3)
ax.clabel(ct, inline=True, fontsize=8.5, fmt="%.0f kT",
          manual=[(6.45, 41.5), (6.45, 29.5)])

ax.set_yscale("log")
ax.set_xlabel("pH", fontsize=10, labelpad=2.0)
ax.set_ylabel("ionic strength (mM)", fontsize=10, labelpad=2.0)
ax.tick_params(labelsize=9, length=2.5, pad=1.5)
ax.tick_params(which='minor', length=1.5)
for sp in ['top', 'right']:
    ax.spines[sp].set_visible(False)

# food target conditions (label offsets keep text inside the axes)
LBL = {"whole milk": (-9, 5, "right"), "yogurt": (8, 4, "left"),
       "bread-release": (8, 4, "left")}
for _, r in tg.iterrows():
    dx, dy, ha = LBL.get(r.condition, (8, 4, "left"))
    ax.scatter(r.pH, r.ionic_strength_mM, marker="*", s=170,
               color=BLUE, edgecolor="white", linewidth=0.7, zorder=5)
    ax.annotate(r.condition, (r.pH, r.ionic_strength_mM), textcoords="offset points",
                xytext=(dx, dy), fontsize=9, color=BLUE, ha=ha, zorder=6)

# measured condition that aggregated (radius outlier)
for _, r in cs[cs.radius_outlier].iterrows():
    ax.scatter(r.pH, r.ionic_strength_mM, marker="X", s=130,
               color=RED, edgecolor="white", linewidth=0.7, zorder=6)
    ax.annotate("large aggregates\n(near CCC)", (r.pH, r.ionic_strength_mM),
                textcoords="offset points", xytext=(9, -13), fontsize=8.5, color=RED)

cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.055)
cbar.set_label(r"energy barrier $V_{max}$ ($k_BT$)", fontsize=10, labelpad=3)
cbar.ax.tick_params(labelsize=9, length=2.5, pad=1.5)
cbar.outline.set_linewidth(0.6)

ax.set_title("Food targets sit outside the stable window",
             fontsize=11, fontweight='bold', color=DARK, pad=5)
ax.text(0.0, -0.175, f"predicted from measured ζ(pH, I) · a = {info['radius_nm_used_for_sweep']:.0f} nm · "
        r"$A$ = 10$^{-20}$ J · DLS 08/10/2026",
        transform=ax.transAxes, fontsize=7.5, color=GREY, va='top', ha='left')

fig.savefig(DEST, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved", DEST)
