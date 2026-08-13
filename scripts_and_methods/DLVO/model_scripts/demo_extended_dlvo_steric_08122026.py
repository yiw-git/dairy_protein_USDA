"""
demo_extended_dlvo_steric_08122026.py
=======================================

DEMO / illustrative run of extended DLVO (classical V_EDL + V_vdW, plus a
new V_steric term for the sodium caseinate brush layer) on the 08/12/2026
dataset, using LITERATURE-TYPICAL values for the two new inputs the
steric term needs (layer thickness L, adsorbed surface density Gamma) --
NOT measured on this system. See dlvo_model.py's new steric-term comment
block for the physics, derivation, and validation of V_steric(D) itself.

Literature ranges used (sodium caseinate adsorbed layer at an interface;
see conversation / dlvo_model.py for full citation):
    L (layer thickness)      : 10 - 20 nm   -> low/mid/high = 10/15/20 nm
    Gamma (surface density)  : 1 - 3 mg/m^2 -> low/mid/high = 1/2/3 mg/m^2
"low/mid/high" pairs L and Gamma together (both scale with "more protein
adsorbed" in the same direction), giving three illustrative scenarios
rather than a full 3x3 sweep.

This is explicitly a sensitivity demo, not a final result -- swap in your
own measured L and Gamma (see dlvo_model.py docstring for how) once
available.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = "/sessions/youthful-gifted-dijkstra/mnt/dairy_protein_USDA/scripts_and_methods/DLVO"
sys.path.insert(0, os.path.join(HERE, "model_scripts"))
import dlvo_model as dlvo

SRC_DIR = os.path.join(HERE, "outputs", "real_data_08122026")
OUTDIR = os.path.join(HERE, "outputs", "real_data_08122026", "extended_dlvo_demo")
os.makedirs(OUTDIR, exist_ok=True)

HAMAKER_J = dlvo.HAMAKER_ZEIN_DEFAULT_J
STABLE_KT, MARGINAL_KT = 15.0, 10.0

SCENARIOS = {
    "low":  {"L_nm": 10.0, "gamma_mg_m2": 1.0},
    "mid":  {"L_nm": 15.0, "gamma_mg_m2": 2.0},
    "high": {"L_nm": 20.0, "gamma_mg_m2": 3.0},
}
for name, sc in SCENARIOS.items():
    sc["s_nm"] = dlvo.anchor_spacing_nm(sc["gamma_mg_m2"])
    print(f"Scenario '{name}': L={sc['L_nm']} nm, Gamma={sc['gamma_mg_m2']} mg/m^2 "
          f"-> anchor spacing s={sc['s_nm']:.2f} nm")

# ---------------------------------------------------------------------------
# Load 08/12 measured conditions + target food conditions
# ---------------------------------------------------------------------------
cs = pd.read_csv(f"{SRC_DIR}/condition_summary.csv")
tg = pd.read_csv(f"{SRC_DIR}/target_condition_predictions.csv")

rows = []
for _, r in cs.iterrows():
    rows.append({"label": f"pH {r['pH']:.1f}, {r['ionic_strength_mM']:.0f} mM" + (" *" if r["radius_outlier"] else ""),
                 "pH": r["pH"], "I_mM": r["ionic_strength_mM"], "radius_nm": r["radius_nm"],
                 "zeta_mV": r["zeta_mV"], "group": "measured grid", "radius_outlier": r["radius_outlier"],
                 "V_classic_kT": r["V_max_kT"]})
# representative radius for target conditions (same one used to generate target_condition_predictions.csv)
import json
with open(f"{SRC_DIR}/zeta_model_info.json") as f:
    model_info = json.load(f)
radius_repr = model_info["radius_nm_used_for_sweep"]
for _, r in tg.iterrows():
    rows.append({"label": r["condition"], "pH": r["pH"], "I_mM": r["ionic_strength_mM"],
                 "radius_nm": radius_repr, "zeta_mV": r["zeta_mV_predicted"], "group": "target food condition",
                 "radius_outlier": False, "V_classic_kT": r["V_max_kT"]})

df = pd.DataFrame(rows)
kT_J = dlvo.K_B * 298.15

# IMPORTANT FINDING (discovered while validating this demo, see chat):
# at these literature-typical L/Gamma values, the steric term is strong
# enough that V(D) is monotonically repulsive all the way down to the
# D_min=0.3 nm hard-wall cutoff for EVERY condition tested here -- i.e.
# find_energy_barrier() never finds a genuine interior local maximum
# (barrier_present=False throughout). The raw "V_max" it returns in that
# case is just V(D_min), an artifact of the arbitrary 0.3 nm cutoff
# (calibrated for the classical Born-repulsion case, not appropriate once
# a 10-20 nm brush layer is present) -- NOT a physically meaningful
# barrier height, and not comparable in magnitude to the classical V_max
# values. We therefore report and classify these cases distinctly:
# "steric-excluded (no accessible primary minimum found >= D_min)" rather
# than quoting the raw D_min-limited number as a real barrier.
D_MIN_ARTIFACT_NOTE = "steric-excluded: V(D) monotonically repulsive to D_min=0.3nm cutoff -- reported value is D_min-limited, not a true interior barrier"

for name, sc in SCENARIOS.items():
    vmax_col, barrier_present_col, note_col = [], [], []
    for _, r in df.iterrows():
        prof = dlvo.energy_profile_extended(
            radius_nm=r["radius_nm"], zeta_mV=r["zeta_mV"], ionic_strength_M=r["I_mM"] / 1000.0,
            hamaker_J=HAMAKER_J, layer_thickness_nm=sc["L_nm"], anchor_spacing_nm_val=sc["s_nm"],
        )
        barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
        vmax_col.append(barrier["V_max_kT"])
        barrier_present_col.append(barrier["barrier_present"])
        note_col.append(D_MIN_ARTIFACT_NOTE if (not barrier["barrier_present"] and barrier["V_max_kT"] > MARGINAL_KT) else "")
    df[f"V_extended_{name}_kT"] = vmax_col
    df[f"barrier_present_{name}"] = barrier_present_col
    if name == "mid":
        df["note_mid"] = note_col

df["stability_classic"] = df["V_classic_kT"].apply(lambda v: dlvo.classify_stability(v, STABLE_KT, MARGINAL_KT))
# For the extended model: classify "steric-excluded" cases as stable (repulsion dominates the
# entire resolvable range) but keep the raw number + note for transparency rather than pretending
# it's a normal barrier height.
df["stability_extended_mid"] = df.apply(
    lambda r: ("stable (steric-excluded)" if r["note_mid"] else dlvo.classify_stability(r["V_extended_mid_kT"], STABLE_KT, MARGINAL_KT)),
    axis=1,
)
# Capped value for plotting only -- keeps the chart readable; raw value stays in the CSV.
PLOT_CAP_KT = 40.0
df["V_extended_mid_kT_plot"] = df["V_extended_mid_kT"].clip(upper=PLOT_CAP_KT)
df["V_extended_low_kT_plot"] = df["V_extended_low_kT"].clip(upper=PLOT_CAP_KT)
df["V_extended_high_kT_plot"] = df["V_extended_high_kT"].clip(upper=PLOT_CAP_KT)

out_csv = f"{OUTDIR}/extended_dlvo_demo_results.csv"
df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")
print(df[["label", "V_classic_kT", "V_extended_low_kT", "V_extended_mid_kT", "V_extended_high_kT",
          "stability_classic", "stability_extended_mid", "note_mid"]].to_string(index=False))

# ---------------------------------------------------------------------------
# Comparison bar chart: classic vs extended (mid, with low-high as range)
# Bars are CAPPED at PLOT_CAP_KT for readability -- raw (often D_min-limited,
# very large) values are in the CSV and annotated on capped bars here.
# ---------------------------------------------------------------------------
plot_df = df.sort_values("V_extended_mid_kT").reset_index(drop=True)
y = np.arange(len(plot_df))
h = 0.35

fig, ax = plt.subplots(figsize=(9, 8))
colors_classic = ["#7f8c8d" for _ in plot_df["group"]]
colors_ext = ["#c0392b" if g == "measured grid" else "#2980b9" for g in plot_df["group"]]

ax.barh(y + h/2, plot_df["V_classic_kT"], height=h, color=colors_classic, label="classical DLVO (EDL+vdW only)")
ax.barh(y - h/2, plot_df["V_extended_mid_kT_plot"], height=h, color=colors_ext,
        label="extended DLVO (+ steric, mid estimate; capped at 40 kT for display)")

for i, r in plot_df.iterrows():
    if r["note_mid"]:
        ax.annotate("steric-excluded\n(D_min-limited)", xy=(PLOT_CAP_KT, i - h/2), xytext=(PLOT_CAP_KT + 1, i - h/2),
                    va="center", ha="left", fontsize=6.5, color="#c0392b" if r["group"] == "measured grid" else "#2980b9")

ax.set_yticks(y)
ax.set_yticklabels(plot_df["label"])
ax.axvline(STABLE_KT, color="green", linestyle="--", linewidth=1.2)
ax.axvline(MARGINAL_KT, color="orange", linestyle="--", linewidth=1.2)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlim(-5, PLOT_CAP_KT + 15)
ax.set_xlabel(r"Energy barrier $V_{max}$ ($k_BT$), capped at 40 kT for display")
ax.set_title("Classical vs. extended (+steric) DLVO -- 08/12/2026 data\n"
             "DEMO: steric term uses literature-typical L=10-20 nm, Gamma=1-3 mg/m^2 (not measured on this system)\n"
             "\"steric-excluded\" = V(D) never turns over before the 0.3 nm cutoff -- not a conventional barrier, see notes",
             fontsize=8.5)
legend_elems = [
    Patch(facecolor="#7f8c8d", label="classical DLVO (EDL + vdW)"),
    Patch(facecolor="#c0392b", label="extended, measured grid condition"),
    Patch(facecolor="#2980b9", label="extended, target food condition"),
    plt.Line2D([0], [0], color="green", linestyle="--", label="stable threshold (15 kT)"),
    plt.Line2D([0], [0], color="orange", linestyle="--", label="marginal threshold (10 kT)"),
]
ax.legend(handles=legend_elems, fontsize=7.5, loc="lower right")
fig.tight_layout()
fig.savefig(f"{OUTDIR}/classic_vs_extended_barchart.png", dpi=150)
print(f"Saved: {OUTDIR}/classic_vs_extended_barchart.png")

# ---------------------------------------------------------------------------
# Example V(D) curve overlay for one representative condition (whole milk)
# ---------------------------------------------------------------------------
wm = tg[tg["condition"] == "whole milk"].iloc[0]
fig2, ax2 = plt.subplots(figsize=(7, 5.5))
prof_classic = dlvo.energy_profile(radius_repr, wm["zeta_mV_predicted"], wm["ionic_strength_mM"] / 1000.0, HAMAKER_J)
ax2.plot(prof_classic["D_nm"], prof_classic["V_total_kT"], color="#7f8c8d", linewidth=2, label="classical DLVO (EDL+vdW)")
for name, sc in SCENARIOS.items():
    prof_ext = dlvo.energy_profile_extended(radius_repr, wm["zeta_mV_predicted"], wm["ionic_strength_mM"] / 1000.0,
                                             HAMAKER_J, sc["L_nm"], sc["s_nm"])
    ax2.plot(prof_ext["D_nm"], prof_ext["V_total_kT"], linewidth=1.6, label=f"extended, {name} (L={sc['L_nm']:.0f} nm)")
ax2.axhline(STABLE_KT, color="green", linestyle="--", linewidth=1)
ax2.axhline(MARGINAL_KT, color="orange", linestyle="--", linewidth=1)
ax2.axhline(0, color="black", linewidth=0.5)
ax2.set_xlim(0, 40)
ax2.set_ylim(-20, 60)
ax2.set_xlabel("Surface-to-surface separation D (nm)")
ax2.set_ylabel(r"Total interaction energy $V(D)$ ($k_BT$)")
ax2.set_title("Whole milk condition (pH 6.6, 80 mM): classical vs. extended DLVO\nDEMO -- literature steric parameters, not measured", fontsize=9.5)
ax2.legend(fontsize=8)
fig2.tight_layout()
fig2.savefig(f"{OUTDIR}/whole_milk_classic_vs_extended_profile.png", dpi=150)
print(f"Saved: {OUTDIR}/whole_milk_classic_vs_extended_profile.png")

print("\nDone.")
