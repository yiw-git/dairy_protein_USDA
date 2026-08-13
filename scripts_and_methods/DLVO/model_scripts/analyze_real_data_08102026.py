"""
analyze_real_data_08102026.py
==============================

Runs the DLVO pipeline (dlvo_model.py / run_dlvo_analysis.py) on the
zein-caseinate DLS dataset DLS_data_zein_caseinate_08102026.csv
(pH 4/5.5/7 x ionic strength 10/50/150 mM, 3 replicates per condition).

UPDATED (this revision): the user re-measured/corrected this CSV. The
pH 4 / 150 mM condition (S019-S021) previously reported a full
instrument-failure sentinel (zeta=0.00 mV exactly, radius=100000 nm
exactly, all 3 replicates) -- treated in the prior revision of this
script as "unmeasurable, exclude entirely from both the zeta(pH) fit and
the representative radius." The new data for this condition instead
reports REAL, small positive zeta values (3.21-6.58 mV) with a large but
non-sentinel radius (8000 nm exactly across all 3 reps -- still ~20-90x
every other condition in this dataset, 90-370 nm elsewhere). This matches
the pattern seen in the 08/11 and 08/12 trials (pH4/high-I: measurable
zeta, real DLS z-average dominated by large aggregates) -- i.e. genuine
evidence of aggregation near the CCC at pH 4, not an unmeasurable sample.

Handling (updated to match the 08/11/08/12 scripts' approach):
  - radius_outlier flag (radius_nm_mean > 2000 nm): EXCLUDED only from
    the representative radius used in the continuous sweep (one 8000 nm
    value would otherwise dominate that average). Its own zeta IS used
    in the zeta(pH, I) model -- unlike the old sentinel row, this zeta is
    a real measurement.
  - INCLUDED in condition_summary.csv with its own real barrier
    calculation and a `radius_outlier=True` flag.
  - Sigmoid zeta(pH) fit is attempted (as in the 08/11/08/12 scripts) and
    checked for sanity (IEP inside the measured pH range, bounded
    plateaus) before trusting it; falls back to piecewise-linear pH
    interpolation per ionic-strength group if it does not converge
    sensibly. (The previous revision of this script skipped straight to
    piecewise-linear because the old all-negative zeta data had no sign
    crossover at all; the updated pH4/150mM point reintroduces a
    crossover, so the sigmoid is retried here.)

Everything else (Debye length, Grahame equation, EDL + vdW energy,
barrier extraction, stability thresholds, Hamaker constant) is unchanged
from dlvo_model.py / run_dlvo_analysis.py -- see README.md for the full
physics writeup.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/sessions/youthful-gifted-dijkstra/mnt/dairy_protein_USDA/scripts_and_methods/DLVO"
sys.path.insert(0, os.path.join(HERE, "model_scripts"))
import dlvo_model as dlvo
import run_dlvo_analysis as pipeline

REAL_CSV = os.path.join(HERE, "data_DLVO", "DLS_data_zein_caseinate_08102026.csv")
OUTDIR = os.path.join(HERE, "outputs", "real_data_08102026")
os.makedirs(OUTDIR, exist_ok=True)

HAMAKER_J = dlvo.HAMAKER_ZEIN_DEFAULT_J
STABLE_KT = pipeline.STABLE_THRESHOLD_KT
MARGINAL_KT = pipeline.MARGINAL_THRESHOLD_KT
RADIUS_OUTLIER_CUTOFF_NM = 2000.0  # genuine readings here are 90-370 nm

# ---------------------------------------------------------------------------
# Load + aggregate (3 real replicates per condition)
# ---------------------------------------------------------------------------
print(f"Loading real DLS data: {REAL_CSV}")
raw = pipeline.load_anchor_data(REAL_CSV)
agg = pipeline.aggregate_replicates(raw)
agg["radius_outlier"] = agg["radius_nm_mean"] > RADIUS_OUTLIER_CUTOFF_NM
print(f"  {len(raw)} raw rows -> {len(agg)} unique (pH, ionic strength) conditions, "
      f"n_replicates={agg['n_replicates'].iloc[0]} each")
outlier_rows = agg[agg["radius_outlier"]]
for _, r in outlier_rows.iterrows():
    print(f"  RADIUS OUTLIER (real zeta, large aggregates): pH={r['pH']}, I={r['ionic_strength_M']*1000:.0f} mM "
          f"(zeta={r['zeta_mV_mean']:.2f} +/- {r['zeta_mV_sd']:.2f} mV, radius={r['radius_nm_mean']:.0f} nm) "
          "-- kept in condition_summary + zeta(pH) model, excluded only from representative-radius calc.")

agg_radius_clean = agg[~agg["radius_outlier"]]

# ---------------------------------------------------------------------------
# zeta(pH, I) model: try the sigmoid, sanity-check it, fall back if needed
# ---------------------------------------------------------------------------
pooled_popt, pooled_pcov = dlvo.fit_zeta_vs_pH(agg["pH"].values, agg["zeta_mV_mean"].values)
zeta_low_pH, zeta_high_pH, pH_iep, slope = pooled_popt
pH_min, pH_max = agg["pH"].min(), agg["pH"].max()
sigmoid_sane = (
    pH_min <= pH_iep <= pH_max
    and abs(zeta_low_pH) < 200 and abs(zeta_high_pH) < 200
    and np.all(np.isfinite(pooled_pcov))
)
print(f"Pooled sigmoid fit: zeta_low_pH={zeta_low_pH:.1f} mV, zeta_high_pH={zeta_high_pH:.1f} mV, "
      f"pH_iep={pH_iep:.2f}, slope={slope:.2f} -- {'SANE, using it' if sigmoid_sane else 'discarded, falling back to linear interpolation'}")

I_GROUPS = sorted(agg["ionic_strength_M"].unique())
PH_BY_I, ZETA_BY_I = {}, {}
for I_val, g in agg.groupby("ionic_strength_M"):
    g = g.sort_values("pH")
    PH_BY_I[I_val] = g["pH"].values
    ZETA_BY_I[I_val] = g["zeta_mV_mean"].values

per_I_sigmoid = {}
if sigmoid_sane:
    for I_val in I_GROUPS:
        try:
            popt_i, _ = dlvo.fit_zeta_vs_pH_fixed_shape(PH_BY_I[I_val], ZETA_BY_I[I_val], pH_iep, slope)
            per_I_sigmoid[I_val] = tuple(popt_i)
        except Exception as e:
            print(f"  per-I fixed-shape fit failed for I={I_val*1000:.0f} mM ({e}); using pooled plateaus.")
            per_I_sigmoid[I_val] = (zeta_low_pH, zeta_high_pH)


def predict_zeta(pH: float, ionic_strength_M: float) -> float:
    def zeta_at_I_linear(I_val):
        return float(np.interp(pH, PH_BY_I[I_val], ZETA_BY_I[I_val]))

    def zeta_at_I_sigmoid(I_val):
        zlo, zhi = per_I_sigmoid[I_val]
        return float(dlvo.sigmoid_zeta_ph(pH, zlo, zhi, pH_iep, slope))

    zeta_at_I = zeta_at_I_sigmoid if sigmoid_sane else zeta_at_I_linear

    if ionic_strength_M <= I_GROUPS[0]:
        return zeta_at_I(I_GROUPS[0])
    elif ionic_strength_M >= I_GROUPS[-1]:
        return zeta_at_I(I_GROUPS[-1])
    else:
        lo = max(v for v in I_GROUPS if v <= ionic_strength_M)
        hi = min(v for v in I_GROUPS if v >= ionic_strength_M)
        if lo == hi:
            return zeta_at_I(lo)
        t = (np.log(ionic_strength_M) - np.log(lo)) / (np.log(hi) - np.log(lo))
        return zeta_at_I(lo) + t * (zeta_at_I(hi) - zeta_at_I(lo))


# ---------------------------------------------------------------------------
# Evaluate DLVO on every measured condition
# ---------------------------------------------------------------------------
print("Evaluating DLVO on measured conditions...")
rows = []
for _, r in agg.iterrows():
    prof = dlvo.energy_profile(r["radius_nm_mean"], r["zeta_mV_mean"], r["ionic_strength_M"], HAMAKER_J)
    barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
    stability = dlvo.classify_stability(barrier["V_max_kT"], STABLE_KT, MARGINAL_KT)
    rows.append({
        "pH": r["pH"], "ionic_strength_mM": r["ionic_strength_M"] * 1000.0,
        "zeta_mV_mean": r["zeta_mV_mean"], "zeta_mV_sd": r["zeta_mV_sd"],
        "radius_nm_mean": r["radius_nm_mean"], "radius_nm_sd": r["radius_nm_sd"],
        "n_replicates": r["n_replicates"], "kappa_inv_nm": prof["kappa_inv_nm"],
        "V_max_kT": barrier["V_max_kT"], "D_at_barrier_nm": barrier["D_max_nm"],
        "barrier_present": barrier["barrier_present"], "stability_call": stability,
        "radius_outlier": bool(r["radius_outlier"]),
    })
measured_results = pd.DataFrame(rows)
measured_out = os.path.join(OUTDIR, "condition_summary.csv")
measured_results.to_csv(measured_out, index=False)
print(f"  -> {measured_out}")
print(measured_results.to_string(index=False))

# ---------------------------------------------------------------------------
# Continuous sweep -- radius-clean data only for representative radius
# ---------------------------------------------------------------------------
radius_repr = agg_radius_clean["radius_nm_mean"].mean()
print(f"Representative radius for sweep: {radius_repr:.1f} nm (mean of {len(agg_radius_clean)} "
      "conditions excluding the pH4/150mM radius outlier)")

pH_range = (float(agg["pH"].min()), float(agg["pH"].max()))
I_range_mM = (float(agg["ionic_strength_M"].min() * 1000), float(agg["ionic_strength_M"].max() * 1000))


def sweep_stability_map(hamaker_J, radius_nm, pH_range, I_range_mM, n_pH=60, n_I=60):
    pH_grid = np.linspace(*pH_range, n_pH)
    I_grid_mM = np.geomspace(*I_range_mM, n_I)
    rows = []
    for I_mM in I_grid_mM:
        I_M = I_mM / 1000.0
        for pH in pH_grid:
            zeta = predict_zeta(pH, I_M)
            prof = dlvo.energy_profile(radius_nm, zeta, I_M, hamaker_J)
            barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
            rows.append({"pH": pH, "ionic_strength_mM": I_mM, "zeta_mV_predicted": zeta,
                         "V_max_kT": barrier["V_max_kT"]})
    return pd.DataFrame(rows)


print("Sweeping continuous pH x ionic-strength stability map...")
sweep_df = sweep_stability_map(HAMAKER_J, radius_repr, pH_range, I_range_mM)
sweep_out = os.path.join(OUTDIR, "stability_sweep_grid.csv")
sweep_df.to_csv(sweep_out, index=False)
print(f"  -> {sweep_out}")

# ---------------------------------------------------------------------------
# Target food conditions
# ---------------------------------------------------------------------------
print("Evaluating target conditions (whole milk / yogurt / bread-release)...")
target_rows = []
for cond in pipeline.TARGET_CONDITIONS:
    zeta = predict_zeta(cond["pH"], cond["ionic_strength_M"])
    prof = dlvo.energy_profile(radius_repr, zeta, cond["ionic_strength_M"], HAMAKER_J)
    barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
    stability = dlvo.classify_stability(barrier["V_max_kT"], STABLE_KT, MARGINAL_KT)
    target_rows.append({"condition": cond["label"], "pH": cond["pH"],
                        "ionic_strength_mM": cond["ionic_strength_M"] * 1000.0,
                        "zeta_mV_predicted": zeta, "V_max_kT": barrier["V_max_kT"],
                        "stability_call": stability})
target_df = pd.DataFrame(target_rows)
target_out = os.path.join(OUTDIR, "target_condition_predictions.csv")
target_df.to_csv(target_out, index=False)
print(f"  -> {target_out}")
print(target_df.to_string(index=False))

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot_energy_profiles(agg, hamaker_J, out_path):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    cmap = plt.get_cmap("viridis")
    conditions = agg.sort_values(["ionic_strength_M", "pH"]).reset_index(drop=True)
    n = len(conditions)
    for i, r in conditions.iterrows():
        prof = dlvo.energy_profile(r["radius_nm_mean"], r["zeta_mV_mean"], r["ionic_strength_M"], hamaker_J)
        label = f"pH {r['pH']:.1f}, {r['ionic_strength_M']*1000:.0f} mM"
        if r["radius_outlier"]:
            label += " (radius outlier)"
        ls = "--" if r["radius_outlier"] else "-"
        ax.plot(prof["D_nm"], prof["V_total_kT"], color=cmap(i / max(n - 1, 1)), label=label, linestyle=ls)
    ax.axhline(STABLE_KT, color="green", linestyle="--", linewidth=1, label=f"stable threshold ({STABLE_KT:.0f} kT)")
    ax.axhline(MARGINAL_KT, color="orange", linestyle="--", linewidth=1, label=f"marginal threshold ({MARGINAL_KT:.0f} kT)")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xlim(0, 20)
    ax.set_ylim(-30, 60)
    ax.set_xlabel("Surface-to-surface separation D (nm)")
    ax.set_ylabel(r"Total interaction energy $V(D)$ ($k_BT$)")
    ax.set_title("DLVO energy profiles -- real DLS, 08/10/2026 (updated data)", fontsize=10)
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_stability_map(sweep_df, agg, out_path):
    pivot = sweep_df.pivot(index="ionic_strength_mM", columns="pH", values="V_max_kT")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    vmax = max(30.0, np.nanpercentile(pivot.values, 99))
    im = ax.pcolormesh(pivot.columns.values, pivot.index.values, pivot.values,
                        shading="auto", cmap="RdYlGn", vmin=-10, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"Energy barrier $V_{max}$ ($k_BT$)")
    cs = ax.contour(pivot.columns.values, pivot.index.values, pivot.values,
                     levels=[MARGINAL_KT, STABLE_KT], colors=["orange", "black"], linewidths=1.2)
    ax.clabel(cs, inline=True, fontsize=8, fmt="%.0f kT")
    ax.set_yscale("log")
    ax.set_xlabel("pH")
    ax.set_ylabel("Ionic strength (mM, log scale)")
    ax.set_title("Predicted colloidal stability map -- real DLS, 08/10/2026 (updated data)", fontsize=10)
    for cond in pipeline.TARGET_CONDITIONS:
        ax.scatter(cond["pH"], cond["ionic_strength_M"] * 1000.0, marker="*", s=180,
                   color="blue", edgecolor="white", zorder=5)
        ax.annotate(cond["label"], (cond["pH"], cond["ionic_strength_M"] * 1000.0),
                    textcoords="offset points", xytext=(6, 6), fontsize=8, color="blue")
    for _, r in agg[agg["radius_outlier"]].iterrows():
        ax.scatter(r["pH"], r["ionic_strength_M"] * 1000.0, marker="X", s=140,
                   color="red", edgecolor="white", zorder=6)
        ax.annotate("large aggregates\n(near CCC)", (r["pH"], r["ionic_strength_M"] * 1000.0),
                    textcoords="offset points", xytext=(8, -14), fontsize=7, color="red")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


print("Generating plots...")
plot_energy_profiles(agg, HAMAKER_J, os.path.join(OUTDIR, "energy_profiles.png"))
plot_stability_map(sweep_df, agg, os.path.join(OUTDIR, "stability_map.png"))
print(f"  -> {os.path.join(OUTDIR, 'energy_profiles.png')}")
print(f"  -> {os.path.join(OUTDIR, 'stability_map.png')}")

model_out = os.path.join(OUTDIR, "zeta_model_info.json")
with open(model_out, "w") as f:
    json.dump({
        "n_conditions": int(len(agg)),
        "sigmoid_fit_used": bool(sigmoid_sane),
        "pooled_sigmoid": {"zeta_low_pH_mV": zeta_low_pH, "zeta_high_pH_mV": zeta_high_pH,
                            "pH_iep": pH_iep, "slope": slope} if sigmoid_sane else None,
        "hamaker_J": HAMAKER_J, "radius_nm_used_for_sweep": radius_repr,
        "radius_outlier_conditions": outlier_rows[["pH", "ionic_strength_M", "radius_nm_mean", "zeta_mV_mean"]].to_dict("records"),
    }, f, indent=2)
print(f"  -> {model_out}")
print("\nDone.")
