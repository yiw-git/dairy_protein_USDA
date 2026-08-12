"""
analyze_real_data_08122026.py
===============================

DLVO analysis of the second real DLS dataset (08/12/2026): expanded pH
range (3, 4, 5.5, 7) x finer ionic-strength scan at the pH 4 instability
region (25, 50, 100 mM) x SINGLE replicate only (replicates 2-3 not yet
collected -- user will fill these in later).

Handling notes
--------------
1. SINGLE REPLICATE: rows for replicate 2/3 exist in the CSV with pH/I
   filled in but zeta_mV/radius_nm blank -- these are placeholders for
   data not yet collected, not missing/failed measurements. They are
   dropped before aggregation (not treated as an error). Every reported
   number here is therefore a single measurement, not a replicate mean --
   NO error bars / SD are available yet, and single noisy points cannot
   be distinguished from real trends. Flagged clearly in all outputs.

2. pH 4 / 100 mM: zeta was measured fine (9.05 mV) but radius_nm=8404 nm,
   ~20-100x every other condition in this dataset (all others are
   60-250 nm). Unlike the previous dataset's pH4/150mM sentinel
   (zeta=0, radius=100000 -- an instrument placeholder for a fully failed
   read), this looks like a genuine DLS z-average reading dominated by
   large aggregates/clusters -- i.e. real evidence this condition is
   close to the critical coagulation concentration (CCC) at pH 4, but
   the z-average size is skewed by a broad/bimodal distribution and is
   not representative of "the particle size" in the same sense as the
   other conditions. It is INCLUDED in condition_summary (with its own,
   real barrier calculation) but EXCLUDED from the representative radius
   used for the continuous sweep (one 8404 nm value would otherwise
   dominate that average).

3. Sigmoid zeta(pH) fit: this dataset (unlike 08/10) actually shows a
   sign crossover -- positive zeta at pH 3 and (weakly) pH 4, negative at
   pH 5.5 and 7 -- consistent with a true isoelectric point somewhere
   around pH 4-5. The 4-parameter sigmoid is fit here and checked for
   sanity (bounded, IEP inside the measured range) before being trusted;
   falls back to piecewise-linear interpolation if it does not converge
   sensibly.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
HERE = "/sessions/youthful-gifted-dijkstra/mnt/dairy_protein_USDA/scripts_and_methods/DLVO"
sys.path.insert(0, os.path.join(HERE, "model_scripts"))
import dlvo_model as dlvo
import run_dlvo_analysis as pipeline

CSV_PATH = os.path.join(HERE, "data_DLVO", "DLS_data_zein_caseinate_08122026.csv")
OUTDIR = os.path.join(HERE, "outputs", "real_data_08122026")
os.makedirs(OUTDIR, exist_ok=True)

HAMAKER_J = dlvo.HAMAKER_ZEIN_DEFAULT_J
STABLE_KT = pipeline.STABLE_THRESHOLD_KT
MARGINAL_KT = pipeline.MARGINAL_THRESHOLD_KT
RADIUS_OUTLIER_CUTOFF_NM = 2000.0

# ---------------------------------------------------------------------------
# Load -- drop unfilled replicate placeholder rows and fully-blank rows
# ---------------------------------------------------------------------------
print(f"Loading: {CSV_PATH}")
raw = pd.read_csv(CSV_PATH)
raw = raw.dropna(subset=["pH"])  # drop trailing fully-blank rows
n_total = len(raw)
raw_filled = raw.dropna(subset=["zeta_mV", "radius_nm"]).copy()
n_placeholder = n_total - len(raw_filled)
print(f"  {n_total} rows in CSV; {n_placeholder} are unfilled replicate-2/3 "
      f"placeholders (dropped); {len(raw_filled)} rows with real data (all replicate 1).")
raw_filled["ionic_strength_M"] = raw_filled["ionic_strength_mM"] / 1000.0

agg = (
    raw_filled.groupby(["pH", "ionic_strength_M"])
    .agg(zeta_mV_mean=("zeta_mV", "mean"), zeta_mV_sd=("zeta_mV", "std"),
         radius_nm_mean=("radius_nm", "mean"), radius_nm_sd=("radius_nm", "std"),
         n_replicates=("zeta_mV", "count"))
    .reset_index().sort_values(["pH", "ionic_strength_M"])
)
agg["radius_outlier"] = agg["radius_nm_mean"] > RADIUS_OUTLIER_CUTOFF_NM
print(f"  {len(agg)} unique (pH, I) conditions, n_replicates={agg['n_replicates'].iloc[0]} each "
      "(single replicate -- no SD available).")
outlier_rows = agg[agg["radius_outlier"]]
for _, r in outlier_rows.iterrows():
    print(f"  RADIUS OUTLIER (real zeta, but z-avg radius {r['radius_nm_mean']:.0f} nm dominated by "
          f"large aggregates): pH={r['pH']}, I={r['ionic_strength_M']*1000:.0f} mM -- kept in "
          "condition_summary with its own real barrier calc, excluded from representative-radius avg.")

agg_radius_clean = agg[~agg["radius_outlier"]]

# ---------------------------------------------------------------------------
# zeta(pH) model: try the real sigmoid now that there's a sign crossover
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
        "zeta_mV": r["zeta_mV_mean"], "radius_nm": r["radius_nm_mean"],
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
      "conditions excluding the pH4/100mM radius outlier)")

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
# Target food conditions (clamped to measured pH/I envelope)
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
    ax.set_title("DLVO energy profiles -- real DLS, 08/12/2026 (single replicate, preliminary)", fontsize=10)
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
    ax.set_title("Predicted colloidal stability map -- real DLS, 08/12/2026 (single replicate, preliminary)", fontsize=10)
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
        "note": "SINGLE REPLICATE dataset -- every number here is one measurement, not a mean. "
                "No SD/error bars available.",
        "n_conditions": int(len(agg)), "n_placeholder_rows_dropped": int(n_placeholder),
        "sigmoid_fit_used": bool(sigmoid_sane),
        "pooled_sigmoid": {"zeta_low_pH_mV": zeta_low_pH, "zeta_high_pH_mV": zeta_high_pH,
                            "pH_iep": pH_iep, "slope": slope} if sigmoid_sane else None,
        "hamaker_J": HAMAKER_J, "radius_nm_used_for_sweep": radius_repr,
        "radius_outlier_conditions": outlier_rows[["pH", "ionic_strength_M", "radius_nm_mean", "zeta_mV_mean"]].to_dict("records"),
    }, f, indent=2)
print(f"  -> {model_out}")
print("\nDone.")
