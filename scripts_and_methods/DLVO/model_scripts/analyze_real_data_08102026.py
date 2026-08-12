"""
analyze_real_data_08102026.py
==============================

Runs the DLVO pipeline (dlvo_model.py / run_dlvo_analysis.py) on the newly
acquired real DLS dataset (DLS_data_zein_caseinate_08102026.csv), with one
deliberate deviation from the stock pipeline:

FLAGGED CONDITION: pH 4 / 150 mM NaCl (samples S019-S021)
-----------------------------------------------------------
All three replicates report zeta = 0.00 mV exactly and radius = 100000 nm
exactly. These are not plausible DLS readings (real particle radii across
every other condition are 90-370 nm; a "0.0 mV" zeta measured 3/3 times to
two decimal places, alongside a suspiciously round 100000 nm size, is the
signature of an instrument placeholder/off-scale value reported when a
sample has visibly flocculated/precipitated and the correlation function
can no longer be fit) -- NOT a real quantitative measurement of particle
size or charge.

Physically this is a plausible real event: pH 4 is far from the zein-
caseinate IEP (~5.8) so the particles should be highly charged, but 150 mM
is also the highest ionic strength tested, and double-layer screening at
high salt is exactly the condition expected to destabilize a
charge-stabilized colloid. So the *qualitative* result (this condition
aggregated) is likely real and consistent with DLVO physics -- but the
literal numbers (a=100000 nm, zeta=0 mV) cannot be fed into the
mechanistic model as if they were a normal anchor point:
  - Feeding zeta=0 into the sigmoid zeta(pH) fit for the 150 mM group
    distorts the fitted low-pH plateau for that ionic strength, which
    then contaminates the log-linear interpolation used for target
    conditions between 50 and 150 mM (e.g. whole milk at 80 mM).
  - Averaging radius_nm_mean across all 9 measured conditions to get the
    single representative radius used for the continuous stability sweep
    means one 100000 nm outlier would dominate that average and corrupt
    the ENTIRE stability map, not just the pH 4/150 mM cell.

Handling used here (flagged, not silently dropped):
  - EXCLUDED from: the sigmoidal zeta(pH) fit (pooled + per-I), and the
    representative radius used for the continuous pH x ionic-strength
    sweep.
  - INCLUDED in: condition_summary.csv (reported as its own row, with
    stability_call computed from the actual measured a/zeta -- which
    correctly comes out "unstable" since zeta=0 means no electrostatic
    barrier at all, just pure van der Waals attraction) and marked with
    a `flagged_aggregated=True` column.
  - Marked with an "X" annotation on the stability map at its (pH, I)
    coordinates, and left out of the energy-profile overlay plot (its own
    V(D) curve is degenerate/uninformative -- zero barrier at all D --
    and would just be a flat attractive line dominating the y-axis).

Everything else in the pipeline (Debye length, Grahame equation, EDL +
vdW energy, energy barrier extraction, stability thresholds, Hamaker
constant) is unchanged from dlvo_model.py / run_dlvo_analysis.py -- see
README.md in this folder for the full physics writeup.
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
sys.path.insert(0, "/sessions/youthful-gifted-dijkstra/mnt/dairy_protein_USDA/scripts_and_methods/DLVO")

import dlvo_model as dlvo
import run_dlvo_analysis as pipeline

HERE = "/sessions/youthful-gifted-dijkstra/mnt/dairy_protein_USDA/scripts_and_methods/DLVO"
REAL_CSV = os.path.join(HERE, "DLS_data_zein_caseinate_08102026.csv")
OUTDIR = os.path.join(HERE, "outputs", "real_data_08102026")
os.makedirs(OUTDIR, exist_ok=True)

HAMAKER_J = dlvo.HAMAKER_ZEIN_DEFAULT_J
STABLE_KT = pipeline.STABLE_THRESHOLD_KT
MARGINAL_KT = pipeline.MARGINAL_THRESHOLD_KT

# Radii above this are treated as DLS-instrument placeholder/off-scale
# values (aggregated sample), not real size measurements. Every genuine
# reading in this dataset is 90-370 nm.
AGGREGATION_RADIUS_CUTOFF_NM = 2000.0

print(f"Loading real DLS data: {REAL_CSV}")
raw = pipeline.load_anchor_data(REAL_CSV)
agg = pipeline.aggregate_replicates(raw)
agg["flagged_aggregated"] = agg["radius_nm_mean"] > AGGREGATION_RADIUS_CUTOFF_NM
print(f"  {len(raw)} raw rows -> {len(agg)} unique (pH, ionic strength) conditions")
flagged = agg[agg["flagged_aggregated"]]
for _, r in flagged.iterrows():
    print(f"  FLAGGED as aggregated/off-scale: pH={r['pH']}, I={r['ionic_strength_M']*1000:.0f} mM "
          f"(zeta={r['zeta_mV_mean']:.2f} mV, radius={r['radius_nm_mean']:.0f} nm) "
          "-- excluded from sigmoid fit and representative-radius calc, reported separately.")

agg_clean = agg[~agg["flagged_aggregated"]].reset_index(drop=True)

# ---------------------------------------------------------------------------
# zeta(pH, I) interpolation model
# ---------------------------------------------------------------------------
# The stock pipeline assumes a sigmoidal, isoelectric-point (IEP) crossing
# zeta(pH) curve (positive at low pH, negative at high pH, per the bare-zein
# behavior in Su et al. 2024). The REAL data does not follow that shape:
# zeta is negative at every measured pH (4, 5.5, 7) and every ionic
# strength -- there is no sign change / IEP in the measured range at all
# (likely because the caseinate coating, not bare zein, dominates surface
# charge here). Fitting the 4-parameter sigmoid to data with no zero
# crossing is unconstrained: it diverged to pH_iep=14.6, zeta_high_pH=-2394 mV
# (see OptimizeWarning above) -- physically meaningless and unsafe to
# interpolate/extrapolate with.
#
# Given only 3 pH anchor points per ionic-strength group anyway, a
# piecewise-linear interpolation in pH (per I group, clamped/flat outside
# the measured pH range -- i.e. no extrapolation) is the more honest,
# non-parametric choice for this dataset. Ionic strength is still
# interpolated log-linearly between measured groups, as in the stock
# pipeline.
print("zeta(pH) does not show a sigmoidal IEP crossing in the real data "
      "(negative at all measured pH 4-7) -- using piecewise-linear pH "
      "interpolation per ionic-strength group instead of the sigmoid fit.")

I_GROUPS = sorted(agg_clean["ionic_strength_M"].unique())
PH_BY_I = {}
ZETA_BY_I = {}
for I_val, g in agg_clean.groupby("ionic_strength_M"):
    g = g.sort_values("pH")
    PH_BY_I[I_val] = g["pH"].values
    ZETA_BY_I[I_val] = g["zeta_mV_mean"].values


def predict_zeta(pH: float, ionic_strength_M: float, fit=None) -> float:
    """Piecewise-linear in pH (per I group, clamped at range ends),
    log-linear in I between the nearest measured groups."""
    def zeta_at_I(I_val):
        return float(np.interp(pH, PH_BY_I[I_val], ZETA_BY_I[I_val]))

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


fit = None  # sigmoid fit not used for interpolation; kept as None sentinel

# ---------------------------------------------------------------------------
# Evaluate DLVO on every measured condition (including the flagged one,
# using its own real -- if degenerate -- a/zeta values)
# ---------------------------------------------------------------------------
print("Evaluating DLVO on measured conditions (all, including flagged)...")
rows = []
for _, r in agg.iterrows():
    prof = dlvo.energy_profile(
        radius_nm=r["radius_nm_mean"], zeta_mV=r["zeta_mV_mean"],
        ionic_strength_M=r["ionic_strength_M"], hamaker_J=HAMAKER_J,
    )
    barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
    stability = dlvo.classify_stability(barrier["V_max_kT"], STABLE_KT, MARGINAL_KT)
    rows.append({
        "pH": r["pH"], "ionic_strength_mM": r["ionic_strength_M"] * 1000.0,
        "zeta_mV_mean": r["zeta_mV_mean"], "zeta_mV_sd": r["zeta_mV_sd"],
        "radius_nm_mean": r["radius_nm_mean"], "radius_nm_sd": r["radius_nm_sd"],
        "n_replicates": r["n_replicates"],
        "kappa_inv_nm": prof["kappa_inv_nm"], "V_max_kT": barrier["V_max_kT"],
        "D_at_barrier_nm": barrier["D_max_nm"], "barrier_present": barrier["barrier_present"],
        "stability_call": stability, "flagged_aggregated": bool(r["flagged_aggregated"]),
    })
measured_results = pd.DataFrame(rows)
measured_out = os.path.join(OUTDIR, "condition_summary.csv")
measured_results.to_csv(measured_out, index=False)
print(f"  -> {measured_out}")
print(measured_results.to_string(index=False))

# ---------------------------------------------------------------------------
# Continuous stability sweep -- clean data only for fit + representative radius
# ---------------------------------------------------------------------------
print("Sweeping continuous pH x ionic-strength stability map (clean data)...")
radius_repr = agg_clean["radius_nm_mean"].mean()
print(f"  representative radius for sweep: {radius_repr:.1f} nm "
      f"(mean of {len(agg_clean)} clean conditions; flagged condition excluded)")
def sweep_stability_map(hamaker_J, radius_nm, pH_range=(4.0, 7.0), I_range_mM=(10, 150),
                         n_pH=60, n_I=60):
    """Same grid-sweep logic as pipeline.sweep_stability_map, but using the
    local piecewise-linear predict_zeta() and restricted to the actual
    measured pH/I envelope (4-7, 10-150 mM) -- no extrapolation beyond
    measured conditions, consistent with the clamped linear interpolation."""
    pH_grid = np.linspace(*pH_range, n_pH)
    I_grid_mM = np.geomspace(*I_range_mM, n_I)
    rows = []
    for I_mM in I_grid_mM:
        I_M = I_mM / 1000.0
        for pH in pH_grid:
            zeta = predict_zeta(pH, I_M)
            prof = dlvo.energy_profile(radius_nm, zeta, I_M, hamaker_J)
            barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
            rows.append({"pH": pH, "ionic_strength_mM": I_mM,
                         "zeta_mV_predicted": zeta, "V_max_kT": barrier["V_max_kT"]})
    return pd.DataFrame(rows)


sweep_df = sweep_stability_map(HAMAKER_J, radius_repr)
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
    target_rows.append({
        "condition": cond["label"], "pH": cond["pH"],
        "ionic_strength_mM": cond["ionic_strength_M"] * 1000.0,
        "zeta_mV_predicted": zeta, "V_max_kT": barrier["V_max_kT"],
        "stability_call": stability,
    })
target_df = pd.DataFrame(target_rows)
target_out = os.path.join(OUTDIR, "target_condition_predictions.csv")
target_df.to_csv(target_out, index=False)
print(f"  -> {target_out}")
print(target_df.to_string(index=False))

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot_energy_profiles(agg_clean, flagged, hamaker_J, out_path):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    cmap = plt.get_cmap("viridis")
    conditions = agg_clean.sort_values(["ionic_strength_M", "pH"]).reset_index(drop=True)
    n = len(conditions)
    for i, r in conditions.iterrows():
        prof = dlvo.energy_profile(r["radius_nm_mean"], r["zeta_mV_mean"], r["ionic_strength_M"], hamaker_J)
        label = f"pH {r['pH']:.1f}, {r['ionic_strength_M']*1000:.0f} mM"
        ax.plot(prof["D_nm"], prof["V_total_kT"], color=cmap(i / max(n - 1, 1)), label=label)

    ax.axhline(STABLE_KT, color="green", linestyle="--", linewidth=1, label=f"stable threshold ({STABLE_KT:.0f} kT)")
    ax.axhline(MARGINAL_KT, color="orange", linestyle="--", linewidth=1, label=f"marginal threshold ({MARGINAL_KT:.0f} kT)")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xlim(0, 20)
    ax.set_ylim(-30, 60)
    ax.set_xlabel("Surface-to-surface separation D (nm)")
    ax.set_ylabel(r"Total interaction energy $V(D)$ ($k_BT$)")
    title = "DLVO interaction energy profiles -- measured conditions (real DLS, 08/10/2026)"
    if len(flagged):
        title += "\n(pH 4 / 150 mM omitted -- aggregated sample, zeta=0 mV/off-scale size; see note)"
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_stability_map(sweep_df, flagged, out_path):
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
    ax.set_title("Predicted colloidal stability map -- real DLS data (08/10/2026)")

    for cond in pipeline.TARGET_CONDITIONS:
        ax.scatter(cond["pH"], cond["ionic_strength_M"] * 1000.0, marker="*", s=180,
                   color="blue", edgecolor="white", zorder=5)
        ax.annotate(cond["label"], (cond["pH"], cond["ionic_strength_M"] * 1000.0),
                    textcoords="offset points", xytext=(6, 6), fontsize=8, color="blue")

    for _, r in flagged.iterrows():
        ax.scatter(r["pH"], r["ionic_strength_M"] * 1000.0, marker="X", s=140,
                   color="red", edgecolor="white", zorder=6)
        ax.annotate("observed aggregation\n(off-scale DLS)", (r["pH"], r["ionic_strength_M"] * 1000.0),
                    textcoords="offset points", xytext=(8, -14), fontsize=7, color="red")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


print("Generating plots...")
plot_energy_profiles(agg_clean, flagged, HAMAKER_J, os.path.join(OUTDIR, "energy_profiles.png"))
plot_stability_map(sweep_df, flagged, os.path.join(OUTDIR, "stability_map.png"))
print(f"  -> {os.path.join(OUTDIR, 'energy_profiles.png')}")
print(f"  -> {os.path.join(OUTDIR, 'stability_map.png')}")

fit_params_out = os.path.join(OUTDIR, "zeta_interp_model.json")
with open(fit_params_out, "w") as f:
    json.dump({
        "interpolation_method": "piecewise-linear in pH per ionic-strength group "
            "(clamped at range ends), log-linear in ionic strength between groups. "
            "Sigmoidal IEP fit was attempted but diverged (no sign change in real "
            "zeta(pH) data over pH 4-7) and was discarded -- see script docstring.",
        "measured_anchor_points_mV": {
            f"{I*1000:.0f}mM": {"pH": PH_BY_I[I].tolist(), "zeta_mV": ZETA_BY_I[I].tolist()}
            for I in I_GROUPS
        },
        "hamaker_J": HAMAKER_J, "radius_nm_used_for_sweep": radius_repr,
        "flagged_conditions_excluded_from_fit": flagged[["pH", "ionic_strength_M"]].to_dict("records"),
    }, f, indent=2)
print(f"  -> {fit_params_out}")
print("\nDone.")
