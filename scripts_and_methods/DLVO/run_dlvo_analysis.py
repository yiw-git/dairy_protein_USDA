"""
run_dlvo_analysis.py
=====================

End-to-end DLVO analysis pipeline for Section 3.1.3 (zein-caseinate NP,
+/- WPI / casein exposure). Implements the full workflow described in
Section 3A of the memo:

    1. Load the measured DLS anchor matrix (pH x ionic strength grid,
       zeta potential + particle radius, with replicates).
    2. Average replicates -> condition-level table.
    3. Fit the sigmoidal zeta(pH) curve:
         (a) pooled across all ionic strengths, to pin down the shared
             shape parameters (pH_iep, slope);
         (b) per ionic-strength group (zeta_low_pH, zeta_high_pH only), to
             capture double-layer-screening magnitude changes with I.
    4. For every MEASURED condition: compute the full V(D) DLVO profile,
       extract the energy barrier V_max, and classify stability.
    5. Sweep DLVO across a CONTINUOUS pH x ionic-strength grid (using the
       fitted sigmoid to interpolate/extrapolate zeta), producing a
       stability map -- including conditions you cannot test directly
       right now (e.g. skim milk, or the specific dairy/food conditions
       used by the curcumin ML manuscript: whole milk, yogurt, and the
       bread-release condition).
    6. Save all outputs (CSVs + plots) to outputs/.

USAGE
-----
    # Run on the bundled synthetic demo data (no real DLS data needed --
    # good for testing the pipeline / sanity-checking equations before
    # your wet-lab window):
    python run_dlvo_analysis.py --demo

    # Run on your real measured anchor matrix once you have it:
    python run_dlvo_analysis.py --anchor-csv anchor_data_template.csv

See README.md for the full methodology writeup, equation references, and
a description of every output file.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless-safe backend
import matplotlib.pyplot as plt

import dlvo_model as dlvo

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(HERE, "outputs")
DEMO_CSV = os.path.join(HERE, "demo_data_zein_caseinate.csv")

# Hamaker constant to use for the sweep. Default: Su et al. 2024's fitted
# value for zein NPs specifically. Override with --hamaker if you have a
# system-specific estimate (e.g. after fitting your own CCC data the same
# way Su et al. did -- see README "Extending this script").
DEFAULT_HAMAKER_J = dlvo.HAMAKER_ZEIN_DEFAULT_J

# Specific food-relevant conditions from Section 4 of the memo that we
# cannot test directly right now but want DLVO predictions for. These are
# overlaid as annotated points on the stability map.
TARGET_CONDITIONS = [
    {"label": "whole milk", "pH": 6.6, "ionic_strength_M": 0.080},
    {"label": "yogurt", "pH": 4.5, "ionic_strength_M": 0.020},
    {"label": "bread-release", "pH": 5.3, "ionic_strength_M": 0.050},
]

STABLE_THRESHOLD_KT = 15.0
MARGINAL_THRESHOLD_KT = 10.0


# ---------------------------------------------------------------------------
# Step 1-2: load + aggregate anchor data
# ---------------------------------------------------------------------------
def load_anchor_data(csv_path: str) -> pd.DataFrame:
    """
    Expected columns (see anchor_data_template.csv):
        sample_id, protein_system, pH, ionic_strength_mM, zeta_mV,
        radius_nm, temperature_C, replicate
    """
    df = pd.read_csv(csv_path)
    required = {"pH", "ionic_strength_mM", "zeta_mV", "radius_nm"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Anchor CSV is missing required columns: {missing}")

    n_missing_zeta = df["zeta_mV"].isna().sum()
    n_missing_radius = df["radius_nm"].isna().sum()
    if n_missing_zeta or n_missing_radius:
        raise ValueError(
            f"Anchor CSV has {n_missing_zeta} blank zeta_mV and {n_missing_radius} "
            "blank radius_nm cells. This looks like the unfilled "
            "anchor_data_template.csv -- fill in your DLS measurements for every "
            "row (or delete rows you did not run) before analyzing. Use --demo to "
            "test the pipeline with synthetic data in the meantime."
        )

    df["ionic_strength_M"] = df["ionic_strength_mM"] / 1000.0
    return df


def aggregate_replicates(df: pd.DataFrame) -> pd.DataFrame:
    """Average replicates within each (pH, ionic_strength) condition."""
    group_cols = ["pH", "ionic_strength_M"]
    agg = (
        df.groupby(group_cols)
        .agg(
            zeta_mV_mean=("zeta_mV", "mean"),
            zeta_mV_sd=("zeta_mV", "std"),
            radius_nm_mean=("radius_nm", "mean"),
            radius_nm_sd=("radius_nm", "std"),
            n_replicates=("zeta_mV", "count"),
        )
        .reset_index()
        .sort_values(group_cols)
    )
    return agg


# ---------------------------------------------------------------------------
# Step 3: sigmoid zeta(pH) fitting -- pooled shape + per-I magnitude
# ---------------------------------------------------------------------------
def fit_zeta_ph_model(agg: pd.DataFrame) -> dict:
    """
    Returns a dict:
        {
          "pooled": (zeta_low_pH, zeta_high_pH, pH_iep, slope),
          "per_I": {ionic_strength_M: (zeta_low_pH, zeta_high_pH), ...},
          "I_values_sorted": [...],
        }
    per_I uses the pooled (pH_iep, slope) held fixed -- see
    dlvo_model.fit_zeta_vs_pH_fixed_shape docstring for why.
    """
    pooled_popt, _ = dlvo.fit_zeta_vs_pH(agg["pH"].values, agg["zeta_mV_mean"].values)
    pH_iep, slope = pooled_popt[2], pooled_popt[3]

    per_I = {}
    for I_val, group in agg.groupby("ionic_strength_M"):
        if len(group) < 2:
            # Not enough points even for the 2-parameter fixed-shape fit;
            # fall back to the pooled zeta_low_pH/zeta_high_pH.
            per_I[I_val] = (pooled_popt[0], pooled_popt[1])
            continue
        popt_i, _ = dlvo.fit_zeta_vs_pH_fixed_shape(
            group["pH"].values, group["zeta_mV_mean"].values, pH_iep, slope
        )
        per_I[I_val] = tuple(popt_i)

    return {
        "pooled": tuple(pooled_popt),
        "per_I": per_I,
        "I_values_sorted": sorted(per_I.keys()),
    }


def predict_zeta(pH: float, ionic_strength_M: float, fit: dict) -> float:
    """
    Predict zeta(pH, I) from the fitted model:
      - shape (pH_iep, slope) is shared (pooled fit)
      - magnitude (zeta_low_pH, zeta_high_pH) is interpolated log-linearly in I
        between the nearest measured ionic-strength groups (extrapolated
        with the nearest group's parameters if outside the measured range
        -- flagged in the README as an assumption to revisit with more
        data).
    """
    _, _, pH_iep, slope = fit["pooled"]
    I_values = fit["I_values_sorted"]

    if ionic_strength_M <= I_values[0]:
        zeta_low_pH, zeta_high_pH = fit["per_I"][I_values[0]]
    elif ionic_strength_M >= I_values[-1]:
        zeta_low_pH, zeta_high_pH = fit["per_I"][I_values[-1]]
    else:
        # log-linear interpolation between bracketing measured I groups
        lo = max(v for v in I_values if v <= ionic_strength_M)
        hi = min(v for v in I_values if v >= ionic_strength_M)
        if lo == hi:
            zeta_low_pH, zeta_high_pH = fit["per_I"][lo]
        else:
            t = (np.log(ionic_strength_M) - np.log(lo)) / (np.log(hi) - np.log(lo))
            zlow_lo, zhigh_lo = fit["per_I"][lo]
            zlow_hi, zhigh_hi = fit["per_I"][hi]
            zeta_low_pH = zlow_lo + t * (zlow_hi - zlow_lo)
            zeta_high_pH = zhigh_lo + t * (zhigh_hi - zhigh_lo)

    return float(dlvo.sigmoid_zeta_ph(pH, zeta_low_pH, zeta_high_pH, pH_iep, slope))


# ---------------------------------------------------------------------------
# Step 4: DLVO on every measured condition
# ---------------------------------------------------------------------------
def evaluate_measured_conditions(agg: pd.DataFrame, hamaker_J: float) -> pd.DataFrame:
    rows = []
    for _, r in agg.iterrows():
        prof = dlvo.energy_profile(
            radius_nm=r["radius_nm_mean"],
            zeta_mV=r["zeta_mV_mean"],
            ionic_strength_M=r["ionic_strength_M"],
            hamaker_J=hamaker_J,
        )
        barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
        stability = dlvo.classify_stability(
            barrier["V_max_kT"], STABLE_THRESHOLD_KT, MARGINAL_THRESHOLD_KT
        )
        rows.append({
            "pH": r["pH"],
            "ionic_strength_mM": r["ionic_strength_M"] * 1000.0,
            "zeta_mV_mean": r["zeta_mV_mean"],
            "radius_nm_mean": r["radius_nm_mean"],
            "kappa_inv_nm": prof["kappa_inv_nm"],
            "V_max_kT": barrier["V_max_kT"],
            "D_at_barrier_nm": barrier["D_max_nm"],
            "barrier_present": barrier["barrier_present"],
            "stability_call": stability,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 5: continuous pH x ionic-strength stability sweep
# ---------------------------------------------------------------------------
def sweep_stability_map(fit: dict, hamaker_J: float, radius_nm: float,
                         pH_range=(3.5, 7.5), I_range_mM=(5, 200),
                         n_pH=60, n_I=60) -> pd.DataFrame:
    """
    Sweep DLVO across a continuous pH x ionic-strength grid using the
    fitted sigmoid to supply zeta(pH, I) at every grid point. radius_nm
    is held at a single representative value (the mean measured radius)
    since particle size varies far less than zeta/ionic strength across
    conditions in the designed grid -- revisit if your DLS matrix shows
    otherwise.
    """
    pH_grid = np.linspace(*pH_range, n_pH)
    I_grid_mM = np.geomspace(*I_range_mM, n_I)  # log-spaced, matches how ionic strength is usually varied

    rows = []
    for I_mM in I_grid_mM:
        I_M = I_mM / 1000.0
        for pH in pH_grid:
            zeta = predict_zeta(pH, I_M, fit)
            prof = dlvo.energy_profile(radius_nm, zeta, I_M, hamaker_J)
            barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
            rows.append({
                "pH": pH,
                "ionic_strength_mM": I_mM,
                "zeta_mV_predicted": zeta,
                "V_max_kT": barrier["V_max_kT"],
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_energy_profiles(agg: pd.DataFrame, hamaker_J: float, out_path: str):
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("viridis")
    conditions = agg.sort_values(["ionic_strength_M", "pH"]).reset_index(drop=True)
    n = len(conditions)
    for i, r in conditions.iterrows():
        prof = dlvo.energy_profile(
            r["radius_nm_mean"], r["zeta_mV_mean"], r["ionic_strength_M"], hamaker_J
        )
        label = f"pH {r['pH']:.1f}, {r['ionic_strength_M']*1000:.0f} mM"
        ax.plot(prof["D_nm"], prof["V_total_kT"], color=cmap(i / max(n - 1, 1)), label=label)

    ax.axhline(STABLE_THRESHOLD_KT, color="green", linestyle="--", linewidth=1,
               label=f"stable threshold ({STABLE_THRESHOLD_KT:.0f} kT)")
    ax.axhline(MARGINAL_THRESHOLD_KT, color="orange", linestyle="--", linewidth=1,
               label=f"marginal threshold ({MARGINAL_THRESHOLD_KT:.0f} kT)")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xlim(0, 20)
    ax.set_ylim(-30, 60)
    ax.set_xlabel("Surface-to-surface separation D (nm)")
    ax.set_ylabel(r"Total interaction energy $V(D)$ ($k_BT$)")
    ax.set_title("DLVO interaction energy profiles -- measured conditions")
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_stability_map(sweep_df: pd.DataFrame, out_path: str):
    pivot = sweep_df.pivot(index="ionic_strength_mM", columns="pH", values="V_max_kT")

    fig, ax = plt.subplots(figsize=(8, 5.5))
    vmax = max(30.0, np.nanpercentile(pivot.values, 99))
    im = ax.pcolormesh(
        pivot.columns.values, pivot.index.values, pivot.values,
        shading="auto", cmap="RdYlGn", vmin=-10, vmax=vmax,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"Energy barrier $V_{max}$ ($k_BT$)")

    # Stability threshold contour
    cs = ax.contour(pivot.columns.values, pivot.index.values, pivot.values,
                     levels=[MARGINAL_THRESHOLD_KT, STABLE_THRESHOLD_KT],
                     colors=["orange", "black"], linewidths=1.2)
    ax.clabel(cs, inline=True, fontsize=8, fmt="%.0f kT")

    ax.set_yscale("log")
    ax.set_xlabel("pH")
    ax.set_ylabel("Ionic strength (mM, log scale)")
    ax.set_title("Predicted colloidal stability map (DLVO energy barrier)")

    # Overlay target food conditions we can't test directly
    for cond in TARGET_CONDITIONS:
        ax.scatter(cond["pH"], cond["ionic_strength_M"] * 1000.0,
                    marker="*", s=180, color="blue", edgecolor="white", zorder=5)
        ax.annotate(cond["label"], (cond["pH"], cond["ionic_strength_M"] * 1000.0),
                    textcoords="offset points", xytext=(6, 6), fontsize=8, color="blue")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="DLVO analysis for zein-caseinate NP (Sec. 3.1.3).")
    parser.add_argument("--anchor-csv", type=str, default=None,
                         help="Path to your measured DLS anchor matrix CSV.")
    parser.add_argument("--demo", action="store_true",
                         help="Use the bundled synthetic demo dataset instead of real data.")
    parser.add_argument("--hamaker", type=float, default=DEFAULT_HAMAKER_J,
                         help=f"Hamaker constant, J (default: {DEFAULT_HAMAKER_J:.2e}, "
                              "Su et al. 2024 zein NP value).")
    parser.add_argument("--outdir", type=str, default=OUTPUT_DIR,
                         help="Directory to write outputs to.")
    args = parser.parse_args()

    if not args.anchor_csv and not args.demo:
        parser.error("Provide --anchor-csv <path> or --demo.")

    csv_path = DEMO_CSV if args.demo else args.anchor_csv
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading anchor data from: {csv_path}")
    raw = load_anchor_data(csv_path)
    agg = aggregate_replicates(raw)
    print(f"  {len(raw)} raw rows -> {len(agg)} unique (pH, ionic strength) conditions")

    print("Fitting sigmoidal zeta(pH) model (pooled shape + per-I magnitude)...")
    fit = fit_zeta_ph_model(agg)
    zeta_low_pH, zeta_high_pH, pH_iep, slope = fit["pooled"]
    print(f"  pooled shape: zeta_low_pH={zeta_low_pH:.1f} mV, zeta_high_pH={zeta_high_pH:.1f} mV, "
          f"pH_iep={pH_iep:.2f}, slope={slope:.2f}")

    print("Evaluating DLVO on measured conditions...")
    measured_results = evaluate_measured_conditions(agg, args.hamaker)
    measured_out = os.path.join(args.outdir, "condition_summary.csv")
    measured_results.to_csv(measured_out, index=False)
    print(f"  -> {measured_out}")
    print(measured_results.to_string(index=False))

    print("Sweeping continuous pH x ionic-strength stability map...")
    radius_repr = agg["radius_nm_mean"].mean()
    sweep_df = sweep_stability_map(fit, args.hamaker, radius_repr)
    sweep_out = os.path.join(args.outdir, "stability_sweep_grid.csv")
    sweep_df.to_csv(sweep_out, index=False)
    print(f"  -> {sweep_out}")

    print("Evaluating specific target conditions (whole milk / yogurt / bread-release)...")
    target_rows = []
    for cond in TARGET_CONDITIONS:
        zeta = predict_zeta(cond["pH"], cond["ionic_strength_M"], fit)
        prof = dlvo.energy_profile(radius_repr, zeta, cond["ionic_strength_M"], args.hamaker)
        barrier = dlvo.find_energy_barrier(prof["D_nm"], prof["V_total_J"])
        stability = dlvo.classify_stability(barrier["V_max_kT"], STABLE_THRESHOLD_KT, MARGINAL_THRESHOLD_KT)
        target_rows.append({
            "condition": cond["label"], "pH": cond["pH"],
            "ionic_strength_mM": cond["ionic_strength_M"] * 1000.0,
            "zeta_mV_predicted": zeta, "V_max_kT": barrier["V_max_kT"],
            "stability_call": stability,
        })
    target_df = pd.DataFrame(target_rows)
    target_out = os.path.join(args.outdir, "target_condition_predictions.csv")
    target_df.to_csv(target_out, index=False)
    print(f"  -> {target_out}")
    print(target_df.to_string(index=False))

    print("Generating plots...")
    plot_energy_profiles(agg, args.hamaker, os.path.join(args.outdir, "energy_profiles.png"))
    plot_stability_map(sweep_df, os.path.join(args.outdir, "stability_map.png"))
    print(f"  -> {os.path.join(args.outdir, 'energy_profiles.png')}")
    print(f"  -> {os.path.join(args.outdir, 'stability_map.png')}")

    fit_params_out = os.path.join(args.outdir, "zeta_ph_fit_parameters.json")
    with open(fit_params_out, "w") as f:
        json.dump({
            "pooled_zeta_low_pH_mV": zeta_low_pH, "pooled_zeta_high_pH_mV": zeta_high_pH,
            "pooled_pH_iep": pH_iep, "pooled_slope": slope,
            "per_I_mM": {f"{k*1000:.0f}": v for k, v in fit["per_I"].items()},
            "hamaker_J": args.hamaker,
            "radius_nm_used_for_sweep": radius_repr,
        }, f, indent=2)
    print(f"  -> {fit_params_out}")

    print("\nDone.")


if __name__ == "__main__":
    main()
