"""
dlvo_model.py
=============

Core DLVO / extended-DLVO physics for protein-based food nanoparticles
(zein-caseinate NP, +/- WPI / casein exposure), written to implement the
workflow described in Section 3A ("DLVO / extended-DLVO modeling") of
"Dry-Lab Strategy for Section 3.1.3 - v3".

WHAT THIS MODULE DOES
----------------------
For a colloidal particle characterized by:
    - zeta potential  zeta   (mV, from electrophoretic light scattering / DLS instrument)
    - particle radius a      (nm, from DLS)
    - solution ionic strength I (mol/L, from your buffer/electrolyte recipe)
    - a literature Hamaker constant A (J, for protein-based colloids)

...this module computes the total pairwise interaction energy profile

    V(D) = V_EDL(D) + V_vdW(D)

as a function of surface-to-surface separation D, and reads off the
energy barrier V_max that determines kinetic (in)stability.

THIS IS A MECHANISTIC MODEL, NOT A REGRESSION. Each measured condition
(pH, I, zeta, a) is used to *parameterize* a physical equation; it is not
a training example for curve-fitting in the ML sense. See Section 4 of
the memo for the separate, genuinely data-driven XGBoost track.

GOVERNING EQUATIONS (all SI unless noted; see README.md for full derivation
notes, assumptions, and literature citations)
--------------------------------------------------------------------------
1. Debye (electric double layer) screening length, 1:1 (monovalent-equivalent)
   electrolyte:

        kappa = sqrt( 2 * N_A * e^2 * I_SI / (epsilon_0 * epsilon_r * k_B * T) )
        kappa_inverse = 1 / kappa                                    [Debye length]

   I_SI is ionic strength in mol/m^3 (= 1000 * I_molar).
   Reference: standard Debye-Huckel theory; used with this exact form in
   Su et al., J. Colloid Interface Sci. 2024, 656, 2321-2329 (zein NPs).

2. Zeta potential -> effective surface charge density (Grahame equation).
   Two options are provided:

     (a) Linear / Debye-Huckel approximation (valid for |zeta| below
         roughly 50-60 mV; this is the form used by Su et al. 2024 for
         zein NPs and is the one referenced in the memo):

             sigma = epsilon_0 * epsilon_r * kappa * zeta

     (b) Full nonlinear Grahame equation for a 1:1 electrolyte (more
         accurate at higher |zeta|, needed because some of our pH
         extremes push zein zeta beyond +-50 mV):

             sigma = sqrt(8 * epsilon_0 * epsilon_r * k_B * T * I_SI)
                     * sinh(e * zeta / (2 * k_B * T))

3. Electrostatic double-layer (EDL) repulsion between two equal spheres,
   Hogg-Healy-Fuerstenau (1966) constant-potential, weak-overlap
   approximation combined with the Derjaguin approximation. Valid when
   kappa*a >> 1, which holds here (a ~ 50-150 nm vs kappa^-1 ~ 0.8-3 nm
   across our 10-150 mM NaCl range):

        V_EDL(D) = 2 * pi * epsilon_0 * epsilon_r * a * zeta^2
                   * ln(1 + exp(-kappa * D))

   Reference: Hogg, Healy & Fuerstenau, Trans. Faraday Soc. 1966, 62,
   1638; reproduced in Israelachvili, "Intermolecular and Surface
   Forces," 3rd ed., Eq. 14.46.

4. Van der Waals attraction between two equal spheres, exact
   (non-retarded) Hamaker expression -- valid at all separations D, not
   just D << a:

        V_vdW(D) = -(A/6) * [ 2a^2/(D*(D+4a))
                               + 2a^2/(D+2a)^2
                               + ln( D*(D+4a) / (D+2a)^2 ) ]

   Reference: Israelachvili, "Intermolecular and Surface Forces," 3rd
   ed., Eq. 13.14. A is the Hamaker constant (J); literature values for
   protein-based colloids are roughly 1e-20 to 5e-21 J (memo Section 3A);
   Su et al. 2024 fit A = 1.0e-20 +/- 0.4e-20 J specifically for zein
   NPs by regression against measured critical coagulation concentrations.

5. Total interaction energy and stability criterion:

        V(D) = V_EDL(D) + V_vdW(D)
        V_max = height of the primary energy barrier (local maximum of V)

        V_max > ~15 kT  -> kinetically stable (slow aggregation)
        V_max ~10-15 kT -> marginal / borderline
        V_max < ~10 kT  -> fast aggregation expected

   (kT = k_B * T; "kT units" below always means V_max / (k_B * T).)

6. Because zeta(pH) follows the protein's known sigmoidal isoelectric-
   point curve, a handful of measured zeta values (per ionic strength)
   let us interpolate/extrapolate the full pH response and then sweep
   DLVO across a continuous pH x ionic-strength grid -- including
   conditions we cannot test directly right now (e.g. skim milk).

        zeta(pH) = zeta_max + (zeta_min - zeta_max) / (1 + exp((pH - pH_iep) / slope))

LIMITATIONS / ASSUMPTIONS (read before trusting numbers -- see README.md
for the full discussion):
  - Spheres of uniform radius `a` (DLS z-average radius) and uniform
    surface potential are assumed; real zein-caseinate NPs are not
    perfectly spherical or homogeneous.
  - Constant-potential boundary condition for V_EDL (as opposed to
    constant-charge); the two bracket the true (charge-regulated)
    behavior and typically differ by a factor of order unity near
    contact, negligibly at the barrier itself for kappa*a >> 1.
  - Non-retarded van der Waals (no zero-frequency retardation
    correction); adequate at the D ~ 1-10 nm separations relevant to
    the primary barrier for these particle sizes.
  - Single, system-wide Hamaker constant A; does not capture protein
    layer heterogeneity or bound-water ("hydration force") effects that
    the memo's Molina-Bolivar et al. 2001 reference documents at high
    ionic strength. If your measured stability disagrees with the model
    at high I, that paper's short-range hydration repulsion is the
    most likely missing term -- flag it, don't silently force-fit A.
  - epsilon_r (water relative permittivity) is treated as a constant at
    its ~298 K value unless you override it; adjust if you run the
    63 degC pasteurization-relevant condition.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

# ---------------------------------------------------------------------------
# Physical constants (SI units, CODATA values)
# ---------------------------------------------------------------------------
K_B = 1.380649e-23          # Boltzmann constant, J/K
E_CHARGE = 1.602176634e-19  # elementary charge, C
N_A = 6.02214076e23         # Avogadro's number, 1/mol
EPSILON_0 = 8.8541878128e-12  # vacuum permittivity, F/m
EPSILON_R_WATER = 78.5      # relative permittivity of water at ~25 degC (dimensionless)

# Literature Hamaker constant range for protein-based colloids (J)
# (memo Section 3A; see also Su et al. 2024 fitted value of 1.0e-20 +/- 0.4e-20 J
# for zein NPs specifically -- a good default if you have no better estimate).
HAMAKER_RANGE_PROTEIN_J = (5e-21, 1e-20)
HAMAKER_ZEIN_DEFAULT_J = 1.0e-20


# ---------------------------------------------------------------------------
# 1. Debye length
# ---------------------------------------------------------------------------
def debye_length_nm(ionic_strength_M: float, temperature_K: float = 298.15,
                     epsilon_r: float = EPSILON_R_WATER) -> float:
    """
    Debye (electric double layer) screening length for a 1:1-equivalent
    electrolyte.

    Parameters
    ----------
    ionic_strength_M : float
        Ionic strength, mol/L (e.g. 0.05 for 50 mM NaCl).
    temperature_K : float
        Absolute temperature, K.
    epsilon_r : float
        Relative permittivity of the solvent (water ~78.5 at 25 degC).

    Returns
    -------
    float
        Debye length (1/kappa), in nm.
    """
    if ionic_strength_M <= 0:
        raise ValueError("ionic_strength_M must be > 0")
    I_SI = ionic_strength_M * 1000.0  # mol/L -> mol/m^3
    kappa = np.sqrt(
        (2.0 * N_A * E_CHARGE**2 * I_SI)
        / (EPSILON_0 * epsilon_r * K_B * temperature_K)
    )  # 1/m
    return 1.0 / kappa * 1e9  # m -> nm


def debye_kappa_per_nm(ionic_strength_M: float, temperature_K: float = 298.15,
                        epsilon_r: float = EPSILON_R_WATER) -> float:
    """Inverse Debye length kappa, in 1/nm (convenience wrapper)."""
    return 1.0 / debye_length_nm(ionic_strength_M, temperature_K, epsilon_r)


# ---------------------------------------------------------------------------
# 2. Zeta potential -> surface charge density (Grahame equation)
# ---------------------------------------------------------------------------
def surface_charge_density(zeta_mV: float, ionic_strength_M: float,
                            temperature_K: float = 298.15,
                            epsilon_r: float = EPSILON_R_WATER,
                            method: str = "linear") -> float:
    """
    Convert zeta potential to an effective surface charge density sigma
    (C/m^2) via the Grahame equation.

    Parameters
    ----------
    zeta_mV : float
        Zeta potential, mV (can be negative).
    ionic_strength_M : float
        Ionic strength, mol/L.
    method : {"linear", "grahame"}
        "linear"  -- Debye-Huckel linear approximation, sigma = eps0*epsR*kappa*zeta.
                     Valid for |zeta| below ~50-60 mV. Matches Su et al. 2024 (zein NPs).
        "grahame" -- full nonlinear Grahame equation for a 1:1 electrolyte.
                     Use this if |zeta| exceeds ~50-60 mV (e.g. zein at pH 4 or 9).

    Returns
    -------
    float
        Surface charge density, C/m^2.
    """
    zeta_V = zeta_mV * 1e-3
    if method == "linear":
        kappa = debye_kappa_per_nm(ionic_strength_M, temperature_K, epsilon_r) * 1e9  # 1/m
        return EPSILON_0 * epsilon_r * kappa * zeta_V
    elif method == "grahame":
        I_SI = ionic_strength_M * 1000.0
        prefactor = np.sqrt(8.0 * EPSILON_0 * epsilon_r * K_B * temperature_K * I_SI)
        return prefactor * np.sinh(E_CHARGE * zeta_V / (2.0 * K_B * temperature_K))
    else:
        raise ValueError('method must be "linear" or "grahame"')


# ---------------------------------------------------------------------------
# 3. Electrostatic double-layer (EDL) repulsion energy
# ---------------------------------------------------------------------------
def edl_energy_J(D_nm, radius_nm: float, zeta_mV: float, ionic_strength_M: float,
                  temperature_K: float = 298.15, epsilon_r: float = EPSILON_R_WATER):
    """
    Electrostatic double-layer repulsion energy between two equal spheres
    (Hogg-Healy-Fuerstenau constant-potential, weak-overlap / Derjaguin
    approximation). See module docstring, eq. 3.

    Parameters
    ----------
    D_nm : float or np.ndarray
        Surface-to-surface separation(s), nm. Must be > 0.
    radius_nm : float
        Particle radius, nm (DLS z-average).
    zeta_mV : float
        Zeta potential, mV.
    ionic_strength_M : float
        Ionic strength, mol/L.

    Returns
    -------
    float or np.ndarray
        V_EDL(D), Joules (always >= 0, repulsive).
    """
    D_nm = np.asarray(D_nm, dtype=float)
    a_m = radius_nm * 1e-9
    D_m = D_nm * 1e-9
    zeta_V = zeta_mV * 1e-3
    kappa_per_m = debye_kappa_per_nm(ionic_strength_M, temperature_K, epsilon_r) * 1e9

    return (2.0 * np.pi * EPSILON_0 * epsilon_r * a_m * zeta_V**2
            * np.log(1.0 + np.exp(-kappa_per_m * D_m)))


# ---------------------------------------------------------------------------
# 4. Van der Waals attraction energy
# ---------------------------------------------------------------------------
def vdw_energy_J(D_nm, radius_nm: float, hamaker_J: float):
    """
    Non-retarded van der Waals attraction energy between two equal
    spheres, exact Hamaker sphere-sphere expression (valid for all D, not
    just D << a). See module docstring, eq. 4.

    Parameters
    ----------
    D_nm : float or np.ndarray
        Surface-to-surface separation(s), nm. Must be > 0.
    radius_nm : float
        Particle radius, nm.
    hamaker_J : float
        Hamaker constant, J (positive number; the function returns a
        negative, i.e. attractive, energy).

    Returns
    -------
    float or np.ndarray
        V_vdW(D), Joules (always <= 0, attractive).
    """
    D_nm = np.asarray(D_nm, dtype=float)
    a = radius_nm  # keep everything in nm here; it's a dimensionless ratio anyway
    D = D_nm

    term1 = 2.0 * a**2 / (D * (D + 4.0 * a))
    term2 = 2.0 * a**2 / (D + 2.0 * a)**2
    term3 = np.log((D * (D + 4.0 * a)) / (D + 2.0 * a)**2)

    return -(hamaker_J / 6.0) * (term1 + term2 + term3)


# ---------------------------------------------------------------------------
# 5. Total interaction energy + energy-barrier extraction
# ---------------------------------------------------------------------------
def total_energy_J(D_nm, radius_nm: float, zeta_mV: float, ionic_strength_M: float,
                    hamaker_J: float, temperature_K: float = 298.15,
                    epsilon_r: float = EPSILON_R_WATER):
    """V_total(D) = V_EDL(D) + V_vdW(D), in Joules."""
    v_edl = edl_energy_J(D_nm, radius_nm, zeta_mV, ionic_strength_M, temperature_K, epsilon_r)
    v_vdw = vdw_energy_J(D_nm, radius_nm, hamaker_J)
    return v_edl + v_vdw


def energy_profile(radius_nm: float, zeta_mV: float, ionic_strength_M: float,
                    hamaker_J: float, D_min_nm: float = 0.3, D_max_nm=None,
                    n_points: int = 3000, temperature_K: float = 298.15,
                    epsilon_r: float = EPSILON_R_WATER) -> dict:
    """
    Compute the full V(D) interaction profile for one (pH-implicit)
    condition, i.e. one (radius, zeta, ionic strength, Hamaker constant)
    combination.

    D_min_nm=0.3 nm is used as a hard-wall / Born-repulsion cutoff: the
    vdW term formally diverges as D -> 0, which is unphysical (real
    surfaces have finite roughness/hydration layers); 0.3 nm is a
    standard, conservative choice in the DLVO literature for this cutoff.

    D_max_nm defaults to max(20 * Debye length, 50 nm) so the barrier
    search window comfortably contains the barrier for typical
    food-relevant ionic strengths (10-150 mM NaCl).

    Returns
    -------
    dict with keys:
        D_nm, V_EDL_J, V_vdW_J, V_total_J, V_total_kT (all np.ndarray),
        kappa_inv_nm (float), temperature_K (float)
    """
    kappa_inv_nm = debye_length_nm(ionic_strength_M, temperature_K, epsilon_r)
    if D_max_nm is None:
        D_max_nm = max(20.0 * kappa_inv_nm, 50.0)

    D_nm = np.linspace(D_min_nm, D_max_nm, n_points)
    v_edl = edl_energy_J(D_nm, radius_nm, zeta_mV, ionic_strength_M, temperature_K, epsilon_r)
    v_vdw = vdw_energy_J(D_nm, radius_nm, hamaker_J)
    v_total = v_edl + v_vdw
    kT = K_B * temperature_K

    return {
        "D_nm": D_nm,
        "V_EDL_J": v_edl,
        "V_vdW_J": v_vdw,
        "V_total_J": v_total,
        "V_total_kT": v_total / kT,
        "kappa_inv_nm": kappa_inv_nm,
        "temperature_K": temperature_K,
    }


def find_energy_barrier(D_nm: np.ndarray, V_total_J: np.ndarray,
                         temperature_K: float = 298.15) -> dict:
    """
    Locate the primary energy barrier (local maximum) in a V(D) profile.

    Returns
    -------
    dict with keys:
        D_max_nm       : separation at the barrier, nm
        V_max_J        : barrier height, J
        V_max_kT       : barrier height, kT units
        barrier_present: bool -- True if a genuine local maximum was found
                         strictly inside the scanned window (i.e. not just
                         the profile still rising at the D_max boundary,
                         and not just monotonically falling from D_min).
    """
    kT = K_B * temperature_K
    idx_max = int(np.argmax(V_total_J))
    n = len(V_total_J)

    interior_max = 0 < idx_max < (n - 1)
    if not interior_max:
        # Either monotonically decreasing from D_min (no barrier at all),
        # or still rising at D_max (window too narrow -- caller should widen it).
        return {
            "D_max_nm": float(D_nm[idx_max]),
            "V_max_J": float(V_total_J[idx_max]),
            "V_max_kT": float(V_total_J[idx_max] / kT),
            "barrier_present": False,
        }

    return {
        "D_max_nm": float(D_nm[idx_max]),
        "V_max_J": float(V_total_J[idx_max]),
        "V_max_kT": float(V_total_J[idx_max] / kT),
        "barrier_present": True,
    }


def classify_stability(V_max_kT: float, stable_threshold_kT: float = 15.0,
                        marginal_threshold_kT: float = 10.0) -> str:
    """
    Apply the memo's stability heuristic:
        V_max > ~15 kT   -> "stable"     (kinetically stable, slow aggregation)
        10-15 kT         -> "marginal"   (borderline; flag for wet-lab check)
        V_max < ~10 kT   -> "unstable"   (fast aggregation expected)
    """
    if V_max_kT >= stable_threshold_kT:
        return "stable"
    elif V_max_kT >= marginal_threshold_kT:
        return "marginal"
    else:
        return "unstable"


# ---------------------------------------------------------------------------
# 6. Sigmoidal zeta(pH) model -- lets a handful of measured points be
#    interpolated/extrapolated across a continuous pH grid.
# ---------------------------------------------------------------------------
def sigmoid_zeta_ph(pH, zeta_low_pH, zeta_high_pH, pH_iep, slope):
    """
    Boltzmann-sigmoid model of zeta potential vs. pH around an
    isoelectric point (IEP).

    zeta(pH) = zeta_high_pH + (zeta_low_pH - zeta_high_pH) / (1 + exp((pH - pH_iep) / slope))

    As pH -> -infinity (i.e. deep below the IEP), zeta -> zeta_low_pH
    (typically POSITIVE for proteins like zein, since amino groups are
    protonated at low pH).
    As pH -> +infinity (deep above the IEP), zeta -> zeta_high_pH
    (typically NEGATIVE, carboxyl groups deprotonated at high pH).
    `slope` controls the steepness of the transition (smaller = steeper).

    NOTE on naming: parameters are named by WHICH SIDE of the IEP they
    describe (low_pH / high_pH), not by numeric sign or magnitude --
    zeta_low_pH is usually the larger, positive plateau and zeta_high_pH
    the more negative one, but the function does not assume this; it will
    fit whatever sigmoid shape the data show.
    """
    pH = np.asarray(pH, dtype=float)
    return zeta_high_pH + (zeta_low_pH - zeta_high_pH) / (1.0 + np.exp((pH - pH_iep) / slope))


def fit_zeta_vs_pH(pH_data, zeta_data, p0=None):
    """
    Fit `sigmoid_zeta_ph` to measured (pH, zeta) anchor points.

    Parameters
    ----------
    pH_data, zeta_data : array-like
        Measured anchor points (e.g. the pH-grid rows of your DLS matrix,
        typically averaged across replicates, for ONE fixed ionic
        strength -- fit a separate curve per ionic-strength group).
    p0 : tuple, optional
        Initial guess (zeta_min, zeta_max, pH_iep, slope). If omitted, a
        reasonable guess is derived from the data.

    Returns
    -------
    (popt, pcov) : fitted parameters and covariance matrix from
    scipy.optimize.curve_fit. popt = (zeta_low_pH, zeta_high_pH, pH_iep, slope).
    """
    pH_data = np.asarray(pH_data, dtype=float)
    zeta_data = np.asarray(zeta_data, dtype=float)

    if p0 is None:
        # Order anchor points by pH to guess which plateau is "low_pH" vs
        # "high_pH" from the data itself, rather than assuming sign.
        order = np.argsort(pH_data)
        zeta_low_pH_guess = float(zeta_data[order][0])
        zeta_high_pH_guess = float(zeta_data[order][-1])
        pH_iep_guess = float(pH_data[np.argmin(np.abs(zeta_data))]) if len(pH_data) else float(np.mean(pH_data))
        slope_guess = 0.5
        p0 = (zeta_low_pH_guess, zeta_high_pH_guess, pH_iep_guess, slope_guess)

    popt, pcov = curve_fit(sigmoid_zeta_ph, pH_data, zeta_data, p0=p0, maxfev=10000)
    return popt, pcov


def fit_zeta_vs_pH_fixed_shape(pH_data, zeta_data, pH_iep: float, slope: float, p0=None):
    """
    Fit ONLY (zeta_min, zeta_max) with pH_iep and slope held fixed.

    Why this exists: the designed anchor grid (Section 2 of the memo) has
    only 3 pH points per ionic-strength group. The full 4-parameter
    sigmoid (zeta_min, zeta_max, pH_iep, slope) is under-determined with
    3 points. Standard practice here: fit the full 4-parameter sigmoid
    once on the POOLED data (all ionic strengths together, via
    fit_zeta_vs_pH) to pin down the shape (pH_iep, slope) -- which mostly
    reflects the protein's intrinsic protonation chemistry and should not
    depend strongly on ionic strength -- then use this function to refit
    only the 2 remaining, genuinely I-dependent parameters (zeta_min,
    zeta_max, which capture double-layer screening magnitude) separately
    for each ionic-strength group. 2 free parameters vs. 3 data points is
    well-posed.

    Parameters
    ----------
    pH_data, zeta_data : array-like
        Measured (pH, zeta) anchor points for ONE ionic-strength group.
    pH_iep, slope : float
        Fixed shape parameters, typically taken from a pooled
        fit_zeta_vs_pH() call across all groups (or from literature,
        e.g. zein IEP ~5.8; Su et al. 2024).
    p0 : tuple, optional
        Initial guess (zeta_min, zeta_max).

    Returns
    -------
    (popt, pcov) : popt = (zeta_low_pH, zeta_high_pH)
    """
    pH_data = np.asarray(pH_data, dtype=float)
    zeta_data = np.asarray(zeta_data, dtype=float)

    def _model(pH, zeta_low_pH, zeta_high_pH):
        return sigmoid_zeta_ph(pH, zeta_low_pH, zeta_high_pH, pH_iep, slope)

    if p0 is None:
        order = np.argsort(pH_data)
        p0 = (float(zeta_data[order][0]), float(zeta_data[order][-1]))

    popt, pcov = curve_fit(_model, pH_data, zeta_data, p0=p0, maxfev=10000)
    return popt, pcov


if __name__ == "__main__":
    # Minimal self-test / sanity check (run: python dlvo_model.py)
    # Uses parameters representative of a zein-caseinate NP at pH 4,
    # 10 mM NaCl (near the low end of our designed grid).
    a_nm = 75.0
    zeta = -35.0
    I = 0.010
    A = HAMAKER_ZEIN_DEFAULT_J

    prof = energy_profile(a_nm, zeta, I, A)
    barrier = find_energy_barrier(prof["D_nm"], prof["V_total_J"])
    stability = classify_stability(barrier["V_max_kT"])

    print(f"Debye length: {prof['kappa_inv_nm']:.2f} nm")
    print(f"Barrier: D = {barrier['D_max_nm']:.2f} nm, "
          f"V_max = {barrier['V_max_kT']:.1f} kT, present={barrier['barrier_present']}")
    print(f"Stability call: {stability}")
