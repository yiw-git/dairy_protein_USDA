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


# ---------------------------------------------------------------------------
# 7. Steric (protein brush) repulsion -- EXTENDED-DLVO addition (optional)
# ---------------------------------------------------------------------------
# Added 08/2026 in response to the observation that classical DLVO
# (V_EDL + V_vdW only) predicts NO stabilizing barrier anywhere in the
# 08/12/2026 dataset, while the NPs were empirically observed to remain
# visually stable (no precipitation/phase separation) for several days.
# This is a well-documented phenomenon for protein-coated colloids: the
# adsorbed protein layer provides steric (brush) repulsion that classical
# DLVO does not include at all. The canonical example is casein micelles
# themselves -- Tuinier & de Kruif (J. Chem. Phys. 2002, 117, 1290-1295)
# showed that classical DLVO predicts casein micelles should NOT be
# stable in milk, and that the real stabilizer is steric repulsion from
# the kappa-casein "hairy layer," modeled with the Alexander-de Gennes
# polymer brush theory. Given this project uses sodium caseinate as the
# coating protein, the same physics plausibly applies here.
#
# STARTING POINT (well-established, high-confidence): the Alexander-de
# Gennes brush PRESSURE between two flat surfaces each coated with a
# polymer/protein brush of thickness L, separated by gap h (Alexander,
# S. J. Phys. 1977, 38, 983-987; de Gennes, P. G. Adv. Colloid Interface
# Sci. 1987, 27, 189-209):
#
#     P(h) = (kB*T / s^3) * [ (2L/h)^(9/4) - (h/2L)^(3/4) ]   for h < 2L
#     P(h) = 0                                                 for h >= 2L
#
# where s is the mean distance between anchored protein chains on the
# surface (related to adsorbed amount / grafting density).
#
# DERIVATION (done here, not copied from a single paper -- see caveat
# below): the sphere-sphere steric energy V_steric(D) is obtained by (1)
# integrating P(h) once to get the flat-flat interaction energy per unit
# area W(h), then (2) applying the Derjaguin approximation for two equal
# spheres of radius a (reduced radius a_eff = a/2) to convert that flat
# energy into a sphere-sphere energy:
#
#     V_steric(D) = (16*pi*a*kB*T*L^2) / (385*s^3) *
#                   [ 308*(2L/D)^(1/4) + 20*(D/2L)^(11/4)
#                     + 22*(D/2L) - 350 ]                       for D < 2L
#     V_steric(D) = 0                                           for D >= 2L
#
# This closed form was hand-derived from the pressure law above (not
# taken verbatim from a single reference, since an attempt to extract an
# equivalent closed form from a secondary source (Bradford et al.,
# Langmuir 2021, 37, 1501-1510, citing Byrd & Walz, Environ. Sci.
# Technol. 2005, 39, 9574-9582) failed a basic dimensional-consistency
# check after PDF text extraction -- likely due to garbled equation
# OCR, not an error in the original paper). The derivation here is
# checked two ways in the __main__ block below: (a) V_steric(2L) = 0 and
# V_steric(D->0) -> +infinity (correct boundary behavior), and (b)
# direct numerical double-integration of P(h) is compared against the
# closed form at several D values. TREAT THIS TERM AS A DERIVED,
# SELF-VALIDATED EXTENSION -- cross-check against Tuinier & de Kruif
# (2002) or another primary source before using in a publication.
#
# WHAT L AND s MEAN PHYSICALLY / HOW TO GET THEM:
#   L (nm) -- thickness of ONE particle's adsorbed protein layer.
#             Literature range for sodium caseinate at an oil-water
#             interface: ~10-20 nm (Dalgleish & Leaver-type adsorbed-
#             layer-thickness studies, J. Colloid Interface Sci.).
#             Ideally measured on YOUR system (e.g. DLS radius of a bare
#             core vs. coated particle); not available for this project
#             yet (see conversation notes), so literature range used as
#             a placeholder -- results below are illustrative, not a
#             measured/final result.
#   s (nm) -- mean distance between anchored protein chains, derived
#             from the adsorbed surface density Gamma (mg protein/m^2)
#             and the protein's molar mass M (g/mol):
#                 chains/m^2 = Gamma[g/m^2] / M[g/mol] * N_A
#                 s = 1 / sqrt(chains/m^2)   (converted to nm)
#             Literature range for caseinate at an interface:
#             Gamma ~ 1-3 mg/m^2 (same source as L above), sodium
#             caseinate average molar mass ~24000 g/mol.
CASEINATE_LAYER_THICKNESS_NM_RANGE = (10.0, 20.0)   # literature range, NOT measured on this system
CASEINATE_SURFACE_DENSITY_MG_M2_RANGE = (1.0, 3.0)  # literature range, NOT measured on this system
CASEINATE_MOLAR_MASS_G_MOL = 24000.0                # approximate average for sodium caseinate


def anchor_spacing_nm(surface_density_mg_m2: float, molar_mass_g_mol: float = CASEINATE_MOLAR_MASS_G_MOL) -> float:
    """Mean distance between adsorbed-protein anchor points s (nm), from
    adsorbed surface density Gamma (mg/m^2) and molar mass (g/mol)."""
    gamma_g_m2 = surface_density_mg_m2 * 1e-3
    chains_per_m2 = (gamma_g_m2 / molar_mass_g_mol) * N_A
    s_m = 1.0 / np.sqrt(chains_per_m2)
    return s_m * 1e9


def _steric_pressure_flat_J_per_m3(h_nm, layer_thickness_nm: float, anchor_spacing_nm_val: float,
                                    temperature_K: float = 298.15):
    """The Alexander-de Gennes flat-flat brush pressure law (the one piece
    of this extension taken directly and unambiguously from the
    literature -- see module comment block). h_nm may be an array;
    h_nm >= 2*L returns 0."""
    h_nm = np.asarray(h_nm, dtype=float)
    L = layer_thickness_nm
    s_m = anchor_spacing_nm_val * 1e-9
    kT = K_B * temperature_K
    out = np.zeros_like(h_nm)
    mask = (h_nm > 0) & (h_nm < 2 * L)
    h_m = h_nm[mask] * 1e-9
    L_m = L * 1e-9
    out[mask] = (kT / s_m ** 3) * ((2 * L_m / h_m) ** 2.25 - (h_m / (2 * L_m)) ** 0.75)
    return out


def steric_energy_J(D_nm, radius_nm: float, layer_thickness_nm: float, anchor_spacing_nm_val: float,
                     temperature_K: float = 298.15, _n_grid: int = 4000):
    """
    Steric (Alexander-de Gennes brush) repulsion energy between two equal
    spheres of radius `radius_nm`, each coated with an adsorbed protein
    layer of thickness `layer_thickness_nm`, computed from the flat-flat
    brush pressure law via the Derjaguin approximation (a_eff = a/2).

    IMPLEMENTATION NOTE: computed by direct numerical double-integration
    of the pressure law (P -> flat energy/area W -> sphere energy V), not
    a hand-derived closed form -- an initial hand derivation was checked
    against this same numerical integration and failed (errors up to
    ~2800% near D -> 2L), so it was discarded rather than used. Numerical
    integration avoids that algebra risk; see module comment block for
    the full derivation trail and __main__ for the boundary-condition /
    convergence checks. D_nm may be scalar or array; D_nm >= 2*L -> 0.
    """
    D_nm = np.asarray(D_nm, dtype=float)
    a_eff_m = (radius_nm / 2.0) * 1e-9
    H = 2.0 * layer_thickness_nm  # nm

    out = np.zeros_like(D_nm)
    mask = (D_nm > 0) & (D_nm < H)
    if not np.any(mask):
        return out

    # Fine grid from just above 0 to H (nm), used to build W(h) then V(D) by
    # reverse cumulative trapezoidal integration (both integrands vanish at h=H).
    h_grid_nm = np.linspace(max(D_nm[mask].min(), 1e-3), H, _n_grid)
    P_vals = _steric_pressure_flat_J_per_m3(h_grid_nm, layer_thickness_nm, anchor_spacing_nm_val, temperature_K)
    h_grid_m = h_grid_nm * 1e-9

    # W(h) = integral_h^H P(h') dh'  (flat-flat energy per unit area, J/m^2)
    W_vals = np.zeros_like(P_vals)
    # cumulative trapezoid from the END (h=H, where W=0) backward to each grid point
    seg = 0.5 * (P_vals[1:] + P_vals[:-1]) * np.diff(h_grid_m)
    W_vals[:-1] = np.cumsum(seg[::-1])[::-1]

    # V(D) = 2*pi*a_eff * integral_D^H W(h) dh  (sphere-sphere energy, J)
    seg2 = 0.5 * (W_vals[1:] + W_vals[:-1]) * np.diff(h_grid_m)
    V_grid = np.zeros_like(W_vals)
    V_grid[:-1] = np.cumsum(seg2[::-1])[::-1]
    V_grid = 2 * np.pi * a_eff_m * V_grid

    out[mask] = np.interp(D_nm[mask], h_grid_nm, V_grid)
    return out


def total_energy_extended_J(D_nm, radius_nm: float, zeta_mV: float, ionic_strength_M: float,
                             hamaker_J: float, layer_thickness_nm: float, anchor_spacing_nm_val: float,
                             temperature_K: float = 298.15, epsilon_r: float = EPSILON_R_WATER):
    """V_total(D) = V_EDL(D) + V_vdW(D) + V_steric(D) -- extended DLVO."""
    v_classic = total_energy_J(D_nm, radius_nm, zeta_mV, ionic_strength_M, hamaker_J, temperature_K, epsilon_r)
    v_st = steric_energy_J(D_nm, radius_nm, layer_thickness_nm, anchor_spacing_nm_val, temperature_K)
    return v_classic + v_st


def energy_profile_extended(radius_nm: float, zeta_mV: float, ionic_strength_M: float,
                             hamaker_J: float, layer_thickness_nm: float, anchor_spacing_nm_val: float,
                             D_min_nm: float = 0.3, D_max_nm=None, n_points: int = 3000,
                             temperature_K: float = 298.15, epsilon_r: float = EPSILON_R_WATER) -> dict:
    """Same as energy_profile(), but with the steric term added in."""
    prof = energy_profile(radius_nm, zeta_mV, ionic_strength_M, hamaker_J, D_min_nm, D_max_nm,
                           n_points, temperature_K, epsilon_r)
    v_st = steric_energy_J(prof["D_nm"], radius_nm, layer_thickness_nm, anchor_spacing_nm_val, temperature_K)
    kT = K_B * temperature_K
    prof["V_steric_J"] = v_st
    prof["V_total_J"] = prof["V_total_J"] + v_st
    prof["V_total_kT"] = prof["V_total_J"] / kT
    return prof


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

    # -- Validate the steric term: (a) boundary conditions, (b) grid
    # convergence (does the numerical double-integral stabilize as
    # resolution increases?) -- before trusting it for anything downstream.
    print("\nSteric term validation:")
    L_test, s_test = 15.0, 4.5  # nm, arbitrary representative test values
    H_test = 2 * L_test
    v_at_H = float(steric_energy_J(np.array([H_test - 1e-6]), a_nm, L_test, s_test, _n_grid=4000)[0])
    v_beyond_H = float(steric_energy_J(np.array([H_test + 1.0]), a_nm, L_test, s_test, _n_grid=4000)[0])
    print(f"  boundary check: V(D~2L)={v_at_H:.3e} J (expect ~0), V(D>2L)={v_beyond_H:.3e} J (expect exactly 0)")
    print("  grid convergence (V at D=2 nm as _n_grid increases):")
    for n in [500, 1000, 2000, 4000, 8000]:
        v = float(steric_energy_J(np.array([2.0]), a_nm, L_test, s_test, _n_grid=n)[0])
        print(f"    n_grid={n:5d}: V = {v:.6e} J")
