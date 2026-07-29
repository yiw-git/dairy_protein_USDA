# DLVO / extended-DLVO modeling -- Section 3.1.3

Implements the "A. DLVO / extended-DLVO modeling" workflow from Section 3A
of *Dry-Lab Strategy for Section 3.1.3 - v3*, for the zein-caseinate NP
(bare and after WPI / casein exposure).

This is a **mechanistic model, not a regression**: each measured DLS
condition (pH, ionic strength, zeta potential, particle radius) is used
to parameterize a physical equation, not as a training row for
curve-fitting in the ML sense. The genuinely data-driven track (XGBoost
zeta-potential model) is a separate, parallel effort described in
Section 4 of the memo -- this script is the complementary, physics-based
leg of the "triangulation" story.

## What it computes

1. **Debye length** from ionic strength.
2. **Surface charge density** from zeta potential (Grahame equation).
3. **Electrostatic double-layer (EDL) repulsion** V_EDL(D) between two
   spheres.
4. **Van der Waals attraction** V_vdW(D) between two spheres (Hamaker
   constant).
5. **Total interaction energy** V(D) = V_EDL(D) + V_vdW(D), and the
   **energy barrier** V_max that governs kinetic stability.
6. A **sigmoidal zeta(pH) fit** so that a handful of measured anchor
   points can be interpolated/extrapolated across a continuous pH x
   ionic-strength grid -- including food conditions you can't test
   directly right now (skim milk, whole milk, yogurt, bread-release).

## Files

| File | Purpose |
|---|---|
| `dlvo_model.py` | Core physics: every equation, fully documented, unit-tested via its own `__main__` block. Import this in your own notebooks/scripts too. |
| `run_dlvo_analysis.py` | End-to-end pipeline: load anchor CSV -> fit zeta(pH) -> evaluate measured conditions -> sweep continuous stability map -> save CSVs + plots. |
| `anchor_data_template.csv` | **Blank template matching the Section 2 designed grid** (3 pH x 3 ionic strengths x bare/protein-exposed x 2 replicates = 36 rows). Fill in `zeta_mV` and `radius_nm` from your DLS run. |
| `demo_data_zein_caseinate.csv` | **Synthetic** data (NOT real measurements) for testing the pipeline before your wet-lab window. Loosely grounded in the zein pH-zeta behavior reported in Su et al. 2024 (IEP ~5.8; see References). |
| `requirements.txt` | Python dependencies. |
| `outputs/` | Created automatically; holds all generated CSVs and plots. |

## Quick start

```bash
pip install -r requirements.txt

# 1. Sanity-check the physics module on its own:
python dlvo_model.py

# 2. Run the full pipeline on synthetic demo data (no real DLS data needed):
python run_dlvo_analysis.py --demo

# 3. Once you have real DLS measurements, fill in anchor_data_template.csv
#    (duplicate it first, e.g. anchor_data_measured.csv) and run:
python run_dlvo_analysis.py --anchor-csv anchor_data_measured.csv
```

Outputs land in `outputs/`:
- `condition_summary.csv` -- DLVO barrier + stability call for every measured condition.
- `stability_sweep_grid.csv` -- full continuous pH x ionic-strength sweep (raw numbers behind the heatmap).
- `target_condition_predictions.csv` -- predictions for whole milk / yogurt / bread-release (the Section 4 conditions).
- `energy_profiles.png` -- V(D) curves for every measured condition, with the 10 kT / 15 kT stability thresholds marked.
- `stability_map.png` -- pH x ionic-strength heatmap of the energy barrier, with the target food conditions overlaid as stars.
- `zeta_ph_fit_parameters.json` -- the fitted sigmoid parameters, for reuse/inspection.

## Governing equations (summary -- see `dlvo_model.py` docstrings for full derivations and exact citations)

**Debye length** (1:1-equivalent electrolyte):

```
kappa = sqrt( 2 * N_A * e^2 * I / (epsilon_0 * epsilon_r * k_B * T) )
```

**Grahame equation** (zeta -> surface charge density), linear (Debye-Huckel) form used by default:

```
sigma = epsilon_0 * epsilon_r * kappa * zeta
```

A full nonlinear form is also implemented (`method="grahame"` in
`surface_charge_density()`) for the higher |zeta| conditions (zein at
pH 4 or pH 7 routinely exceeds the ~50-60 mV validity range of the linear
approximation).

**EDL repulsion**, Hogg-Healy-Fuerstenau constant-potential / Derjaguin
approximation (valid because kappa*a >> 1 here: a ~ 50-150 nm vs.
kappa^-1 ~ 0.8-3 nm across 10-150 mM NaCl):

```
V_EDL(D) = 2*pi*epsilon_0*epsilon_r*a*zeta^2 * ln(1 + exp(-kappa*D))
```

**Van der Waals attraction**, exact (non-retarded) Hamaker sphere-sphere
expression (valid at all D, not just D << a):

```
V_vdW(D) = -(A/6) * [ 2a^2/(D(D+4a)) + 2a^2/(D+2a)^2 + ln(D(D+4a)/(D+2a)^2) ]
```

**Stability call** (memo heuristic):

```
V_max > ~15 kT   ->  stable    (kinetically stable, slow aggregation)
V_max 10-15 kT   ->  marginal  (borderline; flag for wet-lab check)
V_max < ~10 kT   ->  unstable  (fast aggregation expected)
```

## Hamaker constant

Default: **1.0e-20 J**, the value Su et al. (2024) fit specifically for
zein NPs by regressing critical coagulation concentrations against
surface charge density. The broader literature range for protein-based
colloids is **5e-21 to 1e-20 J** (memo Section 3A) -- both are available
as constants in `dlvo_model.py` (`HAMAKER_ZEIN_DEFAULT_J`,
`HAMAKER_RANGE_PROTEIN_J`). Override with `--hamaker <value>` on the
command line, or run the sweep at both ends of the range to see how
sensitive your barrier predictions are to this choice -- it is the single
largest source of uncertainty in the model.

## Key assumptions and limitations (read before presenting results)

- **Geometry**: two equal, perfectly spherical particles of radius `a`
  (the DLS z-average radius). Real zein-caseinate NPs are not perfectly
  spherical or monodisperse; treat `a` as an effective radius.
- **Constant-potential boundary condition** for V_EDL (as opposed to
  constant-charge). The two bracket the true charge-regulated behavior;
  they mostly agree at the barrier position for kappa*a >> 1, which holds
  here.
- **Non-retarded van der Waals**: no zero-frequency retardation
  correction. Adequate at the D ~ 1-10 nm separations where the primary
  barrier sits for these particle sizes; would need revisiting for much
  larger separations.
- **Single, system-wide Hamaker constant**: does not capture protein
  layer heterogeneity or short-range "hydration force" effects.
  Molina-Bolivar, Galisteo-Gonzalez & Hidalgo-Alvarez (2001) (memo
  Section 5) document exactly this kind of anomalous high-ionic-strength
  stability for protein-coated colloids -- if your measured stability at
  high I disagrees with this model's prediction, that hydration-force
  physics (not a bad Hamaker constant) is the most likely explanation.
  Don't force-fit A to hide the discrepancy; report it as a known model
  limitation.
- **Sigmoid zeta(pH) model**: the 3-pH-point-per-ionic-strength design
  (Section 2 of the memo) cannot independently fit all 4 sigmoid
  parameters (zeta at low pH, zeta at high pH, pH at the isoelectric
  point, transition slope) per ionic-strength group -- that's
  under-determined. This pipeline fits the shape (pH_iep, slope) once,
  pooled across all ionic strengths (assumes the protein's intrinsic
  protonation chemistry doesn't change much with simple 1:1 salt
  screening -- a standard, defensible assumption), then refits only the
  2 remaining, genuinely ionic-strength-dependent parameters per group.
  Revisit this if your real data show the IEP itself shifting with
  ionic strength (documented for some systems at high multivalent-ion
  concentrations, e.g. the Ca2+ effects reported by Molina-Bolivar et
  al. 2001) -- the pooled-shape assumption would then need to be relaxed.
- **Particle radius held constant** across the pH x ionic-strength sweep
  (uses the mean measured radius). If your DLS matrix shows radius
  changing substantially with pH/I even at the initial/pre-aggregation
  timepoint, extend `sweep_stability_map()` to interpolate radius the
  same way zeta is interpolated.
- **Temperature**: fixed at 298.15 K (25 degC) by default throughout.
  For the pasteurization-relevant 63 degC point (Section 2), re-run with
  `temperature_K=336.15` passed through -- note this also technically
  requires adjusting `epsilon_r` (water's permittivity drops with
  temperature); the current code exposes `epsilon_r` as a keyword
  argument on every physics function for exactly this reason but does
  not auto-correct it, since accurate epsilon_r(T) requires a lookup
  table or empirical fit not yet included here.

## Extending this script

- **Fit your own Hamaker constant** from measured critical coagulation
  concentrations (CCC), the way Su et al. 2024 did (their eq. 7, using
  the Bjerrum length and Schulze-Hardy-rule scaling) -- a natural
  extension once you have a CCC titration series, not currently
  implemented here.
- **Ion-specific effects** (Hofmeister series, multivalent counterions):
  out of scope for the base DLVO model here; Su et al. 2024 is the
  relevant reference if you need to extend to this.
- **Secondary minimum analysis**: the current `find_energy_barrier()`
  only reports the primary energy barrier. If you need to characterize
  reversible (secondary-minimum) flocculation separately from
  irreversible primary aggregation, extend it to also report the local
  minimum at larger D.

## References (see also Section 5 of the memo)

- Hogg, R.; Healy, T. W.; Fuerstenau, D. W. Mutual coagulation of
  colloidal dispersions. *Trans. Faraday Soc.* **1966**, *62*, 1638-1651.
  (EDL repulsion formula.)
- Israelachvili, J. N. *Intermolecular and Surface Forces*, 3rd ed.;
  Academic Press, 2011. (Eq. 13.14 for sphere-sphere van der Waals; Eq.
  14.46 for sphere-sphere EDL repulsion.)
- Su, J.; et al. Electrolyte-induced aggregation of zein protein
  nanoparticles in aqueous dispersions. *J. Colloid Interface Sci.*
  **2024**, *656*, 2321-2329. (Zein-specific Hamaker constant, Grahame
  equation form used, IEP ~5.8.)
- Molina-Bolivar, J. A.; Galisteo-Gonzalez, F.; Hidalgo-Alvarez, R.
  Specific cation adsorption on protein-covered particles and its
  influence on colloidal stability. *Colloids Surf. B Biointerfaces*
  **2001**, *21*(1-3), 125-135. (Anomalous high-ionic-strength stability
  / hydration-force limitation of classical DLVO.)
- Delahaije, R. J. B. M.; Wierenga, P. A.; van Nieuwenhuijzen, N. H.;
  Giuseppin, M. L. F.; Gruppen, H. Protein Concentration and
  Protein-Exposed Hydrophobicity as Dominant Parameters Determining the
  Flocculation of Protein-Stabilized Oil-in-Water Emulsions. *Langmuir*
  **2013**, *29*(37), 11567-11574.
- Ravindran, S.; Williams, M. A. K.; Ward, R. L.; Gillies, G.
  Understanding how the properties of whey protein stabilized emulsions
  depend on pH, ionic strength and calcium concentration, by mapping
  environmental conditions to zeta potential. *Food Hydrocolloids*
  **2018**, *79*, 572-578.
