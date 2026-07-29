# DLVO Method — Detailed Instructions, Parameter Sourcing, and DLS Quality Control

Companion to `README.md` and `dlvo_model.py` in this folder. Implements
Method **A. DLVO / extended-DLVO modeling** from *Dry-Lab Strategy for
Section 3.1.3 (v3)* for the zein–caseinate NP (bare and after WPI /
casein exposure).

This document answers four things you asked:

1. What the **core equation** is and what **every parameter** in it means.
2. Which parameters you **must measure in the wet lab**, and which you
   can take from **literature / estimated values**.
3. What to do when the **PDI is high (e.g. 0.26)** and the size
   distribution shows **two peaks** — concrete solutions.
4. A **reference / publication** for each of the above.

Throughout, "the memo" = *Dry-Lab Strategy for Section 3.1.3 (v3)*.

---

## Part 1 — The core equation, term by term

### 1.1 The master equation

DLVO describes the **pair interaction energy** between two particles as a
function of their surface-to-surface separation `D`. It is the sum of two
competing terms:

```
V_total(D)  =  V_EDL(D)  +  V_vdW(D)
              (repulsive)   (attractive)
```

- **V_EDL(D)** — electrostatic double-layer repulsion. Two like-charged
  particles repel because their diffuse ion clouds (double layers)
  overlap. This is what keeps the dispersion stable.
- **V_vdW(D)** — van der Waals attraction. Always present, always pulls
  particles together, dominates at very short separation.

The **shape** of `V_total(D)` is what matters (see figure logic below):

```
   V/kT
    |        ___  <- primary maximum  V_max  (the energy barrier)
    |       /   \
  0 |------/-----\------------------ D
    |     /       \___________
    |    /              secondary minimum (weak, reversible flocculation)
    |   |
    |   | <- primary minimum (deep, irreversible aggregation)
```

You compute `V_total(D)` for each condition, find the **height of the
primary barrier `V_max`** (in units of thermal energy `kT`), and read off
stability:

| Barrier height | Interpretation | Physical meaning |
|---|---|---|
| `V_max` > ~15 kT | **stable** | Particles rarely have enough thermal energy to cross the barrier; slow aggregation |
| `V_max` ≈ 10–15 kT | **marginal** | Borderline; flag for wet-lab check |
| `V_max` < ~10 kT | **unstable** | Barrier easily crossed; fast aggregation expected |

The 10–15 kT rule of thumb is standard in the colloid literature
(Israelachvili, *Intermolecular and Surface Forces*, 3rd ed., Ch. 14).

**Foundational references (the theory itself):**
- Derjaguin, B.; Landau, L. Theory of the stability of strongly charged
  lyophobic sols. *Acta Physicochim. URSS* **1941**, *14*, 633–662.
- Verwey, E. J. W.; Overbeek, J. Th. G. *Theory of the Stability of
  Lyophobic Colloids*; Elsevier: Amsterdam, 1948.
- Israelachvili, J. N. *Intermolecular and Surface Forces*, 3rd ed.;
  Academic Press, 2011 (Ch. 13–14 — the modern working reference; source
  of the exact sphere–sphere formulas used in `dlvo_model.py`).

---

### 1.2 The van der Waals term

```
V_vdW(D) = -(A/6) · [ 2a²/(D(D+4a)) + 2a²/(D+2a)² + ln( D(D+4a)/(D+2a)² ) ]
```

This is the **exact (non-retarded) Hamaker sphere–sphere expression**,
valid at all separations (not only D ≪ a). Parameters:

| Symbol | Meaning | Units | Where it comes from |
|---|---|---|---|
| `A` | **Hamaker constant** — sets the strength of the vdW attraction; a property of the particle material across the medium (water) | J | **Literature** (see Part 2) |
| `a` | Particle **radius** | m (nm in code) | **Wet lab** — DLS |
| `D` | Surface-to-surface separation (the independent variable you sweep) | m (nm) | Swept by the code, not measured |

Reference for the formula: Israelachvili 2011, Eq. 13.14. The Hamaker
constant concept: Hamaker, H. C. The London–van der Waals attraction
between spherical particles. *Physica* **1937**, *4*(10), 1058–1072.

> **Why `A` matters so much:** it is the single largest source of
> uncertainty in the whole model. Run the analysis at both ends of the
> plausible range (see Part 2) to bracket your barrier predictions.

---

### 1.3 The electrostatic double-layer term

```
V_EDL(D) = 2π · ε₀ · ε_r · a · ζ² · ln(1 + exp(−κ·D))
```

This is the **Hogg–Healy–Fuerstenau (HHF) constant-potential, weak-overlap
expression** combined with the Derjaguin approximation. It is valid when
`κa ≫ 1`, which holds for your system (a ≈ 50–150 nm vs. Debye length
κ⁻¹ ≈ 0.8–3 nm across 10–150 mM NaCl). Parameters:

| Symbol | Meaning | Units | Where it comes from |
|---|---|---|---|
| `ε₀` | Vacuum permittivity | F/m | **Constant** (8.854×10⁻¹²) |
| `ε_r` | Relative permittivity of water | – | **Literature** (~78.5 at 25 °C) |
| `a` | Particle radius | m | **Wet lab** — DLS |
| `ζ` | **Zeta potential** — the potential at the shear plane; proxy for surface charge | V (mV in) | **Wet lab** — electrophoretic light scattering |
| `κ` | Inverse Debye length (double-layer thickness⁻¹) | 1/m | **Derived** from ionic strength |
| `D` | Separation | m | Swept |

Reference: Hogg, R.; Healy, T. W.; Fuerstenau, D. W. Mutual coagulation of
colloidal dispersions. *Trans. Faraday Soc.* **1966**, *62*, 1638–1651.
Reproduced as Israelachvili 2011, Eq. 14.46.

> **Note on the ζ² dependence:** because V_EDL scales with **ζ squared**,
> the barrier is extremely sensitive to zeta potential. A drop from ±40 mV
> to ±20 mV cuts the repulsion to a quarter. This is why ζ is your single
> most important measured input, and why measuring it well (Part 3
> caveats apply) matters.

---

### 1.4 The two derived quantities feeding the EDL term

**(a) Debye length (double-layer thickness), κ⁻¹**, from ionic strength:

```
κ = sqrt( 2·N_A·e²·I_SI / (ε₀·ε_r·k_B·T) )        κ⁻¹ = 1/κ
```

| Symbol | Meaning | Units | Source |
|---|---|---|---|
| `N_A` | Avogadro's number | 1/mol | Constant |
| `e` | Elementary charge | C | Constant |
| `k_B` | Boltzmann constant | J/K | Constant |
| `T` | Absolute temperature | K | **Set by you** (buffer temp; 298.15 K default) |
| `I_SI` | Ionic strength | mol/m³ (= 1000 × mol/L) | **Wet lab / recipe** — your buffer + salt |

`κ⁻¹` is the thickness of the ionic screening layer: ~3 nm at 10 mM,
~1.4 nm at 50 mM, ~0.8 nm at 150 mM NaCl. Higher salt → thinner layer →
weaker/shorter-ranged repulsion → lower barrier → less stable. This is the
mechanistic origin of salt-induced aggregation. Reference: standard
Debye–Hückel theory; used in this exact form for zein NPs by Su et al.
2024 (below).

**(b) Grahame equation (ζ → surface charge density σ)**, if you want σ:

```
linear (|ζ| ≲ 50–60 mV):   σ = ε₀·ε_r·κ·ζ
full nonlinear (1:1 salt):  σ = sqrt(8·ε₀·ε_r·k_B·T·I_SI) · sinh(e·ζ / (2·k_B·T))
```

Use the nonlinear form when |ζ| exceeds ~50–60 mV (zein at pH 4 or 7
routinely does). Reference: Grahame, D. C. *Chem. Rev.* **1947**, *41*,
441–501; and Israelachvili 2011, Ch. 14.

---

### 1.5 The zeta-versus-pH sigmoid (how a few points cover the whole grid)

Because ζ vs. pH follows the protein's sigmoidal isoelectric-point (IEP)
curve, a handful of measured ζ values lets you interpolate the full curve
and sweep DLVO across a continuous pH × ionic-strength grid — including
food conditions you cannot test directly (skim milk, whole milk, yogurt,
bread-release):

```
ζ(pH) = ζ_high + (ζ_low − ζ_high) / (1 + exp((pH − pH_IEP)/slope))
```

`pH_IEP` is where ζ = 0 (zein ≈ 5.8–6.2; measure your own). This is a
descriptive fit, not new physics; the caveat on how many parameters your
3-pH-point grid can support is in the README ("Sigmoid zeta(pH) model"
limitation). Zein IEP reference: Su, J. et al. 2024 (below).

---

## Part 2 — Which parameters come from the wet lab vs. literature

This is the practical heart of your question. Split into three buckets.

### 2.1 MEASURE in the wet lab (per condition) — non-negotiable

| Parameter | Instrument / method | Notes | Reference for method |
|---|---|---|---|
| **Zeta potential, ζ** | Electrophoretic light scattering (ELS) — same Zetasizer-class instrument as DLS | Most important input (enters as ζ²). Measure at each pH × ionic strength. Report Smoluchowski-based ζ. | ISO 13099-2:2012; Bhattacharjee 2016 |
| **Particle radius, a** | DLS (z-average, then halve diameter → radius) | Use the pre-aggregation / initial-timepoint value. See Part 3 on PDI. | ISO 22412:2017 |
| **Ionic strength, I** | Known from your buffer + added NaCl recipe (not "measured" as such, but set and recorded by you) | Compute I = ½ Σ cᵢzᵢ². For food matrices (milk) estimate from composition. | Standard; memo Section 2 |
| **pH** | pH meter | Sets which point on the ζ(pH) sigmoid you are on. | — |
| **Temperature, T** | Recorded from the DLS/ELS cell | 25 °C default; add a 63 °C pasteurization point if time allows (memo). | — |

Design these as the **Section 2 grid**: 3 pH × 3 ionic strength × {bare,
protein-exposed} × 2–3 replicates. That grid is already laid out in
`anchor_data_template.csv`.

### 2.2 TAKE from literature / use an estimated value

| Parameter | Recommended value | Range / caveat | Reference |
|---|---|---|---|
| **Hamaker constant, A** | **1.0×10⁻²⁰ J** (zein-specific default) | Protein colloids span **5×10⁻²¹ – 1×10⁻²⁰ J**. Run both ends to bracket. | Su et al. 2024 (fitted for zein); memo Section 3A |
| **Relative permittivity, ε_r (water)** | **78.5** at 25 °C | Drops with T: ≈ 78.3 (25 °C) → ≈ 65.9 (63 °C). Use a value for your actual temperature. | Malmberg & Maryott 1956 |
| **Zein isoelectric point, pH_IEP** | ~5.8–6.2 (or fit from your own ζ data) | Prefer measuring; literature only as a prior / sigmoid shape constraint. | Su et al. 2024 |
| **Born-repulsion / minimum-approach cutoff, D_min** | 0.3 nm | Standard hard-wall cutoff to avoid the unphysical vdW divergence at D→0. | Common DLVO practice (see `dlvo_model.py`) |
| Physical constants (`N_A`, `e`, `k_B`, `ε₀`) | CODATA | Exact. | CODATA 2018 |

### 2.3 Optionally FIT from your own data (advanced)

| Parameter | How | When | Reference |
|---|---|---|---|
| **Hamaker constant A** | Regress measured critical coagulation concentrations (CCC) against surface charge, per Schulze–Hardy scaling | Only if you run a salt-titration CCC series | Su et al. 2024 (their eq. 7) |
| **ζ(pH) sigmoid params** | `curve_fit` on your measured ζ vs. pH | You already do this in `run_dlvo_analysis.py` | — |

**Bottom line:** the only truly "literature" number you cannot avoid is the
**Hamaker constant** (and ε_r, which is well-established). Everything that
drives the answer — ζ, a, I, pH, T — you measure or set. That is exactly
why DLVO works with a small N: each condition is an *anchor that
parameterizes physics*, not a training row.

---

## Part 3 — High PDI (e.g. 0.26) and a two-peak size distribution

You observed **PDI ≈ 0.26 with two peaks** in the size distribution.
Here is how to interpret it and, more importantly, what to do.

### 3.1 What PDI 0.26 actually means

PDI is a **dimensionless width** derived from the cumulants fit of the DLS
autocorrelation function (PDI = (σ/mean)², roughly the squared relative
width). Rough bands used across the DLS community:

| PDI | Interpretation |
|---|---|
| < 0.05 | Nearly monodisperse (only seen with standards) |
| 0.05 – 0.1 | Narrow / near-monodisperse |
| 0.1 – 0.2 | Moderately polydisperse — typical, acceptable for many NP systems |
| **0.2 – 0.4** | **Broadly polydisperse — often a sign of a second population or aggregates** (your 0.26 is here) |
| 0.4 – 0.7 | Very broad; cumulants marginally meaningful |
| > 0.7 | Too broad for reliable DLS; distribution algorithm unreliable |

Two key facts that make **0.26 + two peaks** a caution flag, not a
verdict:

1. **The cumulants PDI assumes a single-mode (log-normal-ish)
   distribution.** For a genuinely bimodal sample, the single "z-average +
   PDI" pair is a poor descriptor — the z-average sits *between* the two
   real populations and represents neither. (ISO 22412; ISO 13321.)
2. **DLS scattered intensity scales ~ d⁶** (Rayleigh). A tiny number of
   large particles/aggregates scatters enormously more than the primary
   population, so an intensity-weighted second peak at large size can
   correspond to a *very small mass/number fraction*. A 1% aggregate
   population can dominate the intensity plot. (Malvern "Intensity–Volume–
   Number"; Bhattacharjee 2016.)

So the practical question is not "is 0.26 good or bad" but **"is the
second peak real, and is it a small tail or a genuine second population?"**

### 3.2 Diagnostic decision tree (do these first)

1. **Look at the intensity vs. volume vs. number distributions.**
   - If the large peak is prominent in **intensity** but nearly vanishes
     in **volume/number**, it is a **minor large-size tail** (a few big
     aggregates/dust) — the primary population is fine.
   - If it remains substantial in **volume**, you have a genuine second
     population that must be dealt with.
   - Caution: intensity→number transforms amplify noise; treat number
     distributions from a broad sample skeptically. (Malvern; Brookhaven
     FAQ.)

2. **Check peak separation.** DLS cannot resolve two populations unless
   they differ in size by roughly a factor of ~3 (≥5–6× to be robust at
   90°). If your two peaks are closer than ~3×, the "two peaks" may be an
   artifact of the fitting algorithm, not two real populations. (Wei et
   al. 2018; Malvern.)

3. **Check the correlogram / fit residuals and the intercept.** A clean
   single-exponential decay with a good intercept (~0.9–1.0) and low
   residuals argues for one population; a long tail / poor baseline argues
   for large aggregates or dust. Also inspect the **count rate stability**
   — spikes indicate transient large scatterers (dust/bubbles).

4. **Re-measure.** Run 3+ acquisitions. Dust and bubbles are transient and
   irreproducible; a real second population is reproducible.

### 3.3 Solutions — sample-side (fix or clean the sample)

- **Filter** the dispersion through a 0.22 µm (or 0.1 µm / 0.02 µm Anotop)
  syringe filter to remove dust and large aggregates before measuring; if
  the second peak disappears and the primary size is unchanged, it was
  contamination, not your particles. (Malvern DLS tips; Particle
  Technology Labs.)
- **Centrifuge** briefly (low speed) to pellet large aggregates, then
  measure the supernatant. Compare before/after to quantify the aggregate
  fraction. (Malvern; ISO 22412 sample-prep guidance.)
- **Degas / avoid bubbles**, use clean cuvettes, and equilibrate to
  temperature to kill convection artifacts.
- **Dilute appropriately** in the *same* ionic-strength/pH buffer (not
  water) to avoid multiple scattering *and* to avoid changing the
  colloidal chemistry you are trying to measure.
- **Optimize formulation** if the second peak is real aggregation: it may
  reflect that the condition is genuinely near the IEP or above the
  critical coagulation concentration — which is itself a **DLVO-relevant
  result**, not just a nuisance. Note the condition and where it sits on
  your stability map.

### 3.4 Solutions — analysis-side (report it honestly)

- **Do not report a single z-average + PDI for a truly bimodal sample.**
  Report the **peak table** (peak 1 / peak 2 mean size + % intensity),
  and state that the cumulants z-average is not representative. (ISO
  22412.)
- **For the DLVO model specifically:** use the **primary (smaller) peak's
  radius** as `a` for the intact NP population, and treat the large peak
  as evidence of partial aggregation at that condition rather than folding
  it into `a`. Feeding a z-average inflated by a few aggregates into DLVO
  would understate the true barrier. Add a note/flag in your
  `anchor_data_measured.csv` (`notes` column) for any bimodal condition.
- **Run a sensitivity check:** compute the DLVO barrier using both the
  small-peak radius and the z-average radius to show how much your
  stability call depends on the choice.

### 3.5 Solutions — orthogonal methods (confirm the second population)

Because DLS is intensity-weighted and low-resolution, confirm any
suspected second population with a technique that counts particles or
resolves them:

- **NTA (Nanoparticle Tracking Analysis)** — tracks individual particles;
  far better at resolving/quantifying polydisperse and bimodal samples,
  and gives number concentration. (Filipe, Hawe & Jiskoot 2010.)
- **Electron microscopy (TEM/SEM/cryo-TEM)** — direct sizing of the two
  populations; number-weighted, complements DLS. (Eaton et al. 2017.)
- **Multi-angle DLS / MADLS** — measuring at more than 90° improves
  bimodal resolution and reduces angle-dependent bias. (Wei et al. 2018.)
- **DLS + these together** is the standard way to establish whether a
  DLS "second peak" is real; a single technique is not enough for a
  publishable claim about polydispersity. (Bhattacharjee 2016; Varenne et
  al. 2016.)

### 3.6 One-paragraph summary you can paste into a report

> A PDI of 0.26 indicates a broadly polydisperse sample and, combined with
> a visibly bimodal intensity distribution, suggests either a minor
> population of large aggregates/dust or a genuine second size population.
> Because DLS intensity scales with the sixth power of diameter, a large
> intensity peak can represent a very small number fraction; we therefore
> examined intensity, volume, and number distributions, verified peak
> separation exceeded the ~3× DLS resolution limit, filtered/centrifuged to
> test whether the large peak was removable contamination, and confirmed
> the primary population by [NTA/TEM]. For DLVO modeling we used the
> primary-peak radius as the particle radius and flagged bimodal conditions
> as partially aggregated rather than averaging the two populations into a
> single z-average.

---

## Consolidated reference list

**DLVO theory & the equations**
- Derjaguin, B.; Landau, L. Theory of the stability of strongly charged
  lyophobic sols. *Acta Physicochim. URSS* **1941**, *14*, 633–662.
- Verwey, E. J. W.; Overbeek, J. Th. G. *Theory of the Stability of
  Lyophobic Colloids*; Elsevier, 1948.
- Israelachvili, J. N. *Intermolecular and Surface Forces*, 3rd ed.;
  Academic Press, 2011. (Sphere–sphere vdW Eq. 13.14; EDL Eq. 14.46;
  10–15 kT stability criterion.)
- Hamaker, H. C. The London–van der Waals attraction between spherical
  particles. *Physica* **1937**, *4*(10), 1058–1072.
- Hogg, R.; Healy, T. W.; Fuerstenau, D. W. Mutual coagulation of
  colloidal dispersions. *Trans. Faraday Soc.* **1966**, *62*, 1638–1651.
- Grahame, D. C. The electrical double layer and the theory of
  electrocapillarity. *Chem. Rev.* **1947**, *41*(3), 441–501.

**Parameter values (Hamaker, permittivity, zein IEP)**
- Su, J.; et al. Electrolyte-induced aggregation of zein protein
  nanoparticles in aqueous dispersions. *J. Colloid Interface Sci.*
  **2024**, *656*, 2321–2329. (Zein Hamaker ≈ 1.0×10⁻²⁰ J; Grahame form;
  IEP ~5.8; DLVO CCC prediction.)
- Malmberg, C. G.; Maryott, A. A. Dielectric constant of water from 0° to
  100 °C. *J. Res. Natl. Bur. Stand.* **1956**, *56*(1), 1–8.
- Molina-Bolívar, J. A.; Galisteo-González, F.; Hidalgo-Álvarez, R.
  Specific cation adsorption on protein-covered particles and its
  influence on colloidal stability. *Colloids Surf. B Biointerfaces*
  **2001**, *21*(1–3), 125–135. (Where classical DLVO breaks down at high
  ionic strength — the extended-DLVO / hydration-force caveat.)

**DLS, PDI, and bimodal / multimodal distributions**
- ISO 22412:2017. *Particle size analysis — Dynamic light scattering
  (DLS).* International Organization for Standardization. (Cumulants,
  z-average, PDI, reporting requirements.)
- ISO 13321:1996. *Particle size analysis — Photon correlation
  spectroscopy.* (Original cumulants method.)
- Malvern Panalytical. *Dynamic Light Scattering: Common Terms Defined*
  (whitepaper WP111214). (PDI interpretation bands.)
- Malvern Panalytical. *Intensity–Volume–Number: which size is correct?*
  (Knowledge-center insight; d⁶ intensity weighting and transform
  cautions.)
- Stetefeld, J.; McKenna, S. A.; Patel, T. R. Dynamic light scattering: a
  practical guide and applications in biomedical sciences. *Biophys. Rev.*
  **2016**, *8*, 409–427. (PDI, polydispersity, practical guidance.)
- Bhattacharjee, S. DLS and zeta potential — what they are and what they
  are not? *J. Control. Release* **2016**, *235*, 337–351. (Intensity
  weighting, limits of DLS, need for orthogonal methods.)
- Wei, Y.; et al. Effects of angular dependency of particulate light
  scattering intensity on determination of samples with bimodal size
  distributions using DLS. *Int. J. Mol. Sci.* **2018**, *19*(11), 3247.
  (Bimodal resolution limits; multi-angle benefit.)
- Filipe, V.; Hawe, A.; Jiskoot, W. Critical evaluation of Nanoparticle
  Tracking Analysis (NTA) by NanoSight for the measurement of
  nanoparticles and protein aggregates. *Pharm. Res.* **2010**, *27*,
  796–810. (NTA for polydisperse/bimodal confirmation.)
- Eaton, P.; et al. A direct comparison of experimental methods to measure
  dimensions of synthetic nanoparticles. *Ultramicroscopy* **2017**,
  *182*, 179–190. (DLS vs. EM sizing.)
- Varenne, F.; et al. Standardization and validation of a protocol of size
  measurements by DLS for monodisperse and polydisperse samples. *Talanta*
  **2016**, *160*, 194–203. (Multi-technique validation of polydispersity.)
- Panalytical / Particle Technology Labs. *Dynamic Light Scattering:
  Unexpected Results and False Positives.* (Filtration/centrifugation and
  dust-artifact troubleshooting.)

---

*See `README.md` for the software workflow and `dlvo_model.py` docstrings
for the exact coded form of every equation above.*
