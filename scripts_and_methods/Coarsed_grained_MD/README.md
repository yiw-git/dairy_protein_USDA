# Coarse-Grained MD — Section 3.1.3 (Method 3C)

A step-by-step operational guide to the coarse-grained (CG) molecular dynamics leg
of the Section 3.1.3 dry-lab strategy, written for someone who has already run
Method 3B (protein–protein docking) and needs a concrete build/run/analyze plan for
Method 3C. It implements **"C. Coarse-grained MD (MARTINI 3 + GROMACS)"** from
*Dry-Lab Strategy for Section 3.1.3 — v3*, for the zein-caseinate nanoparticle (NP)
surface interacting with whey (β-lactoglobulin, α-lactalbumin) and casein.

**Status of this folder:** guidance only — no run scripts yet. This README is meant
to be the reference you build the scripts against, mirroring how
[`../Molecular_docking/README.md`](../Molecular_docking/README.md) and
[`../DLVO/README.md`](../DLVO/README.md) work for their methods.

---

## 1. What this method adds that 3A/3B don't (30-second version)

- **3A (DLVO)** is a mechanistic *energy* calculation from bulk parameters (ζ, κ, a).
  It tells you *whether* a barrier to aggregation exists, not *which residues* are
  involved.
- **3B (docking)** gives you a single best-guess *static pose* — which protein binds
  best, and which residues touch. It is a hypothesis generator, not a dynamic answer.
- **3C (this method)** takes the docking output and asks the dynamic question 3B
  cannot: *does that interface actually hold together when both molecules are
  jostled by thermal motion in explicit solvent, at a given pH / ionic strength /
  temperature, over time?* It gives you a trajectory, not a single snapshot —
  contact persistence, aggregation behavior, and (optionally) a binding free energy.

Coarse-graining (MARTINI: roughly 4 heavy atoms → 1 bead) is what makes this
affordable. It smooths the energy landscape, which both lets you use a much larger
integration timestep (~20 fs vs. ~2 fs atomistic) and reduces particle count
severalfold — together this is why hundreds of ns–µs of CG sampling is realistic on
a single GPU workstation, where the equivalent atomistic run would not be.

**What it is not:** a substitute for the DLS wet-lab window's error bars, and — per
the memo's honesty note — not a load-bearing deliverable for a 4–6 week timeline.
Treat "1–2 representative conditions, a clean trajectory, a preliminary trend" as a
full success for this method.

---

## 2. Tools you need

| Tool | Role |
|---|---|
| **GROMACS** ≥ 2023.x | MD engine. No special CG build needed — CG just means different force-field parameters and a bigger timestep. |
| **martinize2** (from the `vermouth` Python package) | Converts an atomistic PDB → CG structure + topology. This is the actual "coarse-graining" step. |
| **Martini 3** force field files (`martini_v3.0.0.itp` + water/ion `.itp`s) | Use v3, not the older Martini 2.2 — v3 is the current, actively maintained parameter set and the only one with disordered-protein (IDP) support (relevant for the casein fragment, Section 4 below). |
| **DSSP** (`mkdssp`) | Secondary-structure assignment that martinize2 needs to place backbone dihedrals/elastic bonds correctly. |
| **insane.py** or `gmx solvate` + `gmx genion` | Builds the solvated, ionized CG box. `insane.py` is normally used for membrane systems, but it also works fine for a plain protein-in-water CG box (`-sol W`, no lipid) and is more robust for CG water/ion placement than the plain GROMACS tools — either works. |

Reference tutorials: the Martini Force Field Initiative's official Martini 3
protein tutorials ([cgmartini.nl/docs/tutorials/Martini3/ProteinsI](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/))
are the best up-to-date walkthrough of martinize2 options and are worth having open
while you build the first system.

---

## 3. Step-by-step workflow

### Step 0 — Carry inputs over from 3B (docking)

- Start from the HDOCK best-pose complexes already produced in
  `../Molecular_docking/HADDOCK/` (`D1_HDOCK.pdb` = zein–β-lactoglobulin,
  `D2_HDOCK.pdb` = zein–α-lactalbumin). You still need to generate an equivalent
  docked pose for a casein peptide fragment before starting 3C for casein (see
  `../Molecular_docking/README.md` for the casein-fragment-selection guidance).
- Note the docking-predicted **interface residues** for each complex — these become
  your sanity-check target in Step 5: does the CG trajectory keep those residues in
  contact, or does the interface drift/dissolve?
- Per the memo, do **not** simulate one zein chain + one milk protein 1:1. Build a
  small multi-chain **patch** (e.g., 4–9 zein copies arranged to approximate a flat
  or gently curved local NP surface, plus 2–4 copies of the target milk protein
  placed above it) — this is a partial, cheap fix for the avidity/multivalency
  concern flagged in 3B, without attempting to simulate the full ~100 nm particle
  (which is not tractable even at CG resolution).

### Step 1 — Atomistic → CG conversion (martinize2)

Run martinize2 once per unique protein/fragment:

```bash
martinize2 -f zein_model.pdb -x zein_cg.pdb -o zein_topol.top \
    -dssp mkdssp -p backbone -ff martini3001 \
    -elastic -ef 700 -el 0.5 -eu 0.9
```

- `-elastic` (classic elastic network, ENM) is the simplest structural bias and the
  right first pass to get the pipeline running end-to-end.
- For a more realistic run once the pipeline works, switch to **GōMartini**
  (`-go`, using a contact map derived from the atomistic structure) instead of a
  plain elastic network. Community guidance as of the current Martini 3 protein
  tutorials treats GōMartini as the default choice going forward, because its
  contacts can break and reform under strain/heat, unlike the rigid, always-on
  harmonic bonds of a classic elastic network — directly relevant to the memo's
  temperature-sweep caveat (see Section 4 below).
- **Casein fragment — do not use a rigid elastic network across the whole chain.**
  Casein is intrinsically disordered; standard Martini 3 protein parameters are
  known to over-compact IDP regions. Use martinize2's `-id-regions` flag to mark
  the disordered residues, combined with `-water-bias`/`-idr-tune`, or use the
  IDP-tuned bonded parameters from the Martini3-IDP branch if available in your
  martinize2/vermouth version. Skipping this step is the single most likely source
  of an artificially over-compact, unrealistic casein conformation.

**pH handling (do this *before* martinize2, not after):** MARTINI has no dynamic
protonation. Decide each ionizable residue's charge state (Asp/Glu/Lys/His) at your
target pH via Henderson–Hasselbalch using standard pKa values (Asp ≈ 3.9, Glu ≈ 4.3,
His ≈ 6.0, Lys ≈ 10.5), or run PROPKA on the atomistic structure first for
structure-specific pKa shifts. Rename the relevant residues in the input PDB
(e.g., neutral Asp/Glu, neutral His, neutral Lys) so the correct charge state is
what gets carried into the CG topology.

### Step 2 — Assemble the solvated, ionized box

```bash
gmx editconf -f complex_cg.pdb -o boxed.gro -d 1.5 -bt dodecahedron
gmx solvate -cp boxed.gro -cs water.gro -o solvated.gro -p system.top
gmx genion -s ions.tpr -o solvated_ions.gro -p system.top \
    -pname NA -nname CL -neutral -conc <target_M>
```

- **Ionic strength** is the easiest condition variable to control directly: convert
  your target mM to an ion count for the box volume and pass it to `gmx genion`
  (or `insane.py -salt`).
- Use standard non-polarizable Martini water (`W`) beads. If you'll run any
  condition at/below room temperature, add ~10% antifreeze beads (`WF`) —
  Martini's standard CG water is known to spuriously freeze in that range without
  them.
- Combine per-molecule topologies (zein × N copies, milk protein × M copies, water,
  ions) into one master `system.top` with correct `[molecules]` counts and includes
  for `martini_v3.0.0.itp`.

### Step 3 — Minimize → equilibrate → produce

Standard four-stage protocol, same order as an atomistic run, different `.mdp`
parameters:

1. **EM** — steepest descent, few thousand steps, until `Fmax` converges.
2. **NVT equilibration** — soft position restraints on the protein, v-rescale
   thermostat, short (~ns).
3. **NPT equilibration** — release restraints gradually, C-rescale (or
   Parrinello–Rahman) barostat, tens of ns.
4. **Production** — `dt = 0.02` (20 fs), v-rescale thermostat, C-rescale barostat.
   Realistic scope per the memo: hundreds of ns to a few µs of CG time for one
   condition within a 4–6 week window.

**Temperature** is a direct `ref_t` setting — compare room temperature (298 K) vs. a
pasteurization-relevant temperature (~336 K / 63 °C). Flag the same caveat as the
memo: with a classic elastic network, heating will *not* visibly unfold the protein,
because the harmonic restraints don't care about temperature. GōMartini partially
addresses this (contacts can break under thermal strain) but still isn't a
quantitative stand-in for real thermal denaturation — report any temperature
comparison as a qualitative probe, not a quantitative unfolding measurement.

### Step 4 — Analysis

| Question | Tool |
|---|---|
| Does the docking-predicted interface persist? | `gmx mindist` / `gmx pairdist` (or MDAnalysis contacts) — track contact-pair count over time; `gmx sasa` for buried interface area. |
| How many milk-protein copies stay attached to the zein patch? | `gmx clustsize`, or a short custom script counting protein copies within a contact-distance cutoff of the patch per frame. |
| Is the complex growing/shrinking with pH/IS/T? | `gmx gyrate` (radius of gyration of the whole complex) — the one metric that plugs directly into a cross-check against the DLS particle-size trend from the wet-lab anchor grid (Section 2 of the memo). |
| (Stretch goal) Binding free energy | `gmx pull` code to pull one milk-protein copy away from the patch along a COM–COM reaction coordinate; ~15–25 umbrella windows at 0.1–0.2 nm spacing, harmonic restraint per window, combine with `gmx wham` for a PMF curve → ΔG_bind. Budget this only after the core contact/Rg trajectories are done and clean. |

---

## 4. Known limitations — carry these into any write-up

- **Elastic-network rigidity.** A classic ENM keeps the protein essentially rigid
  regardless of temperature, so a temperature sweep with `-elastic` alone will not
  show real thermal unfolding. GōMartini is a better (but still approximate) option
  for this; the memo already flags this as a known limitation rather than something
  to paper over.
- **Martini 3 protein–protein interaction strength.** A 2024 study found that
  default Martini 3 parameters can misjudge binding for *flexible* proteins in
  solution and proposed a rescaling correction — directly relevant here, since your
  target proteins (especially casein) are not rigid globular folds. Treat CG-MD
  binding *strength* as qualitative/ranking information, not an absolute number, the
  same caution the memo applies to docking scores.
- **Casein disorder.** Standard Martini 3 protein parameters over-compact
  intrinsically disordered regions; a 2025 IDP-tuned parameter set (Martini3-IDP)
  exists specifically to correct this. Use IDP-aware settings (Step 1) rather than
  treating the casein fragment like a folded globular protein.
- **Single-fragment risk.** Casein's disorder means the "representative surface-active
  fragment" choice from 3B matters a lot. Where time allows, run ≥2 candidate
  fragments rather than betting the whole finding on one pick.
- **Patch, not particle.** Even a multi-chain zein patch is a simplification of the
  full curved ~100 nm NP surface — present it as a local-geometry approximation, not
  a whole-particle simulation.
- **Scope honesty.** This is, per the memo, the method with the steepest learning
  curve of the four. In a 4–6 week window, "1–2 conditions run cleanly with a
  preliminary contact/Rg trend" is the realistic target — a stretch goal to report
  alongside 3A/3B/4, not a load-bearing deliverable on its own.

---

## 5. Suggested build order (practical checklist)

1. **Pick one condition first** — e.g., pH 6.6, ~80 mM NaCl, 298 K (the "whole milk"
   condition already used in the ML manuscript inference, Section 4 of the memo).
   Running all pH × ionic-strength × temperature combinations at CG-MD resolution
   is not realistic in this window; one condition proves the pipeline and gives one
   real comparison point across all four methods.
2. **Convert D1 (zein–β-lactoglobulin) first**, using plain `-elastic`, since it
   already has a docked pose and is the simplest complex. Goal: get EM → NVT → NPT
   → a short production run to complete cleanly end-to-end before adding
   complexity.
3. **Add GōMartini and the casein fragment only after Step 2 works.** Bring in
   IDP-aware settings for casein at this point, not before — debugging force-field
   choices and system setup at the same time makes failures hard to diagnose.
4. **Only then** expand to a second condition (different T or pH) or attempt the
   PMF stretch goal from Step 4.

---

## 6. Key references

- Martini Force Field Initiative — official Martini 3 protein tutorials (martinize2
  options, elastic network vs. GōMartini vs. OLIVES structural-bias comparison):
  [cgmartini.nl/docs/tutorials/Martini3/ProteinsI](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/)
- Vermouth / martinize2 unified topology-generation framework: *eLife* 2024,
  [elifesciences.org/articles/90627](https://elifesciences.org/articles/90627)
- Rescaling protein–protein interactions for flexible proteins in Martini 3:
  *Nature Communications* 2024, [nature.com/articles/s41467-024-50647-9](https://www.nature.com/articles/s41467-024-50647-9)
- Martini3-IDP — improved Martini 3 parameters for disordered proteins (relevant to
  the casein fragment): *Nature Communications* 2025,
  [nature.com/articles/s41467-025-58199-2](https://www.nature.com/articles/s41467-025-58199-2)
- Fast calculation of protein–protein binding free energies via umbrella sampling
  with a coarse-grained model (PMF/WHAM workflow precedent):
  *J. Chem. Theory Comput.* 2018, [pubs.acs.org/doi/10.1021/acs.jctc.7b00660](https://pubs.acs.org/doi/10.1021/acs.jctc.7b00660)
- Coarse-grained MD of κ- and β-casein aggregates with curcumin (direct precedent on
  your protein class + a comparable polyphenol): *PLOS One* 2025 — already listed in
  the memo's Section 5 reference list.
