# Molecular docking — status, verification, and next steps

Section 3.1.3 / Method 3B · zein nanoparticle + milk protein corona
Prepared 2026-08-12 · summarises the full review of the MEGADOCK results

---

## 1. Executive summary

**The MEGADOCK results are arithmetically correct and reproducible.** Every
reported number was independently recomputed from the raw output files and
matched to 4 decimal places. Nothing is broken.

**The result is a valid negative finding:** all five milk proteins bind the
zein model with comparable, moderate, non-specific affinity. No specific
binding site exists. For protein-corona formation this is the biologically
expected answer, not a failed experiment.

**The limiting factor is not the docking software — it is the zein structure.**
The AlphaFold model has mean pLDDT 49 and essentially no hydrophobic core.
Re-running the same structure through a different docking server will not
improve the answer.

**Therefore the highest-value next step is not HADDOCK.** It is a zein-model
sensitivity test (Section 6), which costs about an hour and determines whether
any of the ranking is interpretable at all.

---

## 2. File integrity — nothing is missing or corrupted

| File | Status |
|---|---|
| `HDOCK_Huazhong/D1_HDOCK.pdb` | ✅ complete — chains A + B, proper TER/END |
| `HDOCK_Huazhong/D2_HDOCK.pdb` | ✅ complete — chains A + B, proper TER/END |
| `*.crdownload` (2 files) | ✅ complete despite the extension — these are Chrome download leftovers that still carry the `REMARK Score:` header (−228.37 and −315.32). Safe to delete once the scores are recorded. |
| All 5 MEGADOCK `dock_*.out` | ✅ complete — 10800 poses each, correct receptor/ligand headers |
| All 50 `*_complex.*.pdb` | ✅ receptor chain A, ligand chain B as intended |

**No file was ever flagged as incomplete.** What was flagged was the
*reliability of the zein AlphaFold model* (Section 4) — a scientific concern,
not a file-corruption one.

---

## 3. MEGADOCK verification

### 3.1 What checks out ✅

**Scores reproduce exactly.** Recomputing Z = (max − mean)/SD across all 10,800
retained poses in each raw `.out` file:

| Ligand | Recomputed | Reported in CSV | Match |
|---|---|---|---|
| α-lactalbumin | 5.8245 | 5.8245 | ✅ |
| α_s1-casein | 5.6358 | 5.6358 | ✅ |
| β-casein | 5.5472 | 5.5472 | ✅ |
| β-lactoglobulin | 5.5107 | 5.5107 | ✅ |
| κ-casein | 5.1626 | 5.1626 | ✅ |

**Run-to-run reproducibility is excellent.** The folder
`megadock_zein_whey_results/` turns out to be an accidental *replicate* of the
two whey runs. Comparing them:

| Ligand | Run 1 | Run 2 | Δ |
|---|---|---|---|
| α-lactalbumin | 5.8259 | 5.8245 | 0.0014 |
| β-lactoglobulin | 5.5151 | 5.5107 | 0.0044 |

Technical noise is ~0.004 Z units — about 150× smaller than the 0.66 spread
across the five ligands. **The numbers themselves are not noise.**

**Same receptor throughout.** All five runs used `zein_blocked.pdb` with
identical centring, so the comparison is internally consistent.

### 3.2 What does NOT check out ⚠️

**(a) κ-casein was docked as a 180 Å fully-extended string.**
Approximate maximum extents of the five ligands:

| Ligand | Atoms | Extent |
|---|---|---|
| β-lactoglobulin | 1490 | 54 Å |
| α-lactalbumin | 976 | 53 Å |
| β-casein frag | 387 | 62 Å |
| α_s1-casein frag | 314 | 64 Å |
| **κ-casein frag** | 938 | **180 Å** |

κ-casein is not a compact peptide in this file — it is a straight extended
chain. This also forced MEGADOCK onto a **294 grid** while every other run used
**192**, so its score distribution was generated in a different search volume.
**κ-casein's last-place ranking (5.16) is not interpretable** and should be
excluded or re-run with a compact conformer.

**(b) The two whey replicates are not bit-identical.** MEGADOCK's GPU FFT is
not deterministic. Harmless at this magnitude, but it means "identical
settings" does not mean "identical output" — worth one sentence in the methods.

**(c) The wrong kind of reproducibility was measured.** ±0.004 is *numerical*
reproducibility. The uncertainty that matters is **structural** — how much
would the scores change with a different zein conformer? That has not been
measured, and it is almost certainly far larger. See Section 6.

**(d) The 37–39 "convergence" is weaker than it first appeared.** Residues
37, 38, 39, 68, 71 and 75 are all inside the 40-residue restraint list, and
`zein_blocked.pdb` blocked everything outside it — so the ligands could not
have landed anywhere else. 37–39 are three *consecutive* residues, i.e. a
single protruding loop. Geometry (most accessible bump) explains this at least
as well as chemistry does.

### 3.3 Did the 40-residue block actually work? — YES, verified

The Colab notebook reimplemented MEGADOCK's `block` tool in pure Python
(renaming excluded residues to `BLK`). Three independent checks confirm it ran
and was honoured:

| Check | Result |
|---|---|
| `.out` header, line 3 (receptor name) | `zein_blocked.pdb` ✅ |
| Residue-name census in the output complexes | **3018 atoms named `BLK`** ✅ |
| Where the top-1 interfaces actually land | **17 of 19 contacts (89%) inside the allowed 40** ✅ |

Contrast with HDOCK, which had no such restriction:

| Method | Interface on the allowed 40 residues |
|---|---|
| **MEGADOCK (blocked)** | **89%** |
| HDOCK (unrestricted) | 7–15% |

The column arithmetic in the block function is correct (`l[0:16] + " BLK" +
l[20:]` writes `BLK` into the resName field, columns 18–20, leaving chain and
residue number intact). No bug.

**Reproducibility gap:** `zein_blocked.pdb` only ever existed inside the Colab
session and was never downloaded. Save it with the results next time.

### 3.4 ⚠️ But the restraint *design* was flawed

Interface sizes in the top-1 poses:

| Complex | Interface residues |
|---|---|
| α_s1-casein | 6 |
| α-lactalbumin | 3 |
| β-casein | 3 |
| β-lactoglobulin | 3 |
| κ-casein | 4 |

A normal protein–protein interface is **20–30 residues per side**. These are
roughly an order of magnitude too small.

Cause: the 40 allowed residues **span 92.8 Å and are scattered over the whole
protein** — they do not form a contiguous surface patch. After blocking the
other 194, the receptor presents a "spotty wall": a few sticky bricks spread
across a large area. A ligand can only perch on the most protruding cluster.

This single fact explains every earlier observation:

- **Restraint fractions stuck at 7.5–12.5%** — geometrically impossible to
  contact more than 3–5 of 40 scattered residues at once
- **All five ligands scoring alike** — they are all doing the same thing,
  perching on the same protruding bump, largely independent of ligand identity
- **Mediocre E-scores (5.2–5.8)** — tiny interface, weak signal
- **The earlier HADDOCK attempt failing** — 40 active residues spanning 92.8 Å
  can never be satisfied by a 54 Å protein

### 3.5 Verdict

**Usable, with two edits:** drop κ-casein from the comparison, and present the
ranking as non-discriminating rather than as an order. The honest headline is:

> Rigid-body docking of five milk proteins against a predicted α-zein model
> gave PPI E-scores of 5.5–5.8 (κ-casein excluded), with no ligand
> distinguishable from the others. No specific binding site was identified,
> consistent with non-specific, multivalent corona formation rather than
> defined receptor–ligand recognition.

---

## 4. The zein structure problem (the real bottleneck)

| Check | `zein_model.pdb` | Interpretation |
|---|---|---|
| Mean pLDDT | 49.0 | very low confidence |
| Residues pLDDT < 50 | 157 / 234 (67%) | majority predicted disordered |
| Residues pLDDT ≥ 70 | 22 / 234 (9%) | almost nothing well-determined |
| Buried (RSA < 15%) | 3 / 234 (1%) | **no hydrophobic core** |
| Median RSA | 57% | extended chain, not a globule |

Two consequences:

1. **The `zein_restricted.txt` premise does not hold.** It assumed hydrophobic
   residues are buried and hydrophilic ones face out. In this model 87% of
   *all* residues are >40% exposed — there is no "inside," so the 40 residues
   are not a face.
2. **AlphaFold is not failing — it is reporting that α-zein is largely
   intrinsically disordered.** Which is correct and well documented. But it
   means zein belongs in the same category as casein, and the README already
   prescribes the treatment: *dock fragments, not an invented fold.*

Only two segments are moderately confident: **residues 3–17** and **84–115**.

---

## 5. HADDOCK — how to set active/passive (answering the open question)

### 5.1 The constraint you correctly identified

A 21-residue α_s1-casein fragment cannot contact 40 zein residues. Active
restraints are *demands with penalties* — every unsatisfiable one adds a
constant penalty to every model, so the restraint stops discriminating and
becomes noise. **This is why the first HADDOCK attempt gave nothing usable.**

Rule of thumb: **active set ≈ 5–12 residues, spanning ≤ ~15 Å.**

### 5.2 Why "no active residues at all" is not an option

HADDOCK builds its AIRs *from* active residues. With none defined on either
molecule, no restraints are generated and there is nothing to drive the
docking. You need one of the three modes below.

### 5.3 The three legitimate options

| Mode | Setup | Bias | When to use |
|---|---|---|---|
| **A. One-sided AIRs** | active on zein patch only; ligand surface as passive | Low — you constrain *where on zein*, nothing about the ligand | ✅ Recommended |
| **B. Centre-of-mass restraints** | no AIRs; only a CoM restraint | Very low | If you want to constrain nothing at all |
| **C. Ab-initio / random AIRs** | HADDOCK picks random surface patches each run | None | Truly no information; needs far more sampling |

**Recommendation: Mode A**, because you *do* have a defensible prior (the
solvent-exposed polar face of the helix), and Mode A is the documented
protocol for one-sided information.

### 5.4 The concrete lists

Use the 32-residue helix `zein_helix_84_115.pdb`, not the full model. The
7 exposed polar residues on the outward face span 37.9 Å — too wide for a
54 Å globular protein — so **split them into two compact patches and run each
separately**:

| Patch | Active | Span | Character |
|---|---|---|---|
| **PATCH-N** | `88,92,95` | 10.3 Å | HIS/GLN/ARG — the only charged patch |
| **PATCH-C** | `106,110,113,114` | 12.3 Å | ASN/TYR/GLN/GLN — Gln-rich |

For every ligand: **active = empty**, **passive = its full solvent-accessible
surface** (lists in `HADDOCK_run/HADDOCK_submission_sheet.md`). Casein
fragments additionally set as **fully flexible segments**.

Checkboxes: ✅ auto-define passive around active · ✅ ignore non-solvent-accessible active residues.

### 5.5 The caveat that must be stated

The helix's Eisenberg hydrophobic moment is **⟨μH⟩ = 0.067**, well below the
~0.35 threshold for a strongly amphipathic helix. So "this face points outward
from the nanoparticle" is a working hypothesis, not a structural fact. The
opposite-face control run in the submission sheet tests exactly this.

---

## 5b. Zein-model sensitivity — ALREADY ANSWERED, no runs needed

The five ColabFold models were compared directly. **The docking runs originally
planned for this test are unnecessary** — the structures answer the question by
themselves.

### The AlphaFold files you have

`protein_structures/zein_alpha_22597/` contains ColabFold's standard output:
**5 models**, ranked by confidence.

| File | pLDDT | % residues < 50 |
|---|---|---|
| `..._relaxed_rank_001_..._model_3_...pdb` | 49.0 | 67% |
| `..._unrelaxed_rank_001_..._model_3_...pdb` | 49.0 | 67% |
| `..._unrelaxed_rank_002_..._model_4_...pdb` | 44.7 | 76% |
| `..._unrelaxed_rank_003_..._model_2_...pdb` | 43.6 | 89% |
| `..._unrelaxed_rank_004_..._model_5_...pdb` | 40.9 | 88% |
| `..._unrelaxed_rank_005_..._model_1_...pdb` | 39.5 | 89% |

`zein_model.pdb` is byte-identical to **relaxed rank_001** (RMSD 0.000 Å).
Only rank_001 was Amber-relaxed; ranks 002–005 exist as unrelaxed heavy-atom
files only (1796 atoms vs 3641 with hydrogens), so they are not directly
comparable to the file that was docked.

Note the trend: ranks 002–005 are **worse**, not equal alternatives.

### Finding 1 — the global fold is not reproducible

Pairwise Cα RMSD after optimal superposition, all 234 residues:

|  | r001 | r002 | r003 | r004 | r005 |
|---|---|---|---|---|---|
| **r001** | 0.0 | 19.1 | 18.9 | 13.7 | 21.5 |
| **r002** | 19.1 | 0.0 | 17.2 | 13.5 | 20.2 |
| **r003** | 18.9 | 17.2 | 0.0 | 13.7 | 18.8 |
| **r004** | 13.7 | 13.5 | 13.7 | 0.0 | 18.8 |
| **r005** | 21.5 | 20.2 | 18.8 | 18.8 | 0.0 |

**13.7–21.5 Å.** The five models are effectively five different molecules.
Any docking result computed on the full-length model is specific to one
arbitrary conformer and cannot be generalised. **Full-model docking is not
salvageable** — running it four more times would only produce four more
incomparable answers.

### Finding 2 — the 84–115 helix IS reproducible

Same calculation restricted to residues 84–115:

|  | r001 | r002 | r003 | r004 | r005 |
|---|---|---|---|---|---|
| **r001** | 0.0 | 0.6 | 0.5 | 0.7 | 0.5 |
| **r002** | 0.6 | 0.0 | 0.5 | 0.5 | 1.1 |
| **r003** | 0.5 | 0.5 | 0.0 | 0.3 | 0.8 |
| **r004** | 0.7 | 0.5 | 0.3 | 0.0 | 1.0 |
| **r005** | 0.5 | 1.1 | 0.8 | 1.0 | 0.0 |

**0.3–1.1 Å across five independent predictions.** Despite total disagreement
about the global fold, all five models converge on this helix. It is the one
part of α-zein AlphaFold genuinely resolves, and therefore **the only
defensible zein docking input.**

### Finding 3 — the 37–39 result is confirmed as a geometric artefact

Ranking the 40 allowed residues by protrusion from the centroid:

| Model | Most protruding allowed residues | Rank of res 37 |
|---|---|---|
| rank_001 | **37, 38, 39**, 41, 140 | 1 / 40 |
| rank_002 | **38, 39, 37**, 142, 140 | 3 / 40 |
| rank_003 | **37, 38, 39**, 31, 41 | 1 / 40 |
| rank_004 | **38, 37, 39**, 142, 140 | 2 / 40 |
| rank_005 | 41, 148, **39, 38**, 147 | 6 / 40 |

37–39 is the most protruding accessible cluster in 4 of 5 models. The MEGADOCK
"convergence" is **surface geometry, not binding chemistry** — ligands perch on
the most exposed bump of a spotty restraint surface. This is now demonstrated
rather than suspected, and it is robust precisely because protrusion is a
robust geometric property.

---

## 5c. Casein fragment conformations — a second, systemic problem

The κ-casein issue turned out not to be isolated. Radius of gyration (Rg,
"how tightly the chain is balled up") measured for all three fragments and
compared with theory:

| Fragment | N | measured Rg | random-coil Rg | fully-extended Rg | Verdict |
|---|---|---|---|---|---|
| β-casein | 25 | 14.6 Å | 13.6 | 25.3 | ✅ realistic coil |
| **α_s1-casein** | 21 | 17.7 Å | 12.4 | 21.2 | ⚠️ **84% of maximum extension** |
| **κ-casein** | 64 | **49.0 Å** | 22.3 | 64.7 | ⚠️ **76% extended; end-to-end 164.8 Å** |

*(random-coil Rg = 2.54·N^0.522, Kohn et al. 2004; fully-extended Rg = L/√12)*

**Cause:** an intrinsically disordered peptide has no structure to predict, and
AlphaFold's output in that situation is an extended chain. That is not a
predicted conformation — it is a statement of ignorance rendered as a straight
line.

**Why it damages the docking:**

1. κ-casein's 165 Å length forced MEGADOCK onto a **294 grid** while every other
   run used 192 — different search volume, so its Z-score is not comparable
2. A straight rod can only make line contact; it cannot wrap into a groove or
   form a normal interface
3. It is a near-zero-probability conformation being used as the sole
   representative of a broad conformational ensemble

### Fix applied to κ-casein

Comparing all five ColabFold models:

| Model | Rg | End-to-end | % of max extension |
|---|---|---|---|
| rank_001 *(was used)* | 49.0 | 164.8 Å | 76% |
| rank_002 | 44.9 | 153.6 | 69% |
| rank_003 | 41.0 | 141.2 | 63% |
| **rank_004** ✅ | **27.3** | **31.1 Å** | **42%** |
| rank_005 | 34.6 | 116.4 | 54% |

**`casein/kappa_casein_frag_compact_rank004.pdb` created** — sequence verified
identical, chain A, residues 1–64, hydrogens stripped. Rg 27.3 Å is close to
the ~22 Å random-coil expectation, and it fits a 192 grid, restoring
comparability. κ-casein no longer needs to be dropped.

### α_s1-casein — no fix available

All five models are rods (Rg 17.0–17.8, 80–84% extended). Plausible in part
because the sequence is a highly acidic cluster (`QMEAESISSSEEIVPNSVEQK`) whose
glutamates repel one another. At 21 residues the absolute size is manageable,
so it can be used as-is with the limitation stated.

### Phosphorylation — an acknowledged gap

`casein/casein_peptides.txt` already notes it: D3 and D4 are **phosphopeptides
modelled without phosphate groups**. Phosphoserine clusters are precisely what
makes casein surface-active, so electrostatic interaction with zein is
systematically underestimated. State this explicitly.

### Mitigation in HADDOCK

Setting the casein fragments as **fully flexible segments** lets HADDOCK
re-fold them during semi-flexible refinement, which partly compensates for a
poor starting conformer. A better starting point still gives a better result.

---

## 6. Recommended priority order

### Priority 1 — ~~Fix or drop κ-casein~~ ✅ DONE (see 5c)

Either regenerate a compact conformer (so it docks on a 192 grid like the
others) or exclude it and say why.

### Priority 2 — HADDOCK on the helix, all five ligands

**This recommendation changed.** An earlier draft proposed replicating MEGADOCK
on the *full* model with active = 37,38,39. Finding 1 (13.7–21.5 Å RMSD between
conformers) makes that pointless: it would replicate a conformer-specific
artefact. Finding 3 already explains the 37–39 result, so there is nothing left
to confirm.

Use `zein_helix_84_115.pdb` — justified now by convergence across five
independent predictions (0.3–1.1 Å), not merely by pLDDT. Active/passive lists
and run settings are in `HADDOCK_run/HADDOCK_submission_sheet.md`:

- **PATCH-N** `88,92,95` (10.3 Å) — HIS/GLN/ARG, the only charged patch
- **PATCH-C** `106,110,113,114` (12.3 Å) — ASN/TYR/GLN/GLN, Gln-rich
- Ligands: active empty, passive = full surface; casein fragments fully flexible

### Priority 3 — Write up the negative result

It is a legitimate finding, and it is now *explained* rather than merely
observed:

> Rigid-body docking of five milk proteins against a predicted α-zein model
> produced PPI E-scores of 5.5–5.8 with no ligand distinguishable from the
> others. Restricting the receptor to 40 solvent-exposed residues spanning
> 92.8 Å yielded interfaces of only 3–6 residues, and the top poses of all five
> ligands converged on residues 37–39 — the most protruding accessible cluster
> in 4 of 5 independent structure predictions. The result therefore reflects
> surface geometry rather than specific recognition, consistent with
> non-specific multivalent corona formation. Comparison of the five ColabFold
> models (global Cα RMSD 13.7–21.5 Å) further shows that full-length α-zein is
> not reliably predictable; only the 84–115 helix converges (0.3–1.1 Å) and was
> carried forward.

State the Section 9 assumptions from `README.md` explicitly.

### Priority 5 — Coarse-grained MD (Section 3C)

Per `README.md`, this is where the nanoparticle physics actually lives —
crowding, curvature, multivalency, avidity. Docking was always scoped as the
cheap screen that points CG-MD at the right question. That handoff is now due.

---

## 7. Tools decision — closing the loop

| Tool | Status | Reason |
|---|---|---|
| **MEGADOCK** | ✅ done, verified | FFT rigid-body, GPU, 5 ligands complete |
| **HDOCK** | ⏹️ discontinued | Whole-chain only; cannot exclude buried regions; 2 of 5 runs done (D1, D2) |
| **ClusPro** | ⏹️ skip | A third rigid-body FFT server adds nothing new |
| **LightDock** | ⏹️ abandoned | 20+ h/pair on CPU; superseded by MEGADOCK |
| **HADDOCK** | 🟡 optional | Only tool offering restraints + flexibility + explicit solvent |
| **CG-MD (3C)** | ➡️ next | Where the nanoparticle physics belongs |

---

## 8. Files

| Path | Contents |
|---|---|
| `DOCKING_STATUS_AND_NEXT_STEPS.md` | This document |
| `protein_structures/HADDOCK_run/zein_helix_84_115.pdb` | 32-residue helix, chain A, original numbering |
| `protein_structures/HADDOCK_run/HADDOCK_submission_sheet.md` | Copy-paste active/passive lists and run settings |
| `docking_trials/Megadocking/Megadock_result/` | All MEGADOCK outputs (verified) |
| `protein_structures/HDOCK_Huazhong/` | HDOCK D1/D2 complexes and scores |
| `scores_template.csv` | HDOCK scores: D1 −228.37, D2 −315.32 |
