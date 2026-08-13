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

### 3.3 Verdict

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

## 6. Recommended priority order

### Priority 1 — Zein-model sensitivity test (≈1 hour, do this first)

`protein_structures/zein_alpha_22597/` already contains **five** AlphaFold
models (`rank_001` … `rank_005`). Re-run MEGADOCK for **one** ligand (β-lg)
against ranks 002–005 with otherwise identical settings.

This answers the only question that currently matters:

- If E-scores across the five zein models spread **much more than 0.66**,
  then the ligand ranking is an artefact of one arbitrary conformer, and no
  amount of HADDOCK will fix it. You report the negative result and move on.
- If they spread **much less than 0.66**, the ranking survives structural
  uncertainty and becomes worth defending.

Either outcome is publishable methodology, and it costs four Colab runs.

### Priority 2 — Fix or drop κ-casein

Either regenerate a compact conformer (so it docks on a 192 grid like the
others) or exclude it and say why.

### Priority 3 — Write up the negative result

It is a legitimate finding. Frame it as: docking found no specific recognition
site, which supports a non-specific multivalent corona model and motivates the
CG-MD step. State the Section 9 assumptions from `README.md` explicitly.

### Priority 4 — HADDOCK (optional, low priority)

Worth **one or two runs**, not ten — as a methodological cross-check that adds
flexibility and explicit solvent, on the better-founded 32-residue helix. It is
**not** a replication of MEGADOCK: residues 37–39 are not even present in the
helix input, so the two ask different questions.

Only run all five ligands if Priority 1 shows the ranking is structurally
robust.

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
