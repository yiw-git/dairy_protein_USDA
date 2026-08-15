# HADDOCK 2.4 submission sheet — zein helix vs 5 milk proteins

Server: https://wenmr.science.uu.nl/haddock2.4/ (free, academic registration required — approval can take 1–2 days)

Run all 5 with **identical settings**. Only Molecule 2 changes.

---

## Why the zein input changed

The full AlphaFold model `zein_model.pdb` was rejected as a docking receptor because:

| Check | Result | Meaning |
|---|---|---|
| Mean pLDDT | **49.0** | AlphaFold has very low confidence overall |
| Residues pLDDT < 50 | 157 / 234 (**67%**) | Majority predicted as disordered |
| Residues pLDDT ≥ 70 | 22 / 234 (9%) | Almost no well-determined structure |
| Buried residues (RSA < 15%) | **3 / 234 (1%)** | No hydrophobic core — nothing is "inside" |
| Median RSA | 57% | Extended/open chain, not a globule |

Because there is effectively no buried core, the original premise behind
`zein_restricted.txt` ("hydrophobic residues are buried, hydrophilic ones face
out") does not hold in this model — 87% of *all* residues are >40% exposed.

Only two stretches are moderately confident: **residues 3–17** and **84–115**.
Region **84–115** was selected (longer, higher pLDDT 55–78).

**Input file:** `zein_helix_84_115.pdb` (32 residues, chain A, original numbering preserved)
**Sequence:** `LPLVHLLAQNIRAQQLQQLVLANLAAYSQQQQ` — a single α-helix, the canonical
Leu/Gln-rich α-zein repeat.

---

## Molecule 1 — zein helix (receptor)

- **File:** `zein_helix_84_115.pdb`
- **Chain:** A · **Residue range:** 84–115 (numbering kept from the parent model)

### ⚠️ REVISED — do not use all 7 residues in one run

The original 7-residue set (88, 92, 95, 106, 110, 113, 114) spans **37.9 Å**
along the helix. β-lactoglobulin and α-lactalbumin are only ~54 Å across, so a
single molecule cannot contact all seven at once. Demanding it does recreates
the exact failure of the 40-active attempt: permanently unsatisfiable
restraints that stop discriminating between poses.

**Split into two separate runs per ligand, each with a compact patch:**

| Patch | Active residues | Span | Composition |
|---|---|---|---|
| **PATCH-N** | `88,92,95` | 10.3 Å | HIS, GLN, ARG — the only charged patch on the helix |
| **PATCH-C** | `106,110,113,114` | 12.3 Å | ASN, TYR, GLN, GLN — the Gln-rich patch |

Both are polar, solvent-exposed (RSA 45–74%), and on the same helical face.

| Res | AA | RSA | Angle from outward face | Patch |
|---|---|---|---|---|
| 88 | HIS | 48% | 6° | N |
| 92 | GLN | 57% | 27° | N |
| 95 | ARG | 66% | 38° | N |
| 106 | ASN | 45% | 27° | C |
| 110 | TYR | 54% | 8° | C |
| 113 | GLN | 74% | 50° | C |
| 114 | GLN | 62% | 46° | C |

Running both patches is itself informative: if PATCH-N (charged) consistently
beats PATCH-C (Gln-rich), that is a mechanistic statement about electrostatics
driving corona formation — and it is testable against the pH/ionic-strength
work in 3A.

**Checkboxes:**
- ✅ Automatically define passive residues around the active residues
- ✅ Do not use active residues that are not solvent accessible

### ⚠️ Caveat to record in the methods

The Eisenberg hydrophobic moment of this helix is **⟨μH⟩ = 0.067 per residue**,
well below the ~0.35 threshold for a strongly amphipathic helix. So the
helix is only **weakly** amphipathic — the "this face points out of the
nanoparticle" assignment is a reasonable working hypothesis, **not** a
structurally determined fact. Leu84, Leu99, Leu102 and Val103 sit on the same
face as the chosen active residues.

State this explicitly rather than implying a clean amphipathic geometry.

---

## Molecule 2 — the milk protein (ligand), one per run

**Active residues: leave EMPTY for all five.** You have no prior knowledge of
which face of the milk protein binds; declaring one would manufacture the
answer. Defining only the partner's surface as passive is the documented
HADDOCK protocol for one-sided information.

**Uncheck** "automatically define passive residues" for Molecule 2 — paste the
lists below instead.

### D1 — β-lactoglobulin · `whey/beta_whey_clean.pdb` · chain A · 156 res (5–160)

Passive (94 surface residues, RSA ≥ 15%):

```
5,6,8,9,11,13,14,16,18,20,28,29,30,33,34,35,36,38,44,45,47,48,49,50,51,52,53,55,60,61,62,63,64,65,66,67,68,69,70,72,74,75,76,77,78,79,83,85,87,88,90,91,96,98,99,100,101,109,110,111,112,113,114,115,124,125,126,127,128,129,130,131,134,135,137,138,139,141,142,143,144,145,146,148,149,150,151,152,153,154,155,157,158,159
```

### D2 — α-lactalbumin · `whey/alpha_whey_clean.pdb` · chain A · 122 res (1–122)

Passive (75 surface residues):

```
1,2,4,5,6,7,9,10,11,13,14,16,17,18,19,20,22,31,32,33,35,37,39,41,42,43,44,45,46,49,58,59,62,64,65,66,67,68,70,72,74,75,76,78,79,81,82,83,84,86,87,90,93,94,97,98,99,100,101,102,103,105,108,109,110,112,113,114,115,116,117,119,120,121,122
```

### D3 — β-casein fragment · `casein/beta_casein_frag.pdb` · chain A · 25 res

Passive (all 25 — every residue is surface-exposed in a short peptide):

```
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
```

Also set **fully flexible segment: 1–25**.

### D4 — α_s1-casein fragment · `casein/alphaS1_casein_frag.pdb` · chain A · 21 res

Passive (all 21):

```
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21
```

Also set **fully flexible segment: 1–21**.

### D5 — κ-casein fragment · `casein/kappa_casein_frag_compact_rank004.pdb` · chain A · 64 res

> ⚠️ **Input file changed.** The original `kappa_casein_frag.pdb` (AlphaFold
> rank_001) is a nearly straight rod — Rg 49.0 Å, end-to-end 164.8 Å, i.e. 76%
> of maximum possible extension, versus ~22 Å expected for a disordered
> 64-residue random coil. That conformation forced MEGADOCK onto a 294 grid
> instead of 192 and is a physically improbable single snapshot.
>
> **AlphaFold rank_004** of the same sequence is far more realistic —
> **Rg 27.3 Å, end-to-end 31.1 Å** — and is the recommended input. Sequence
> verified identical; chain A, residues 1–64, hydrogens stripped.

Passive (all 64):

```
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64
```

Also set **fully flexible segment: 1–64**.

> The casein fragments are intrinsically disordered. Marking them fully
> flexible is the whole reason for using HADDOCK over MEGADOCK — it lets the
> peptide adapt its conformation instead of docking a frozen arbitrary snapshot.

---

## Sampling / run parameters

Keep the defaults unless you have a reason not to:

| Parameter | Value |
|---|---|
| Structures for rigid-body docking (it0) | 1000 |
| Structures for semi-flexible refinement (it1) | 200 |
| Structures for final water refinement (itw) | 200 |
| Clustering method | FCC |
| Clustering cutoff | 0.60 |
| Minimum cluster size | 4 |

---

## Recommended control run (do at least one)

**Restraint-sensitivity test.** Run one extra job — say β-lactoglobulin —
with the active set moved to the *opposite* helical face:

```
97,98,100,101,112,115
```

Then compare the two HADDOCK scores.

- If the opposite-face run scores **similarly**, the restraints are dictating
  the result and the docking has no independent discriminating power. This is
  a real possibility given ⟨μH⟩ = 0.067, and you want to find it out yourself
  rather than have a reviewer find it.
- If the chosen face scores **clearly better**, the face assignment is
  supported by the energetics and you can say so.

This single control converts "I restrained it and it bound there" into a
defensible statement.

---

## What to bring back for analysis

From the results page, for every run:

1. Cluster table: **cluster size**, **HADDOCK score ± SD**, **Z-score**
2. Energy decomposition per cluster: **E_vdW**, **E_elec**, **E_desolv**, **E_AIR**
3. Interface RMSD of each cluster from the overall best structure
4. The top-cluster representative PDB

Reminders when reading them:
- HADDOCK score: **more negative = better** (opposite sign to MEGADOCK's E-score)
- Z-score: **more negative = better**
- Overlapping SDs between clusters = **not distinguishable**
- Cluster size matters as much as score — a large cluster means convergence
- The score is arbitrary units, **not** kcal/mol

---

## Cross-validation scope vs MEGADOCK

| Independent between MEGADOCK and HADDOCK? | |
|---|---|
| Search algorithm | ✅ yes (FFT exhaustive vs restraint-guided MD) |
| Flexibility | ✅ yes (rigid vs semi-flexible) |
| Solvent | ✅ yes (none vs explicit water) |
| Scoring function | ✅ yes |
| **Zein structure and surface definition** | ❌ **no — and now they differ** |

Note the change: MEGADOCK docked the **full 234-residue model** restrained to 40
residues; HADDOCK will dock the **32-residue helix 84–115**. Residues 37–39, the
patch that dominated the MEGADOCK results, are **not present** in the HADDOCK
input. So these are **not** two tests of the same hypothesis, and the two
rankings should be compared with that stated plainly — HADDOCK is a fresh,
better-founded run, not a replication of MEGADOCK.
