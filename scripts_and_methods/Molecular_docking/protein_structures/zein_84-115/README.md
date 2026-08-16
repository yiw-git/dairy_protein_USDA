> ## ⏹️ SUPERSEDED AS A DOCKING RECEPTOR
>
> This folder was created to serve as the zein receptor for HADDOCK/MEGADOCK.
> **That use was retired.** The justification below — five ColabFold models
> converging on this region at 0.3–1.1 Å — does not hold up:
>
> | Check | Value |
> |---|---|
> | Secondary structure of 84–115 | **100% α-helix** (φ/ψ SD 2.7°/3.2°) |
> | RMSD to a *generic* ideal poly-Ala helix | **0.59 Å mean** |
> | Cross-model RMSD (the original argument) | 0.3–1.1 Å |
>
> Same magnitude. Any two α-helices of equal length superimpose to ~0.5 Å
> regardless of sequence, so the agreement is **secondary-structure**
> confidence and carries no tertiary-structure information. A bare helix is a
> smooth cylinder — no pocket, no groove — leaving shape-complementarity
> docking rotationally and translationally degenerate.
>
> The files remain valid and correctly built; only the *interpretation* changed.
> See `../../DOCKING_STATUS_AND_NEXT_STEPS.md` §5b and §6.

---

# zein_84-115 — the high-confidence α-zein docking receptor

**File:** `zein_helix_84_115.pdb`
**Content:** α-zein residues **84–115** (32 residues), chain A, original
numbering from the parent model preserved.
**Sequence:** `LPLVHLLAQNIRAQQLQQLVLANLAAYSQQQQ`
— a single α-helix, the canonical Leu/Gln-rich α-zein repeat.

## Provenance

Extracted from `../zein_model.pdb`, which is byte-identical (RMSD 0.000 Å) to
`../zein_alpha_22597/..._relaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb`
(ColabFold, UniProt P02859, *Zea mays* α-zein).

## Why this region, and not the full 234-residue model

**1. The full-length fold is not reproducible.** Pairwise Cα RMSD between the
five ColabFold models, after optimal superposition:

|  | r001 | r002 | r003 | r004 | r005 |
|---|---|---|---|---|---|
| **r001** | 0.0 | 19.1 | 18.9 | 13.7 | 21.5 |
| **r002** | 19.1 | 0.0 | 17.2 | 13.5 | 20.2 |
| **r003** | 18.9 | 17.2 | 0.0 | 13.7 | 18.8 |
| **r004** | 13.7 | 13.5 | 13.7 | 0.0 | 18.8 |
| **r005** | 21.5 | 20.2 | 18.8 | 18.8 | 0.0 |

13.7–21.5 Å is not "some uncertainty" — the five predictions are effectively
five different molecules. Any docking result on the full model is specific to
one arbitrary conformer.

**2. This helix IS reproducible.** Same calculation, residues 84–115 only:

|  | r001 | r002 | r003 | r004 | r005 |
|---|---|---|---|---|---|
| **r001** | 0.0 | 0.6 | 0.5 | 0.7 | 0.5 |
| **r002** | 0.6 | 0.0 | 0.5 | 0.5 | 1.1 |
| **r003** | 0.5 | 0.5 | 0.0 | 0.3 | 0.8 |
| **r004** | 0.7 | 0.5 | 0.3 | 0.0 | 1.0 |
| **r005** | 0.5 | 1.1 | 0.8 | 1.0 | 0.0 |

**0.3–1.1 Å across five independent predictions.** Despite total disagreement
about the global fold, all five converge here. This is the one part of α-zein
that AlphaFold genuinely resolves.

**3. Supporting evidence:** pLDDT 55–78 in this window (vs 49.0 mean for the
full model, with 67% of residues below 50). Unrestrained HDOCK runs also placed
41–53% of their interface contacts inside residues 3–17 and 84–115, which
together are only 20% of the protein — a ~2.5× enrichment, because these are
the only segments with real secondary structure for a globular partner to pack
against.

## How to use it

Receptor for HADDOCK 2.4. Active/passive residue lists and full run settings:
`../../HADDOCK_run/HADDOCK_submission_sheet.md`

Two compact polar patches on the same helical face, run separately:

| Patch | Active residues | Span |
|---|---|---|
| **PATCH-N** | `88,92,95` | 10.3 Å |
| **PATCH-C** | `106,110,113,114` | 12.3 Å |

Do **not** use all seven in one run — they span 37.9 Å, wider than a 54 Å
globular ligand can contact, which reproduces the failure mode of the original
40-active-residue attempt.

## Limitations to state when presenting

- Eisenberg hydrophobic moment **⟨μH⟩ = 0.067 per residue**, well below the
  ~0.35 threshold for a strongly amphipathic helix. The "this face points out
  of the nanoparticle" assignment is a working hypothesis, not a structural
  fact — Leu84, Leu99, Leu102 and Val103 sit on the same face as the chosen
  active residues.
- A 32-residue helix is a proxy for one repeat unit of a crowded, curved,
  multi-chain nanoparticle surface. Avidity and multivalency are not captured;
  that physics belongs to the coarse-grained MD step (Section 3C).
- The helix is excised from a low-confidence parent model. Its internal
  geometry is reproducible, but its orientation relative to the rest of the
  protein is not.
