# ⏹️ SUPERSEDED — these 5 notebooks are prepared but were never run

Retained for provenance. They are complete and correct — the build bug
(`FileNotFoundError: ./megadock-gpu`, caused by the missing NVIDIA
`cuda-samples` clone) is fixed in all five, and every residue list is verified.
Only the **scientific justification** for running them was withdrawn.

## Why they were not run

These notebooks dock against `zein_helix_84_115.pdb`. That receptor was chosen
because the five ColabFold models of α-zein agree on residues 84–115 to
0.3–1.1 Å Cα RMSD, while disagreeing on the global fold by 13.7–21.5 Å.

That argument is circular:

| Check | Value |
|---|---|
| Secondary structure of region 84–115 | **100% α-helix** |
| φ/ψ spread | SD 2.7° / 3.2° |
| Rise per residue / radius | 1.45 Å / 2.32 Å (ideal: 1.50 / 2.30) |
| **RMSD to a generic ideal poly-Ala helix** | **0.59 Å mean (0.42–0.98)** |
| Cross-model RMSD (the original argument) | 0.3–1.1 Å |

The agreement with a *generic* helix is as good as the agreement *between
models*. Any two α-helices of equal length superimpose to ~0.5 Å regardless of
sequence. AlphaFold's confidence here is **secondary-structure** confidence and
carries no tertiary-structure information.

For rigid-body docking specifically this is fatal: a 32-residue helix is a
smooth cylinder with no pocket, groove or concavity, so shape-complementarity
scoring — which dominates MEGADOCK's function — is **rotationally degenerate**
about the helix axis and **translationally degenerate** along it. The score
landscape is flat, and the expected outcome is the same "no discrimination
between ligands" seen in the full-length run, for a new reason.

## If you ever do want to run them

They are ready. Upload to Drive and go:

| Notebook | Drive folder | Files to upload |
|---|---|---|
| `..._alpha_whey_...` | `zein84_115_whey` | helix PDB + `zein_helix_restricted.txt` + `alpha_whey_clean.pdb` |
| `..._beta_whey_...` | `zein84_115_whey` | + `beta_whey_clean.pdb` |
| `..._alphaS1_casein_...` | `zein84_115_casein` | + `alphaS1_casein_frag.pdb` |
| `..._beta_casein_...` | `zein84_115_casein` | + `beta_casein_frag.pdb` |
| `..._kappa_casein_...` | `zein84_115_casein` | + `kappa_casein_frag_compact_rank004.pdb` |

Receptor blocking keeps the 15-residue solvent-facing stripe
(`84,85,88,91,92,95,96,99,102,103,106,107,110,113,114`); `APPLY_BLOCK = False`
gives the unrestricted control.

The honest framing if you run them would be a **documented control**: evidence
that the helix is degenerate, rather than an argument that it must be. Look for
interface sizes staying at 3–6 residues and PATCH-N/PATCH-C calls splitting at
random — both would confirm the prediction above.

See `../../../DOCKING_STATUS_AND_NEXT_STEPS.md` §5b Finding 2 and §6.
