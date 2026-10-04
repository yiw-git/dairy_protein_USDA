# Glucose--alpha-lactalbumin pilot inputs

## Downloaded source structures

| File | Source | Purpose |
| --- | --- | --- |
| `human_alpha_lactalbumin_P00709_1B9O.pdb` | RCSB PDB **1B9O**; UniProt **P00709**; human alpha-lactalbumin; 1.15 A X-ray structure | **Recommended protein input** for human-milk-focused exploratory work |
| `bovine_alpha_lactalbumin_1F6S.pdb` | RCSB PDB **1F6S**, holo-bovine alpha-lactalbumin, 2.2 A X-ray structure | Experimental bovine alpha-lactalbumin starting structure |
| `D_glucose_PubChem_CID_5793_3D.sdf` | PubChem CID **5793**, 3D SDF record | Initial D-glucose conformer |
| `lactose_synthase_glucose_control_1NF5_original.pdb` | RCSB PDB **1NF5**, X-ray lactose synthase--glucose complex, 2.0 A | Immutable experimental positive-control source |

Source URLs:

- https://www.rcsb.org/structure/1B9O
- https://www.rcsb.org/structure/1F6S
- https://www.rcsb.org/structure/1NF5
- https://pubchem.ncbi.nlm.nih.gov/compound/5793

## 1NF5 redocking control

`1NF5` contains two crystallographic copies. For one experimentally observed complex, use alpha-lactalbumin chain **A** and beta1,4-galactosyltransferase chain **B**. Its bound ligand is **beta-D-glucopyranose** (`BGC`, chain `B`, residue `403`; 12 heavy atoms). The matching copy is `BGC D 527`.

Run the reproducible preparation script from this folder:

```powershell
.\prepare_1NF5_redocking_control.ps1
```

It creates the following prepared files:

| File | Content | Intended use |
| --- | --- | --- |
| `lactose_synthase_1NF5_AB_receptor_no_glucose.pdb` | Alpha-lactalbumin chain A + beta4Gal-T1 chain B + structural `CA A 124`; no `BGC B 403`, waters, or crystallization additives | Receptor for glucose redocking |
| `lactose_synthase_1NF5_bound_beta_D_glucose_reference.pdb` | Experimental `BGC B 403` only | Reference pose for post-docking heavy-atom RMSD |

### Validation workflow

1. Prepare an **independent beta-D-glucopyranose** ligand with validated stereochemistry, protonation, and carbohydrate parameters. Do not use the reference-pose coordinates as the docking input.
2. Dock that ligand to `lactose_synthase_1NF5_AB_receptor_no_glucose.pdb`.
3. Run both a blind protocol (whole receptor) and a focused protocol centered on the known beta4Gal-T1 acceptor-sugar site. Label the latter as *focused redocking*, not blind prediction.
4. Compare each predicted pose/cluster with `lactose_synthase_1NF5_bound_beta_D_glucose_reference.pdb` using glucose heavy-atom RMSD, pocket recovery, and key contacts. A docking score by itself is not validation.
5. Transfer an unchanged, documented ligand-preparation and docking protocol to the HMO study only after the control reproducibly recovers the experimental site and pose.

This validates pose recovery for this known lactose-synthase system, not binding affinity, catalytic turnover, or universal performance for other protein--HMO pairs. Classical docking/MD cannot model the lactose-forming bond reaction; that requires a reaction-aware method such as QM/MM.

## Important interpretation

These are **source inputs, not HADDOCK-ready files**. The recommended human structure, PDB 1B9O, contains one alpha-lactalbumin chain (A), one calcium ion, and crystallographic waters. Use chain A after structure preparation; retain the structural calcium only if its coordination and parameterization are handled consistently. PDB 1F6S remains only a bovine comparison structure and contains six crystallographic protein chains (A--F).

Glucose is not an established binary alpha-lactalbumin ligand. Alpha-lactalbumin promotes glucose recognition when it forms the **lactose-synthase complex** with beta-1,4-galactosyltransferase (beta4Gal-T1). Therefore, docking glucose to alpha-lactalbumin alone cannot serve as a positive control that validates the HADDOCK protocol. It can only be a negative/low-specificity stress test.

For a structural positive control, use a known lactose-synthase ternary structure (for example, PDB 1O23) or reproduce a glucose-bound beta4Gal-T1/alpha-lactalbumin complex with the enzyme present. For future protein--HMO work, parameterize carbohydrate stereochemistry, protonation, and topology with a carbohydrate-compatible workflow before docking; do not convert this SDF to PDB and treat it as validated.

## Download integrity

- `human_alpha_lactalbumin_P00709_1B9O.pdb` SHA-256: `EE0B78D7BA67CC1E39BBDD76532BF5D9FF20B50BA1015761013CCFEF4478653C`
- `bovine_alpha_lactalbumin_1F6S.pdb` SHA-256: `6D72773A12B5762B901B48DE3199084919BF0437BE4B249B329D904A4383315B`
- `D_glucose_PubChem_CID_5793_3D.sdf` SHA-256: `19E2D7109E9E7CF1E76899D294EBEB023CB432198204F448754CF266228BB15C`
- `lactose_synthase_glucose_control_1NF5_original.pdb` SHA-256: `5A17F312411B62F155F3C4F26C22CC289251A58E235D82087763F3C7FEAF78F5`
- `lactose_synthase_1NF5_AB_receptor_no_glucose.pdb` SHA-256: `B92E209583B0098BFAC438F13D2B94D7F7F2B1558E86C0E677E6DB537B245036`
- `lactose_synthase_1NF5_bound_beta_D_glucose_reference.pdb` SHA-256: `7688C1C40A2281B5DBA834AAA91EFE91722DBADB6FB357E76CE5CEE4C351B02A`
