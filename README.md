# dairy_protein_USDA

Computational (dry-lab) methods for a USDA-funded project on **zein / zein–caseinate
protein nanoparticles and their interaction with dairy proteins** (whey, caseins).

This repository contains **code and methods only**. Manuscripts, internal reports, and
third-party reference material live alongside it locally but are excluded via
`.gitignore`.

## Contents

| Path | What's in it |
|---|---|
| `scripts_and_methods/DLVO/` | Physics-based DLVO / extended-DLVO colloidal stability model. Debye length, Grahame equation, EDL repulsion, van der Waals attraction, energy barriers, and a sigmoidal ζ(pH) fit used to sweep a continuous pH × ionic-strength stability map. |
| `scripts_and_methods/Molecular_docking/` | Protein–protein docking of zein against whey and casein structures — MEGADOCK and LightDock run scripts, predicted structures (AlphaFold/ColabFold outputs), docked complexes, and scoring/analysis. |
| `scripts_and_methods/Coarsed_grained_MD/` | Guidance for the planned MARTINI 3 + GROMACS coarse-grained MD work (Method 3C) — build/run/analyze plan and limitations. Run scripts not yet implemented. |

Each subdirectory has its own `README.md` with equations, assumptions, limitations,
and run instructions. Start with
[`scripts_and_methods/DLVO/README.md`](scripts_and_methods/DLVO/README.md),
[`scripts_and_methods/Molecular_docking/README.md`](scripts_and_methods/Molecular_docking/README.md), and
[`scripts_and_methods/Coarsed_grained_MD/README.md`](scripts_and_methods/Coarsed_grained_MD/README.md).

## Quick start (DLVO)

```bash
cd scripts_and_methods/DLVO
pip install -r requirements.txt

python dlvo_model.py              # sanity-check the physics module
python run_dlvo_analysis.py --demo  # full pipeline on synthetic demo data
```

Outputs (CSVs + plots) are written to `scripts_and_methods/DLVO/outputs/`.

## A note on the data

`demo_data_zein_caseinate.csv` is **synthetic**, not measured — it exists so the
pipeline can be exercised before the wet-lab window. `anchor_data_template.csv` is a
blank template matching the designed DLS grid. Any results generated from the demo
data are illustrative only.

## Status

Work in progress. Structure, parameters, and conclusions are subject to change as
experimental data come in.
