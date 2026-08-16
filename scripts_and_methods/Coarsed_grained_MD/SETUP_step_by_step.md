# CG-MD (Method 3C) — Step-by-Step Setup Guide

**Target system: a flat zein/water interface (a local patch of the nanoparticle surface)
with α-lactalbumin adsorbing onto it from solution. One condition: pH 5.5, ~80 mM NaCl,
298 K. No other proteins.**

Companion to [`README.md`](README.md). The README says *what* Method 3C is and *why*.
This file says **what to type, in what order, starting from zero**, on **WSL2 (Ubuntu)**.

> ### Verification status — read this before trusting a number
>
> **✅ VERIFIED — actually executed against your real files** (martinize2 / vermouth 0.15.0,
> pdb2pqr 3.x): everything in **Stages A, B and C** (Sections 3–5). Bead counts, net
> charges, disulfide detection, the `nan` trap, the force-field capitalisation bug, and
> the starting Rg values are measured.
>
> **📐 DESIGN ONLY — arithmetic and literature, not yet run**: everything in
> **Stages V, 1, 2, 3 and 4** (Sections 6–10). Box sizes, chain counts, densities and
> solvent ratios are calculated; the GROMACS commands are written from standard practice
> but have **not** been executed. Expect to debug them. Section 13 lists what to watch.

---

## 0. The design, in one page

### 0.1 Why a surface patch and not a particle

A real 60 nm zein nanoparticle contains **~3,300 chains ≈ 1.7 million protein beads**, and
with water **10–20 million beads total**. That is not a hardware problem you can solve by
moving to HPC — the sampling you need (µs of adsorption) would take months of cluster
time, and you would still have to *assume* how 3,300 chains pack, which nobody knows.

You don't need the particle. α-lactalbumin is ~4 nm across; it only ever "sees" the 8–10 nm
of surface it touches. Over that window a 60 nm sphere is indistinguishable from a plane:

| Observation window | Deviation of a 60 nm sphere from flat |
|---|---|
| 8 nm | 0.27 nm |
| **10 nm** | **0.42 nm** |
| 15 nm | 0.95 nm |

One Martini bead is **~0.47 nm** in diameter. **At the scale the protein actually samples,
the curvature is smaller than a single bead.** A flat slab is therefore not a compromise
forced by compute — it is physically equivalent at this length scale, and that claim has a
number behind it.

### 0.2 Why zein must be a condensed phase, not a dissolved monomer

α-zein is a prolamin: **insoluble in water**, soluble in ~60–95% aqueous ethanol. A single
zein chain freely dissolved in water is not a state that exists. Any design that starts one
zein chain and one milk protein in water and lets them meet is simulating a non-existent
system.

In water, zein is a **dense amorphous phase**. That is what the slab represents.

We model the **product**, not the **process** — the guide does not simulate antisolvent
precipitation (millisecond timescales, mid-run solvent exchange). This is the same logic as
studying a crystal surface without simulating crystallisation. Stage V tests whether the
force field reproduces the solubility behaviour that justifies this.

### 0.3 The geometry

Cross-section in z; x and y are periodic, so the surface is effectively infinite:

```
   ┌─────────────────────────────────┐  z = 18 nm
   │                                 │
   │   water + 80 mM NaCl            │
   │       ● α-LA      ● α-LA        │  <- milk protein, >= 3 nm off the surface
   │                                 │
   ├─────────────────────────────────┤  <- interface = "the NP surface"
   │#################################│
   │###  dense zein, ~25 chains  ####│     6 nm thick, 1.25 g/cm3
   │#################################│
   ├─────────────────────────────────┤  <- second interface (a free replicate)
   │   water + NaCl                  │
   └─────────────────────────────────┘  z = 0
     |<-------- 12 nm, periodic ------->|
```

Periodic x/y means this 12 × 12 nm patch has **no edges** — which removes the edge
artifacts a finite cluster would suffer. You also get **two independent interfaces** per
box, top and bottom, for free.

### 0.4 The five stages

| Stage | System | Beads | Purpose |
|---|---|---|---|
| **V1/V2** | 1 zein chain in water / in 80% ethanol | ~22,000 each | **Force-field validation.** Does Martini 3 reproduce zein's water-insolubility? |
| **1** | ~25 zein chains, no solvent | ~13,000 | Build the dense amorphous phase |
| **2** | slab + water + ions | ~27,500 | Create and equilibrate the interface |
| **3** | + 2–4 α-lactalbumin | ~29,000 | **The experiment.** Adsorption |
| **4** | — | — | Analysis |

**Do V1/V2 first.** They are the cheapest runs and they can invalidate everything downstream.

---

## 1. Why no docking pose is used

Your [`../Molecular_docking/DOCKING_STATUS_AND_NEXT_STEPS.md`](../Molecular_docking/DOCKING_STATUS_AND_NEXT_STEPS.md)
Section 6 closes the docking leg:

- **MEGADOCK poses are geometric artifacts** — the α-lactalbumin top pose has a
  **3-residue interface** (normal is 20–30), the receptor was `zein_blocked.pdb` so the
  ligand could not land elsewhere, and the convergence on residues 37–39 tracks protrusion,
  not chemistry.
- **HDOCK `D2_HDOCK.pdb`** is the least bad pose available (−315.32, confidence 0.96) but
  sits on **one** zein conformer out of five differing by **13.7–21.5 Å** Cα RMSD.
- Importing any pose and applying a structural bias would hold the interface together *by
  construction* — you would measure your own assumption.

There is also a subtler reason the docking poses cannot transfer here: **docking used the
surface of an isolated chain.** The slab's surface is what ~25 chains expose *after*
rearranging at a water interface — hydrophobic residues retract, polar residues turn
outward. These are very likely different surfaces. Producing the second one is itself a
result (Stage 2).

---

## 2. Stage A — Set up WSL2 ✅

> ### Which terminal am I in?
>
> **Only step A1 runs in Windows PowerShell. Everything else in this entire guide runs in
> the Ubuntu terminal.**
>
> Tell them apart by the prompt:
>
> ```
> PS C:\Users\nnjj1>           <- PowerShell (Windows)
> nnjj1@DESKTOP-ABC123:~$      <- Ubuntu (WSL2)   <-- you want this one
> ```
>
> Open Ubuntu from the Start menu, or type `wsl` in PowerShell.
>
> **What `/mnt/c/` means.** WSL2 automatically mounts your Windows C: drive at `/mnt/c/`,
> so the same folder has two addresses pointing at *the same files*:
>
> ```
> Windows:  C:\Users\nnjj1\UMD-work\dairy_protein_USDA\        (backslashes)
> Ubuntu:   /mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/    (forward slashes)
> ```
>
> Edit in one, the change appears in the other immediately.
>
> **But use `/mnt/c/` only to fetch inputs and deposit results — never to run simulations.**
> WSL2 cross-filesystem I/O is slow enough to matter over hundreds of ns. Copy structures
> to `~/cgmd/` (the Linux-native home directory, i.e. `/home/<your-username>/cgmd/`), run
> there, copy results back at the end.

### A1. Install Ubuntu (Windows PowerShell, as Administrator)

```powershell
wsl --install -d Ubuntu
```

Reboot if prompted, open **Ubuntu** from the Start menu, create your Linux user.

### A2. Check you can see the project

```bash
ls /mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/scripts_and_methods
```

> **Do not run simulations on `/mnt/c/`.** WSL2 cross-filesystem I/O is slow enough to
> matter over hundreds of ns. Work in `~/cgmd/`, copy results back at the end.

### A3. System packages

> **Two layers — don't mix them up.** A3 and A4 install into the Linux system itself; the
> Python virtual environment is not created until A5 and is **not** used here.
>
> | Layer | Command | Installs to | `sudo`? | venv? |
> |---|---|---|---|---|
> | System | `sudo apt install` | the whole Linux system | yes | **no** |
> | Python | `pip install` | the virtual environment | **no** | yes |
>
> `build-essential`, `cmake`, `wget`, `git`, `unzip` are compilers and CLI tools, not Python
> packages. `python3-pip` and `python3-venv` are the tools *used to create* a virtual
> environment, so they must exist system-wide first — you cannot install them inside a venv
> that doesn't exist yet. GROMACS (A4) is likewise a compiled C++ program, unrelated to
> Python.
>
> **If you ever type `sudo pip install`, something has gone wrong.**

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential cmake wget git unzip python3-pip python3-venv
```

### A4. GROMACS

```bash
sudo apt install -y gromacs
gmx --version
nvidia-smi          # if this prints a GPU table, a CUDA build is worth doing later
```

CPU-only, installs in two minutes, sufficient for every stage in this guide.

### A5. Python tools — create the virtual environment here

```bash
python3 -m venv ~/cgmd-env          # create it
source ~/cgmd-env/bin/activate      # activate it
pip install --upgrade pip
pip install vermouth mdtraj MDAnalysis pdb2pqr numpy scipy matplotlib
martinize2 --version     # expect: martinize with vermouth 0.15.x
pdb2pqr30 --version
```

On Ubuntu 23.04+ the venv is **not optional** — the system Python is protected, and
`pip install` outside a venv fails with `error: externally-managed-environment`.

You can tell it is active from the prompt:

```
nnjj1@DESKTOP:~$                 <- not active
(cgmd-env) nnjj1@DESKTOP:~$      <- active
```

**Add `source ~/cgmd-env/bin/activate` to the end of `~/.bashrc`.** The venv deactivates
when you close the terminal, and forgetting to reactivate is the most common cause of
`martinize2: command not found`. Note `gmx` is unaffected — it is a system binary and works
either way; only `martinize2` and `pdb2pqr30` need the venv.

> **Do not install DSSP.** martinize2 only accepts DSSP ≤ 3.1.4; Ubuntu ships 4.x, which
> fails. With `mdtraj` installed, bare `-dssp` works — verified on your files.

### A6. Martini 3 force field

```bash
mkdir -p ~/cgmd/ff && cd ~/cgmd/ff
wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip
unzip martini_v300.zip
find . -name "martini_v3.0.0*"
wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/example_input_files/mdps/martini_v3.0_prod.mdp
```

You need `martini_v3.0.0.itp`, `..._ions_v1.itp`, `..._solvents_v1.itp`. The solvents file
also contains the **ethanol** topology needed for Stage V2.

> **Martini3-IDP needs no extra download** — merged into vermouth in July 2025. It changes
> *bonded* parameters only, so the standard `martini_v3.0.0.itp` still supplies the
> interaction matrix.

### A7. insane ✅

```bash
pip install insane
pip install "setuptools<81"        # required — see below
insane -h | head -3
```

> **The version pin is not optional, and `pip install setuptools` alone will NOT fix it.**
>
> insane 1.2.0 (2017) imports `pkg_resources`, a submodule of `setuptools`. Two things have
> since changed: modern Python/conda environments no longer ship setuptools by default, and
> **setuptools 81 deprecated `pkg_resources` — setuptools 84 removed it entirely.** Verified:
>
> | setuptools | `import pkg_resources` | insane |
> |---|---|---|
> | 84.0.0 (what plain `pip install setuptools` gives you) | ✗ `ModuleNotFoundError` | ✗ |
> | 80.10.2 (`pip install "setuptools<81"`) | ✓ with a deprecation warning | ✓ |
>
> So the usual advice — "just install setuptools" — reproduces the identical error. There is
> also no standalone `pkg_resources` package on PyPI; pinning is the only route.
>
> After the fix, insane prints `UserWarning: pkg_resources is deprecated as an API`. That is
> a **warning, not an error** — it still runs.
>
> Downgrading setuptools is low-risk (it is only a build tool and does not affect GROMACS or
> martinize2). If you would rather not, insane is optional — `gmx solvate` + `gmx genion` is
> the alternative given in Stage 2, at the cost of needing a CG water box and radius tuning.

**Stage A checkpoint:** `gmx --version`, `martinize2 --version`, `pdb2pqr30 --version`,
`python -c "import mdtraj"` all succeed.

---

## 3. Stage B — Prepare the two structures ✅

```bash
mkdir -p ~/cgmd/slab && cd ~/cgmd/slab
P=/mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/scripts_and_methods/Molecular_docking/protein_structures
cp $P/zein_model.pdb .
cp $P/whey/alpha_whey_clean.pdb .
```

### 3.1 The α-lactalbumin trap

`alpha_whey_clean.pdb` looks clean — 122 residues, no gaps, all 8 cysteines — but
**Glu121 and Lys122 are missing their side-chain tips** (only N, CA, C, O, CB). Normal for
a crystal structure: the C-terminal tail is disordered. **No chain in 1F6S has it complete**
(all six checked).

martinize2 does **not** error. It writes a bead with **`nan nan nan` coordinates** for
LYS122 SC2, which fails later with a cryptic message. Verified.

### 3.2 Fix missing atoms and set pH in one command

This step does two jobs:

1. **Rebuilds the missing heavy atoms** (§3.1) — without it, martinize2 emits `nan`.
2. **Assigns protonation states, i.e. fixes every charge.**

Job 2 is the more important one, because **Martini has no dynamic protonation**. Atomistic
MD can do constant-pH simulation, where residues protonate and deprotonate during the run.
CG cannot: each bead's charge is frozen at the moment the topology is written and never
changes again.

**So this command is the only place pH exists in the entire workflow.** Skip it and you
silently get default (≈ pH 7) charges no matter what pH you claim to be simulating.

PROPKA computes **structure-specific pKa shifts** — a buried Asp can have a pKa several
units from an exposed one — then pdb2pqr sets the protonation accordingly, and martinize2
maps the result onto charged or neutral beads.

```bash
pdb2pqr30 --ff=AMBER --keep-chain --with-ph 5.5 --titration-state-method propka \
    --pdb-output ala_ph55.pdb  alpha_whey_clean.pdb ala.pqr

pdb2pqr30 --ff=AMBER --keep-chain --with-ph 5.5 --titration-state-method propka \
    --pdb-output zein_ph55.pdb zein_model.pdb        zein.pqr
```

`WARNING: Missing atom CG in residue GLU A 121` is pdb2pqr **reporting what it rebuilds**,
not failing. Changing pH later is a one-line edit here.

**Checkpoint:** all residues complete; residue names still standard (`HIS`, not `HIE`/`HID`
— martinize2 requires the standard names, and pdb2pqr's PDB output preserves them).

### 3.3 Why pH 5.5, and what it actually changes ✅

**Why 5.5 rather than the memo's 6.6:**

1. **It is a measured experimental condition.** Both DLS datasets in
   `../DLVO/data_DLVO/` use a pH grid of **3, 4, 5.5, 7** — 6.6 appears nowhere, so it
   could not be cross-checked against the wet lab.
2. **6.6 sits on zein's own isoelectric point** (~6.2 in the literature). The CG model is
   *bare* zein, so zein's own pI is the relevant one — not the ~4.5–5 pI of the
   zein–caseinate particle seen in the zeta data, which reflects the caseinate coating.
3. **pH 5.5 is colloidally stable in your own data** (radius 68–110 nm, |ζ| 11–23 mV),
   whereas pH 4 flocculates — radius reaches **8404 nm** at 100 mM ionic strength.

**Net charges, measured from the finished CG topologies at both pH values:**

| | pH 5.5 | pH 6.6 |
|---|---|---|
| zein (521 beads) | **+2** | **+2** |
| α-lactalbumin (293 beads) | **−1** | **−4** |

**zein's charge does not change at all.** Its 234 residues contain only **six ionizable
groups** (1 Glu, 1 Lys, 2 Arg, 2 His) — the classic prolamin composition, Gln- and Pro-rich
and almost devoid of charged residues, which is precisely why zein is water-insoluble. Its
titration curve is nearly flat: net charge moves only from +4.9 to +2.0 across pH 4 → 7.

> **The key consequence: in this system pH acts almost entirely on the milk protein, not on
> zein.** Moving from 6.6 to 5.5 leaves zein untouched and cuts α-lactalbumin's charge
> four-fold.

**Two things this implies for interpretation:**

- **Electrostatic attraction is much weaker at pH 5.5** (+2/−1 vs +2/−4). Adsorption
  observed here is therefore largely **hydrophobically driven** — arguably the more relevant
  question for a corona on a hydrophobic surface, and harder to dismiss as "just
  electrostatics". But be clear that it is a different physical question from pH 6.6.
- **Watch α-lactalbumin's own stability.** Its pI is 4.52, so at pH 5.5 it is only one unit
  above its own isoelectric point with a net charge of −1. In reality α-LA is less soluble
  and prone to self-association there. With 2–4 copies in the box they may **aggregate with
  each other instead of adsorbing to the surface** — so track α-LA↔α-LA contacts separately
  from α-LA↔surface contacts in Stage 4, or you will misread the result.

> **Note on PROPKA vs hand calculation.** A naive Henderson–Hasselbalch estimate gives
> α-LA ≈ −4 at pH 5.5; PROPKA gives **−1**, because several carboxylates in the
> calcium-binding site have substantially elevated pKa values. This gap is exactly why the
> `--titration-state-method propka` step is worth running rather than renaming residues by
> hand.

---

## 4. Stage C — Coarse-grain, with different treatments ✅

The two proteins get deliberately different structural treatments because their structures
are of completely different quality:

| | Source | Quality | Treatment |
|---|---|---|---|
| **α-lactalbumin** | 1F6S crystal | experimental, 4 disulfides, no gaps | Martini 3 + **elastic network** |
| **zein** | ColabFold | mean pLDDT 49, 67% < 50, 3 buried residues of 234 | **Martini3-IDP, no network** |

Restrain what you know; let the rest sample. This asymmetry is easy to defend in a methods
section.

**Run martinize2 once per protein, on a single-chain PDB** — never on a multi-chain file,
or the elastic network can bridge between molecules.

```bash
martinize2 -f zein_ph66.pdb -x zein_cg.pdb -o zein.top \
    -ff martini3IDP -p backbone -dssp -name ZEIN

martinize2 -f ala_ph66.pdb -x ala_cg.pdb -o ala.top \
    -ff martini3001 -p backbone -dssp -elastic -el 0 -eu 0.85 -name ALAC
```

> **Capitalisation: `martini3IDP`, not `martini3idp`.** The Martini3-IDP GitHub README says
> lowercase; vermouth 0.15.0 registers it as `martini3IDP` and exits with
> `ValueError: Unknown force field` otherwise. Check yours with `martinize2 -list-ff`.

**Verified output:**

| Molecule | Force field | Residues | Beads | Net charge | Elastic bonds | Rg (start) |
|---|---|---|---|---|---|---|
| `ZEIN_0.itp` | martini3IDP | 234 | 521 | +2 | **0** ✔ | 3.12 nm |
| `ALAC_0.itp` | martini3001 + ENM | 122 | 293 | −4 | present ✔ | 1.39 nm |

**Checks that must pass every time:**

```bash
grep -c "Rubber band" ZEIN_0.itp     # must be 0 — zein must NOT be restrained
grep -c "Rubber band" ALAC_0.itp     # must be 1 — alpha-LA must be
grep -c "nan" zein_cg.pdb ala_cg.pdb # must both be 0
grep -c "moleculetype" ZEIN_0.itp ALAC_0.itp   # must both be 1
```

The α-lactalbumin run must log **four** disulfides (Cys6–120, 28–111, 61–77, 73–91). Fewer
means the input structure is damaged.

---

## 5. Shared MD parameters 📐

Start from the downloaded `martini_v3.0_prod.mdp` and change only what is listed. **Do not
hand-write the non-bonded block** — Martini's cutoffs, `rlist`, `epsilon_r = 15` and
`verlet-buffer-tolerance = -1` are load-bearing.

**Pressure coupling depends on the stage — this is easy to get wrong:**

| Stage | System shape | `Pcoupltype` |
|---|---|---|
| V1, V2 | single protein in solvent | **`isotropic`** |
| 1 | dense phase, compressing z only | **`semiisotropic`**, x/y `compressibility = 0` |
| 2, 3 | slab spanning x/y with water above/below | **`semiisotropic`** |

The reference `.mdp` ships with `semiisotropic`, which is correct for Stages 1–3 but wrong
for V1/V2.

**The four run types:**

| File | Key settings |
|---|---|
| `em.mdp` | `integrator = steep`, `nsteps = 10000`, `emtol = 100`, `define = -DPOSRES` |
| `nvt.mdp` | `dt = 0.005`, `nsteps = 200000` (1 ns), `define = -DPOSRES`, `tcoupl = v-rescale`, `tc-grps = Solvent Rest`, `tau_t = 1.0 1.0`, `ref_t = 298 298`, `pcoupl = no`, `gen_vel = yes`, `gen_temp = 298` |
| `npt.mdp` | as NVT but `dt = 0.01`, `nsteps = 2000000` (20 ns), `pcoupl = c-rescale`, `gen_vel = no` |
| `prod.mdp` | `dt = 0.02`, `nsteps` per stage, no `define`, `nstxout-compressed = 25000` |

The `dt` ramp 5 → 10 → 20 fs matters: jumping straight to 20 fs is the most common cause of
"LINCS warning / particle out of box" crashes in a fresh system.

> **Position restraints must be OFF in production for zein.** `-DPOSRES` is fine during EM
> and NVT, but leaving it on in production re-freezes the conformer you deliberately chose
> not to freeze.

---

## 6. Stage V — Force-field validation 📐

**Do this first. It is cheap and it can invalidate everything downstream.**

A chain of an insoluble polymer collapses in a poor solvent and stays expanded in a good
one. For zein, water is the poor solvent and aqueous ethanol the good one. If Martini 3
reproduces that contrast, the slab rests on something.

| Run | Solvent | Expectation |
|---|---|---|
| **V1** | water + 80 mM NaCl | Rg **collapses** well below the 3.12 nm starting value |
| **V2** | 80% ethanol / 20% water | Rg **stays expanded** |

V2 is the positive control. Without it, a model that collapses *everything* would pass V1.

### V1 — zein in water

```bash
gmx editconf -f zein_cg.pdb -o zein_box.gro -box 14 14 14 -c
insane -f zein_box.gro -o V1_solvated.gro -p tmp.top -pbc cubic -box 14,14,14 \
       -sol W -salt 0.08 -charge 0 -d 0
sed -i 's/NA+/NA /g; s/CL-/CL /g' V1_solvated.gro tmp.top
```

insane writes Martini 2 ion names (`NA+`/`CL-`); Martini 3 uses `NA`/`CL`, hence the `sed`.
Then EM → NVT → NPT (isotropic) → 500 ns production.

### V2 — zein in 80% ethanol

Ethanol is a **single bead** in Martini 3; its topology is in
`martini_v3.0.0_solvents_v1.itp`. Look up the molecule name there before writing the `.top`.

**Composition arithmetic for 80% v/v ethanol** (per 100 mL):

| | Volume | Density | Mass | Moles |
|---|---|---|---|---|
| ethanol | 80 mL | 0.789 g/mL | 63.1 g | 1.370 mol |
| water | 20 mL | 1.000 g/mL | 20.0 g | 1.110 mol |

One Martini `W` bead = 4 real waters, so 1.110 mol water → 0.278 mol of `W` beads.

**Bead ratio ≈ 4.9 ethanol beads per W bead.** Build the box with roughly that ratio
(insane cannot mix two solvents directly — build the ethanol box first, then replace ~17%
of the ethanol beads with `W`, or solvate with a pre-mixed box).

> ⚠️ **This ratio is arithmetic, not a tested build.** Verify the final composition with
> `gmx density` or by counting molecules in the `.gro` before you trust it. Also confirm
> whether your source specifies 80% v/v or w/w — the bead ratio differs.

### The gate

Plot Rg vs time for both. **Do not proceed to Stage 1 unless V1 collapses and V2 does not.**
If both collapse, the model is simply over-sticky and the adsorption results will be
qualitative at best. If neither collapses, something is wrong with the setup.

---

## 7. Stage 1 — Build the dense zein phase 📐

**Target:** ~25 chains in 12 × 12 × 6 nm = 1.25 g/cm³ (verified arithmetic: 864 nm³ ×
1.25 g/cm³ ÷ 25.7 kDa ≈ 25.3 chains; 25 × 521 = 13,025 beads).

Packing 13,000 beads directly into 864 nm³ will defeat `insert-molecules`. Build loose,
then compress.

**Option A — loose pack, then anisotropic compression (recommended):**

```bash
# 1. place 25 chains loosely in a tall box, correct x/y, no solvent
gmx insert-molecules -ci zein_cg.pdb -nmol 25 -box 12 12 30 -o loose.gro -radius 0.21 -try 2000

# 2. energy minimise
gmx grompp -f em.mdp -c loose.gro -p system.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em

# 3. compress z only: semiisotropic, x/y compressibility = 0
gmx grompp -f compress.mdp -c em.gro -p system.top -o compress.tpr -maxwarn 1
gmx mdrun -deffnm compress
```

`compress.mdp`: `Pcoupltype = semiisotropic`, `compressibility = 0 3e-4`,
`ref_p = 1.0 1.0`, `dt = 0.01`, tens of ns. The zero x/y compressibility holds the lateral
area fixed while z collapses under the chains' own cohesion.

**`-radius 0.21` is not optional.** The default 0.105 nm is an atomistic radius; Martini
beads are roughly twice that.

**Option B — fallback:** if `insert-molecules` cannot place 25 chains even in the tall box,
accept fewer (18–22) and let NPT find the equilibrium density. A slightly thinner slab is
fine; the density is what matters.

**Checkpoints:**

- Box z stops shrinking; volume flat over the last ~20 ns
- Final density ≈ 1.2–1.3 g/cm³ (`gmx energy`, select Density)
- No chain has drifted off on its own — the phase should be one connected blob

---

## 8. Stage 2 — Create the interface 📐

```bash
# extend z, keeping the slab centred
gmx editconf -f compress.gro -o slab_box.gro -box 12 12 18 -c

# solvate the empty space
insane -f slab_box.gro -o S2_solvated.gro -p tmp.top -pbc rectangular \
       -box 12,12,18 -sol W -salt 0.08 -charge 0 -d 0
sed -i 's/NA+/NA /g; s/CL-/CL /g' S2_solvated.gro tmp.top
```

Then EM → NVT → NPT (semiisotropic) → **200–500 ns equilibration**. The surface needs time
to reorganise: hydrophobic side chains retract, polar ones turn outward.

**Checkpoints — these matter more than usual:**

```bash
gmx density -f npt.xtc -s npt.tpr -n index.ndx -o density_z.xvg -d Z
```

1. **Clean step profile along z:** a zein plateau, a sharp-ish interface, a water plateau.
2. **No water trapped inside the slab.** Water density in the slab interior should be near
   zero. If solvation pushed water into cavities, rebuild with a larger insertion radius or
   solvate with `gmx solvate -scale` tuned for CG.
3. **The interface has stopped moving** — surface composition stable over the last ~100 ns.

### This stage is a result, not just preparation

**Which residues does zein expose at a water interface?** Compute the residue composition
within ~1 nm of the interface and compare it with the whole-chain composition, and with the
34 zein residues HDOCK predicted in
`../Molecular_docking/protein_structures/whey/interface_residues.csv`.

Docking used the surface of an *isolated chain*. This is the surface of a *condensed phase*.
If they differ, that is a concrete, publishable statement about why single-chain docking
could not answer the corona question — and it strengthens the negative result already
written up in `DOCKING_STATUS_AND_NEXT_STEPS.md`.

---

## 9. Stage 3 — α-lactalbumin adsorption 📐

Place 2–4 α-LA copies in the water phase, **≥ 3 nm from both interfaces**, before solvating:

```bash
# translate copies to chosen z positions in the empty region, then concatenate with the slab
gmx editconf -f ala_cg.pdb -o ala_1.gro -translate <x> <y> <z>
# ... repeat for each copy, combine, then solvate as in Stage 2
```

Placing them *before* solvation is cleaner than `insert-molecules -replace W`, which can
drop a copy inside the slab.

**Verify the starting separation before you run.** If a copy starts in contact with the
surface you have accidentally assumed a binding site:

```bash
gmx mindist -f start.gro -s start.gro -n index.ndx -od check.xvg
```

Then EM → NVT → NPT → **production 500 ns – 1 µs, three replicas with different
`gen_seed`.**

> **Why replicas and multiple copies:** one adsorption event in one trajectory proves almost
> nothing. 3 replicas × 4 copies gives ~12 semi-independent encounters. This is the
> difference between "we observed binding" and "binding occurred in N of 12 encounters."

---

## 10. Stage 4 — Analysis 📐

| Question | Tool / metric |
|---|---|
| Does it adsorb? How fast? | contact count vs time — `gmx mindist -on numcont.xvg -d 0.6` |
| Where on zein? | residue-resolved contact frequency over the trajectory |
| How much surface is buried? | `gmx sasa` |
| Stable or transient? | count adsorption/desorption events; residence-time distribution |
| Does α-LA deform? | its Rg and RMSD (an ENM-restrained protein should barely change — if it does, something is wrong) |
| **Do the α-LA copies aggregate with each other?** | α-LA↔α-LA contacts, tracked **separately** from α-LA↔surface. At pH 5.5, α-LA carries only −1 and sits near its own pI (4.52), so self-association is a real competing outcome, not a nuisance |
| Consistent with 3B? | compare contact residues with `interface_residues.csv` (34 zein / 38 α-LA residues from HDOCK D2) |
| Consistent with the wet lab? | adsorbed layer thickness vs the DLS size shift in `../DLVO/data_DLVO/` |

**What each outcome means:**

- **Adsorbs and stays** → supports a stable corona. Report *where*, and whether it matches
  3B (independent corroboration) or not (more interesting, and more likely).
- **Adsorbs and releases repeatedly** → a dynamic, exchanging corona. This is a real result,
  not a failed run, and it is consistent with the non-specific binding your docking work
  already concluded.
- **Never adsorbs** → at pH 5.5 the charges are only +2 and −1, so electrostatic attraction
  is weak by design; a null result here is more plausible than it would be at pH 6.6. Before
  concluding anything, check that the α-LA copies have not simply **aggregated with each
  other** (see §3.3) — at pH 5.5 they sit close to their own pI.

---

## 11. Compute budget 📐

**Measure before committing:**

```bash
gmx mdrun -deffnm prod -nsteps 50000 -v
tail -20 prod.log        # read the "Performance: ... ns/day" line
```

| Stage | Beads | Rough CPU-only estimate |
|---|---|---|
| V1 / V2 | ~22,000 each | 1–2 days each |
| 1 | ~13,000 | hours |
| 2 | ~27,500 | several days |
| 3 | ~29,000 | ~1 week per replica |

The whole programme fits on your current machine. **No HPC required.**

| | |
|---|---|
| Week 1 | Stages A–C; V1 and V2 launched |
| Week 2 | V-gate passed; Stage 1 and 2 built |
| Weeks 3–4 | Stage 2 equilibration; Stage 3 replicas launched |
| Weeks 5–6 | Analysis and write-up |

---

## 12. Limitations — carry all of these into the write-up

- **Amorphous packing is assumed, not derived.** The internal structure of a real zein
  nanoparticle is unknown. The slab is a plausible dense phase, not a measured one.
- **No caseinate.** Your particle is zein–**caseinate**, and caseinate is the surface
  stabiliser. Real α-lactalbumin approaches a **caseinate-coated** surface, not bare zein.
  This is the single largest remaining gap in the model and the obvious next step. It is
  also why the simulated and measured electrostatics cannot be compared directly — see next
  point.
- **The model's surface charge does not match the measured particle.** Two independent
  reasons, both documented rather than speculative:
  1. **No caseinate.** The measured ζ-potential is dominated by the caseinate coating
     (pI ≈ 4.6), not by zein.
  2. **No deamidation.** The P02859 α-zein sequence contains **41 Gln (18% of the chain)**
     and essentially no acidic residues, giving a sequence pI of **8.26** — whereas zein
     preparations are commonly reported near pI 6.2. Gln→Glu deamidation closes that gap
     fast: converting just 10% of the Gln drops the calculated pI to 5.36, and 20% to 4.36.
     Commercial zein is typically partially deamidated.
  **Therefore: do not claim agreement between the simulated charges and the DLS ζ-potential
  data.** State the model as bare, non-deamidated zein.
  *Optional cheap sensitivity test:* build a second zein topology with ~15–20% of Gln
  mutated to Glu in the input sequence (re-run pdb2pqr + martinize2, a few minutes) and
  compare adsorption. If the two differ strongly the process is electrostatically driven; if
  not, it is hydrophobically driven. Either answer is a result, and it converts a known
  weakness into a controlled variable.
- **Flat = local.** Valid for what one protein touches; says nothing about whole-particle
  multivalency or avidity.
- **No formation history.** The real surface is kinetically frozen during ethanol→water
  precipitation and may differ from the equilibrium surface simulated here.
- **Martini 3 over-estimates protein–protein stickiness for flexible proteins**
  ([Nat. Commun. 2024](https://www.nature.com/articles/s41467-024-50647-9)). Report
  adsorption as qualitative/comparative — **never as an affinity or a ΔG**.
- **α-lactalbumin's Ca²⁺ was stripped.** 1F6S is the calcium-bound form; apo α-LA is a
  molten globule in reality. The elastic network holds the fold, so the run is stable, but
  the binding-loop charge is wrong by +2 and **no thermal-stability claim about α-LA is
  admissible**.
- **Zein's conformational ensemble is unvalidated.** Removing the elastic network is the
  honest choice given pLDDT 49, but it means there is no reference structure to check
  against. Stage V is the only validation available.
- **Scope.** Per the memo: "V-gate passed, one slab equilibrated, adsorption trend across
  3 replicas" is a **complete success** for 3C in a 4–6 week window. It is a stretch goal
  reported alongside 3A/3B/4, not a load-bearing deliverable.

---

## 13. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `martinize2: command not found` | venv not active — `source ~/cgmd-env/bin/activate` |
| `ValueError: Unknown force field "martini3idp"` | capitalisation — use `martini3IDP` |
| DSSP errors | don't install `mkdssp` (4.x incompatible); `pip install mdtraj`, use bare `-dssp` |
| `insane`: `ModuleNotFoundError: No module named 'pkg_resources'` | `pip install "setuptools<81"` — plain `pip install setuptools` gives 84.x, which **removed** `pkg_resources` and fails identically (§A7) |
| `nvidia-smi: command not found` in WSL2 | **do not** `apt install nvidia-utils-*`. NVIDIA drivers must be installed on **Windows**; WSL2 exposes the GPU via `/usr/lib/wsl/lib/`. Installing a Linux driver inside WSL2 breaks the passthrough. GPU is optional for this guide. |
| `nan` coordinates in the CG pdb | incomplete side chains — run `pdb2pqr30` first (§3.2) |
| martinize2 writes `molecule_0/1.itp` with odd sizes | ran on a multi-chain or gapped PDB — one protein per run |
| `insert-molecules` places far fewer chains than asked | expected at high density — use Option A (loose then compress), keep `-radius 0.21` |
| Slab drifts or breaks up during NPT | wrong `Pcoupltype`; needs `semiisotropic` with x/y compressibility 0 during Stage 1 |
| Water appears inside the slab | solvation inserted into cavities — check `gmx density`, rebuild with a larger insertion radius |
| Box z keeps shrinking and won't converge | not yet equilibrated, or `ref_p` set on the wrong axis |
| `grompp: number of coordinates does not match topology` | `[ molecules ]` counts disagree with the `.gro` |
| `Atomtype W not found` | force-field `#include` lines missing or in the wrong order |
| `Invalid group Solvent` | `index.ndx` not passed (`-n`) or groups misnamed |
| LINCS warnings / crash in first ps | `dt` too large too early — ramp 5 → 10 → 20 fs |
| Water freezes into a lattice | add ~10% `WF` antifreeze beads |
| Runs are painfully slow | you're on `/mnt/c/` — move to `~/cgmd/` |

---

## 14. What this changes in `README.md`

1. **§3 Step 0** — "start from the HDOCK best-pose complexes" is superseded. No docked pose
   is used; the adsorption site is an output.
2. **§3 Step 0 (patch)** — the README's "4–9 zein copies" patch becomes a **~25-chain dense
   slab with periodic x/y**, justified by the curvature calculation in §0.1 rather than by
   compute limits.
3. **§3 Step 1** — `-el 0.5 -eu 0.9` → `-el 0 -eu 0.85`, and the elastic network applies to
   **α-lactalbumin only**. Zein uses `-ff martini3IDP` with no network.
4. **§2 tools table** — DSSP/`mkdssp` is a trap on modern Ubuntu. Add `pdb2pqr`, which does
   missing-atom repair and pH assignment in one step and replaces the manual
   Henderson–Hasselbalch renaming procedure.
5. **§3 Step 2** — pressure coupling is **semiisotropic** for the slab, not isotropic.
6. **§4 limitations** — add the absence of caseinate as the primary structural gap, and the
   Martini 3 flexible-protein stickiness caveat as the primary force-field gap.
7. **New** — the Stage V solvent-quality validation has no counterpart in the README and
   should be added; it is the only available test of whether the CG model captures zein's
   defining physical property.
8. **Paths** — docking outputs are in `../Molecular_docking/HDOCK_Huazhong/`, not
   `../Molecular_docking/HADDOCK/`.

---

## 15. References

- Martini 3 protein models tutorial — [I.I setup and structure bias models](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html) · [I.III IDRs](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut3.html)
- Martini 3 force field — [particle definitions (`martini_v300.zip`)](https://cgmartini.nl/docs/downloads/force-field-parameters/martini3/particle-definitions.html) · [solvents (ethanol)](https://cgmartini.nl/docs/downloads/force-field-parameters/martini3/solvents.html)
- Reference production `.mdp` — [Martini MD parameters](https://cgmartini.nl/docs/downloads/example-input-files/md-parameters.html)
- martinize2 / vermouth — [GitHub](https://github.com/marrink-lab/vermouth-martinize) · [eLife 2024](https://doi.org/10.7554/eLife.90627.2)
- Martini3-IDP — [Nat. Commun. 2025](https://www.nature.com/articles/s41467-025-58199-2)
- Rescaling protein–protein interactions for flexible proteins — [Nat. Commun. 2024](https://www.nature.com/articles/s41467-024-50647-9)
- PDB2PQR / PROPKA — [docs](https://pdb2pqr.readthedocs.io/)
- GROMACS — [install guide](https://manual.gromacs.org/documentation/current/install-guide/index.html)
