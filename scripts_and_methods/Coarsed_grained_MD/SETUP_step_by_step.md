# CG-MD (Method 3C) — Step-by-Step Start Guide

Companion to [`README.md`](README.md). The README says *what* Method 3C is and *why*.
This file says **what to type, in what order, starting from zero**, targeting
**GROMACS + martinize2 running under WSL2 (Ubuntu) on this Windows machine**.

Everything in Sections 1–4 below was **actually executed and verified** against your
real files (`D1_HDOCK.pdb`, `zein_model.pdb`) with martinize2 / vermouth 0.15.0. The
bead counts, charges, warnings and the two input problems in Section 1 are measured,
not assumed.

---

## 0. Where to start — the short answer

Three things, in this order. Do not skip ahead.

1. **Fix the input structures** (Section 1). Your current `D1_HDOCK.pdb` will produce a
   *silently wrong* CG topology. This is the highest-priority item and it costs ~30 min.
2. **Install the environment in WSL2** (Section 2). ~1–2 hours, mostly waiting.
3. **Run the 1:1 D1 complex end-to-end as a pipeline test** (Sections 3–6). Treat this
   as a plumbing test, *not* a scientific result — the science starts at the multi-copy
   patch in Section 8.

Realistic first milestone: **a 50 ns D1 trajectory that doesn't crash, within one week.**

---

## 1. Two problems in your current inputs (fix these first)

### 1.1 `D1_HDOCK.pdb` chain A has a 15-residue hole

Chain A (zein) in both `D1_HDOCK.pdb` and `D2_HDOCK.pdb` covers residues **1–230 with
residues 28–42 missing**, whereas `protein_structures/zein_model.pdb` is complete
(1–234, no gaps). HDOCK dropped that loop.

martinize2 does not error on this. It treats the break as two chain ends, caps both with
charged termini, and **regroups the fragments** — so you get physically wrong molecules.

### 1.2 The consequence: elastic bonds weld zein to β-lactoglobulin

Running the README's Step-1 command directly on `D1_HDOCK.pdb` produced two molecule
topologies — and neither one is "zein" or "β-lactoglobulin":

| Output | Contents | Beads |
|---|---|---|
| `molecule_0.itp` | zein residues **1–27 fused with all of β-lactoglobulin** | 414 |
| `molecule_1.itp` | zein residues 43–230 | 421 |

Because zein 1–27 and β-lg ended up inside one `[ moleculetype ]`, the elastic network
was applied *across the docking interface*: **12 harmonic bonds at 700 kJ mol⁻¹ nm⁻²
directly tie the zein fragment to β-lactoglobulin.** Those bonds never break. The
simulation would report a perfectly persistent interface — the exact question Method 3C
exists to answer — purely as an artifact of the topology.

**Two rules that follow, and they apply to every system you build:**

> **Rule 1.** Repair the chain break before coarse-graining — use the complete
> `zein_model.pdb` geometry, re-fitted onto the docked pose.
>
> **Rule 2.** Run martinize2 **once per protein chain, on a single-chain PDB**, never on
> the multi-chain complex. Assemble the complex afterwards from the CG coordinates. This
> makes cross-chain elastic bonds structurally impossible.

Section 7.1 gives a one-line check to confirm you never violate Rule 2.

### 1.3 A judgement call to make now, not later

`zein_model.pdb` residues 1–~21 (`MAAKIFCLLMLLGLSASAATA`) are the **signal peptide**,
which is cleaved in mature α-zein. It is strongly hydrophobic, so leaving it in gives you
a sticky artificial patch that will preferentially grab milk protein and skew every
contact metric.

Your docking (3B) was run *with* it, so removing it now means 3B and 3C use different
molecules. Neither choice is free:

- **Keep it** → consistent with 3B, but a known artifact you must declare in the write-up.
- **Strip residues 1–21** → physically correct, but you should re-dock to keep 3B/3C aligned.

Decide before building, write the decision down, and keep it fixed across all conditions.
(The steps below keep it, to stay consistent with your existing docking.)

---

## 2. Stage A — Set up WSL2

### A1. Install Ubuntu (Windows PowerShell, as Administrator)

```powershell
wsl --install -d Ubuntu
```

Reboot if prompted, then open **Ubuntu** from the Start menu and create your Linux
username/password. Everything from here is typed in the Ubuntu terminal.

### A2. Confirm you can see your project folder

Windows drives are mounted under `/mnt/c/`:

```bash
ls /mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/scripts_and_methods
```

You should see `Coarsed_grained_MD  DLVO  Molecular_docking`.

> **Performance warning:** do **not** run simulations on `/mnt/c/`. Cross-filesystem I/O
> in WSL2 is slow enough to matter over hundreds of ns. Work in the Linux filesystem
> (`~/cgmd/`) and copy final outputs back to `/mnt/c/...` at the end.

### A3. System packages

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential cmake wget git python3-pip python3-venv
```

### A4. GROMACS

Start with the packaged build — it is CPU-only but takes two minutes instead of an hour,
and it is enough to get through Section 6:

```bash
sudo apt install -y gromacs
gmx --version
```

Check the printed version is **2023.x or newer** and note whether it says GPU support is
enabled (with the apt package it will not be).

**Do you have an NVIDIA GPU?** Check:

```bash
nvidia-smi
```

If that prints a GPU table, WSL2 can use CUDA and it is worth building a GPU-enabled
GROMACS later (see [the install guide](https://manual.gromacs.org/documentation/current/install-guide/index.html),
current release 2026.3). Don't do it now — get the pipeline working on CPU first, then
optimize. Section 9 has the honest performance numbers.

### A5. Python tools (martinize2 + analysis)

```bash
python3 -m venv ~/cgmd-env
source ~/cgmd-env/bin/activate
pip install --upgrade pip
pip install vermouth mdtraj MDAnalysis numpy scipy matplotlib biopython
martinize2 --version
```

Expect `martinize with vermouth 0.15.x`. **Add `source ~/cgmd-env/bin/activate` to
`~/.bashrc`** or you will spend an afternoon debugging "command not found".

> **Skip installing DSSP.** The published tutorials mention a `mkdssp` binary, but
> martinize2 only accepts DSSP ≤ 3.1.4 and Ubuntu now ships 4.x, which fails. With
> `mdtraj` installed (above), plain `-dssp` with no argument works — verified on your
> files, correct secondary structure and both β-lactoglobulin disulfides detected.

### A6. Martini 3 force field files

```bash
mkdir -p ~/cgmd/ff && cd ~/cgmd/ff
wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip
unzip martini_v300.zip
ls
```

You need these three (paths inside the zip may be nested — `find . -name "martini_v3.0.0*"`):

- `martini_v3.0.0.itp` — particle definitions and interaction matrix
- `martini_v3.0.0_ions_v1.itp`
- `martini_v3.0.0_solvents_v1.itp`

Also grab the reference production parameters:

```bash
wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/example_input_files/mdps/martini_v3.0_prod.mdp
```

### A7. insane (box builder)

```bash
pip install insane || (cd ~/cgmd && git clone https://github.com/Tsjerk/Insane.git && pip install ./Insane)
insane -h | head -3
```

If neither works, `gmx solvate` + `gmx genion` is a fine substitute — Section 5 gives both.

**Stage A checkpoint:** `gmx --version`, `martinize2 --version` and `python -c "import mdtraj"`
all succeed.

---

## 3. Stage B — Repair the D1 structure

Working directory for this stage:

```bash
mkdir -p ~/cgmd/D1 && cd ~/cgmd/D1
SRC=/mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/scripts_and_methods/Molecular_docking
cp $SRC/HADDOCK/D1_HDOCK.pdb .
cp $SRC/protein_structures/zein_model.pdb .
```

The repair: superpose the complete `zein_model.pdb` onto the docked chain A using the
shared Cα atoms, then keep the complete zein coordinates plus the docked chain B. Save
this as `repair.py`:

```python
from Bio.PDB import PDBParser, PDBIO, Superimposer, Select, Structure, Model

p = PDBParser(QUIET=True)
dock = p.get_structure('d', 'D1_HDOCK.pdb')[0]
full = p.get_structure('f', 'zein_model.pdb')[0]

dA, fA = dock['A'], full['A']
shared = [r.id[1] for r in dA if r.id[1] in {x.id[1] for x in fA}]
sup = Superimposer()
sup.set_atoms([dA[i]['CA'] for i in shared], [fA[i]['CA'] for i in shared])
print(f"fitted on {len(shared)} Ca atoms, RMSD {sup.rms:.3f} A")
sup.apply(list(full.get_atoms()))

st = Structure.Structure('c'); m = Model.Model(0); st.add(m)
fA.id = 'A'; m.add(fA)
b = dock['B']; b.detach_parent(); b.id = 'B'; m.add(b)

class Sel(Select):
    def accept_residue(self, r): return r.id[0] == ' '

io = PDBIO(); io.set_structure(st); io.save('D1_repaired.pdb', Sel())
```

```bash
python repair.py
```

**Expect:** RMSD well under 0.5 Å (the docked chain A *is* this model, just truncated —
a large RMSD means something is wrong, stop and investigate). Output `D1_repaired.pdb`:
4874 atoms, chain A = 3641 (zein 1–234, no gaps), chain B = 1233 (β-lg 5–160).

Now split into single-chain files (Rule 2):

```bash
awk '/^ATOM/ && substr($0,22,1)=="A"' D1_repaired.pdb > zein.pdb; echo END >> zein.pdb
awk '/^ATOM/ && substr($0,22,1)=="B"' D1_repaired.pdb > blg.pdb;  echo END >> blg.pdb
```

### pH: do this here, before martinize2

Martini has no dynamic protonation — the charge state is baked in at this step and cannot
be changed later. Residue counts in D1:

| | Asp | Glu | Lys | Arg | His | Cys |
|---|---|---|---|---|---|---|
| zein (chain A) | 0 | 1 | 1 | 2 | 2 | 3 |
| β-lactoglobulin (chain B) | 10 | 16 | 15 | 3 | 1 | 5 |

**Zein is essentially uncharged** — that is real prolamin chemistry, and it means the
pH story in this system is carried almost entirely by β-lactoglobulin.

At the recommended starting condition (**pH 6.6**), standard pKa values (Asp 3.9, Glu 4.3,
His 6.0, Lys 10.5) leave every Asp/Glu deprotonated and every Lys protonated. Only the
**His residues sit near their pKa** and need a decision. So for pH 6.6 you can proceed
with default protonation — verified net charges from the CG topology: **zein +4,
β-lactoglobulin −7**, which matches the expected ≈ −8 for β-lg at neutral pH (residues
1–4 are absent).

Only when you move to a **different pH** (e.g. near β-lg's pI ≈ 5.2) do you need to
rename residues in the PDB to neutral forms first, guided by PROPKA3 or the
[H++ server](http://newbiophysics.cs.vt.edu/H++/index.php) on `D1_repaired.pdb`.

---

## 4. Stage C — Coarse-grain, one chain at a time

```bash
martinize2 -f zein.pdb -x zein_cg.pdb -o zein.top \
    -ff martini3001 -p backbone -dssp \
    -elastic -el 0 -eu 0.85 -name ZEIN

martinize2 -f blg.pdb -x blg_cg.pdb -o blg.top \
    -ff martini3001 -p backbone -dssp \
    -elastic -el 0 -eu 0.85 -name BLG
```

**Verified output on your files:**

| Molecule | Residues | Beads | Net charge |
|---|---|---|---|
| `ZEIN_0.itp` | 234 | 521 | +4 |
| `BLG_0.itp` | 156 | 356 | −7 |

The β-lg run should log two disulfide bridges (Cys62–Cys156, Cys102–Cys115 in this
file's numbering) — β-lactoglobulin has exactly two plus one free thiol, so this is a
good structural sanity check that the input is intact.

**Note on flags vs. the README.** README §3 Step 1 suggests `-ef 700 -el 0.5 -eu 0.9`.
The current Martini 3 protein tutorial recommends `-el 0 -eu 0.85`, and warns against
force constants below the 700 default (lower values make proteins artificially sticky —
which would directly corrupt an aggregation study). `-ef 700` is the default, so it can
be omitted. Use `-el 0 -eu 0.85`.

`-p backbone` writes `[ position_restraints ]` guarded by `#ifdef POSRES` with a default
`POSRES_FC` of 1000 kJ mol⁻¹ nm⁻². You switch them on in the `.mdp` with
`define = -DPOSRES`, and can soften them with `-DPOSRES_FC=500`.

---

## 5. Stage D — Build the solvated, ionized box

First combine the two CG chains into one starting structure. `zein_cg.pdb` and
`blg_cg.pdb` retain the docked relative geometry, so a plain concatenation preserves the
pose:

```bash
grep -h "^ATOM" zein_cg.pdb blg_cg.pdb > complex_cg.pdb
echo END >> complex_cg.pdb
```

Measured extents: zein spans **9.0 × 6.3 × 7.5 nm** (Rg 3.12 nm — an elongated model),
β-lg is 4.2 × 3.7 × 3.8 nm (Rg 1.45 nm). A **12 nm cubic box** gives adequate clearance
for the 1:1 test.

**Option A — insane (more robust CG water/ion placement):**

```bash
insane -f complex_cg.pdb -o solvated.gro -p system.top \
       -pbc cubic -box 12,12,12 -sol W -salt 0.08 -charge 0 -d 0
```

`-salt 0.08` is 80 mM NaCl — your "whole milk" reference condition. insane writes ion
names `NA+`/`CL-` (Martini 2 nomenclature); Martini 3 uses `NA`/`CL`, so strip the signs:

```bash
sed -i 's/NA+/NA /g; s/CL-/CL /g' solvated.gro system.top
```

**Option B — plain GROMACS.** Note `gmx solvate` needs a pre-equilibrated *CG* water box
(`water.gro`), which is **not** in `martini_v300.zip` — download it from the Martini
[solvent systems examples](https://cgmartini.nl/docs/downloads/example-applications/solvent-systems.html).
If that's a hassle, just use insane.

```bash
gmx editconf -f complex_cg.pdb -o boxed.gro -d 1.5 -bt cubic
gmx solvate  -cp boxed.gro -cs ~/cgmd/ff/water.gro -o solvated.gro -p system.top
# then grompp a throwaway tpr against system.top, and:
gmx genion -s ions.tpr -o solvated.gro -p system.top -pname NA -nname CL -neutral -conc 0.08
```

**Then write `system.top` by hand.** This is the step that breaks most first attempts —
the includes must be in this exact order, and the `[ molecules ]` counts must match the
`.gro` file exactly:

```
#include "./ff/martini_v3.0.0.itp"
#include "./ff/martini_v3.0.0_ions_v1.itp"
#include "./ff/martini_v3.0.0_solvents_v1.itp"
#include "./ZEIN_0.itp"
#include "./BLG_0.itp"

[ system ]
zein - beta-lactoglobulin, CG Martini 3

[ molecules ]
ZEIN_0   1
BLG_0    1
W        <count from insane/solvate>
NA       <count>
CL       <count>
```

Expected total system size: **~14,000 beads** (877 protein + ~13,500 water).

> If you run at or below ~300 K, replace ~10% of `W` with antifreeze `WF` beads —
> standard Martini water spuriously freezes in that range. At 298 K this is worth doing.

### Index file for temperature coupling

Martini needs separate Solvent/Rest coupling groups, and GROMACS' automatic groups don't
handle CG bead names well. Build them explicitly:

```bash
gmx make_ndx -f solvated.gro -o index.ndx
```

At the prompt, type these one line at a time. `make_ndx` prints the number of each new
group as it is created — substitute those numbers where shown in angle brackets:

```
r W WF NA CL          # creates the solvent group, note its number -> N
name N Solvent
! N                   # everything else, note its number -> M
name M Rest
a 1-521               # zein beads (protein comes first in the .gro), number -> P
name P Zein
a 522-877             # beta-lactoglobulin beads, number -> Q
name Q Milk
q
```

(Select by bead index, not residue name — CG residues are still named `MET`, `ALA`, …,
so `r ZEIN` matches nothing. 521 and 877 come from the verified bead counts in Stage C;
recount them for any other system.)

The `Zein` / `Milk` groups are what `gmx mindist` and `gmx sasa` need in Section 7.

---

## 6. Stage E — Minimize, equilibrate, produce

Four `.mdp` files. Start from the downloaded `martini_v3.0_prod.mdp` and change only what
is listed below — **do not hand-write the non-bonded block**, Martini's cutoffs, `rlist`,
`epsilon_r = 15` and `verlet-buffer-tolerance = -1` are all load-bearing.

**One correction to the reference file:** it ships with `Pcoupltype = semiisotropic`,
which is for membranes. For a protein in water you must use:

```
Pcoupltype = isotropic
tau_p      = 4.0
ref_p      = 1.0
compressibility = 3e-4
```

| File | Key settings |
|---|---|
| `em.mdp` | `integrator = steep`, `nsteps = 10000`, `emtol = 100`, `define = -DPOSRES` |
| `nvt.mdp` | `integrator = md`, `dt = 0.005`, `nsteps = 200000` (1 ns), `define = -DPOSRES`, `tcoupl = v-rescale`, `tc-grps = Solvent Rest`, `tau_t = 1.0 1.0`, `ref_t = 298 298`, `pcoupl = no`, `gen_vel = yes`, `gen_temp = 298` |
| `npt.mdp` | as NVT but `dt = 0.01`, `nsteps = 2000000` (20 ns), `pcoupl = c-rescale` isotropic, `gen_vel = no`, drop `define` in the last part |
| `prod.mdp` | `dt = 0.02`, `nsteps = 2500000` (50 ns to start), `tcoupl = v-rescale`, `pcoupl = c-rescale` isotropic, `ref_t = 298 298`, no `define`, `nstxout-compressed = 25000` |

Note the ramp in `dt`: 5 → 10 → 20 fs. Jumping straight to 20 fs during equilibration is
the most common cause of "LINCS warning / particle out of box" crashes in a fresh system.

Run them:

```bash
gmx grompp -p system.top -c solvated.gro -r solvated.gro -f em.mdp   -o em.tpr   -n index.ndx
gmx mdrun  -deffnm em -v

gmx grompp -p system.top -c em.gro   -r em.gro   -f nvt.mdp  -o nvt.tpr  -n index.ndx
gmx mdrun  -deffnm nvt -v

gmx grompp -p system.top -c nvt.gro  -r nvt.gro  -f npt.mdp  -o npt.tpr  -n index.ndx
gmx mdrun  -deffnm npt -v

gmx grompp -p system.top -c npt.gro  -f prod.mdp -o prod.tpr -n index.ndx -maxwarn 1
gmx mdrun  -deffnm prod -v
```

**Checkpoints:**

- After EM: `Fmax < 1000`, potential energy large and negative. If EM fails, the topology
  or `[ molecules ]` counts are wrong — not the physics.
- After NPT: box volume and density flat over the last few ns (`gmx energy`).
- Production: it simply must not crash. Watch the first 100 ps closely.

---

## 7. Stage F — Analysis

### 7.1 The topology check you should run every single time

Before trusting any result, confirm no elastic bond crosses between the two proteins.
Because you martinized each chain separately (Rule 2), this is guaranteed by construction —
but check anyway, because the failure is silent:

```bash
grep -c "moleculetype" ZEIN_0.itp BLG_0.itp   # must be 1 each
```

Both must return `1`. If a single `.itp` contains two `[ moleculetype ]` entries, or the
bead count of one molecule exceeds its chain length, you have reproduced the Section 1.2
bug and the results are meaningless.

### 7.2 The three metrics that matter

```bash
# Interface persistence — contact count over time
gmx mindist -f prod.xtc -s prod.tpr -n index.ndx -od mindist.xvg -on numcont.xvg -d 0.6

# Buried interface area
gmx sasa -f prod.xtc -s prod.tpr -n index.ndx -o sasa.xvg

# Radius of gyration — the metric that cross-checks against your DLS data
gmx gyrate -f prod.xtc -s prod.tpr -n index.ndx -o gyrate.xvg
```

For `mindist`/`sasa` you need index groups for zein and β-lg separately — add them in
`gmx make_ndx` when you build `index.ndx`.

**The one comparison that ties 3C back to 3B:** take the interface residues HDOCK
predicted for D1, and ask whether those specific residues stay within 0.6 nm across the
trajectory. A persistent interface at the *docked* residues supports 3B. A persistent
interface at *different* residues is arguably the more interesting result — it says
docking found the right partner but the wrong pose.

Rg from `gyrate.xvg` is what you compare, as a trend and not an absolute number, against
the DLS particle-size trend in `../DLVO/data_DLVO/`.

---

## 8. Stage G — Only after Section 6 completes cleanly

In priority order:

1. **The multi-copy patch.** This is where the actual science is — a 1:1 pair cannot
   address avidity. Size estimates for your molecules:

   | System | Box | Protein beads | Total beads |
   |---|---|---|---|
   | 1 zein + 1 β-lg (pipeline test) | 12 nm | 877 | ~14,000 |
   | 4 zein + 2 β-lg | 16 nm | 2,796 | ~34,000 |
   | 9 zein + 4 β-lg | 22 nm | 6,113 | ~88,000 |

   On CPU-only WSL2, **4 zein + 2 β-lg is the realistic ceiling.** Build it by
   duplicating `ZEIN_0` in `[ molecules ]` and placing copies with `gmx insert-molecules`,
   arranged as a rough plane to mimic a local NP surface.

2. **GōMartini instead of `-elastic`.** Generate a contact map with the
   [Gō Contact Map server](http://pomalab.ippt.pan.pl/GoContactMap), then
   `-go contact_map.out -go-low 0.3 -go-up 1.1 -go-eps 9.414`. This adds two extra files
   (`go_atomtypes.itp`, `go_nbparams.itp`) that must be `#include`d into
   `martini_v3.0.0.itp` — follow
   [Tutorial I.I §3.3](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html)
   exactly, the include surgery is fiddly and easy to get wrong.

3. **Casein.** Requires a docked fragment from 3B first (you don't have one yet). Use
   IDP-aware settings — `-ff martini3IDP` with `-id-regions`/`-idr-tune`, or water-bias
   tuning. See [Tutorial I.III](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut3.html).
   Standard Martini 3 will give you an artificially collapsed casein.

4. **Second condition / PMF.** Only if 1–3 are done and clean.

---

## 9. Time and hardware — the honest numbers

**Measure, don't guess.** Run a 5-minute benchmark before committing to anything:

```bash
gmx mdrun -deffnm prod -nsteps 50000 -v
tail -20 prod.log     # read the "Performance: ... ns/day" line
```

For orientation: a ~14,000-bead Martini system at dt = 20 fs on 8 CPU cores typically
lands in the **low hundreds of ns/day**, so 50–200 ns of the 1:1 test is an overnight
job. The ~34,000-bead patch will be roughly 2–3× slower. Scale from your own measured
number, not from these.

Rough calendar for the WSL2 CPU path:

| | |
|---|---|
| Week 1 | Stages A–D, `em.tpr` builds without error |
| Week 2 | D1 1:1 runs 50 ns cleanly; analysis scripts working |
| Weeks 3–4 | 4+2 patch at one condition, few hundred ns |
| Weeks 5–6 | Analysis, write-up, GōMartini if time permits |

Per the memo, that is a **complete success** for 3C. It is a stretch-goal method
reported alongside 3A/3B/4, not a load-bearing deliverable.

---

## 10. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `martinize2: command not found` | venv not active — `source ~/cgmd-env/bin/activate` |
| DSSP errors | Don't install `mkdssp` (Ubuntu ships 4.x, incompatible). `pip install mdtraj` and use bare `-dssp` |
| martinize2 writes `molecule_0/1.itp` with unexpected sizes | You ran it on a multi-chain or gapped PDB — Section 1, Rules 1 & 2 |
| `grompp: number of coordinates does not match topology` | `[ molecules ]` counts disagree with the `.gro` — recount `W`/`NA`/`CL` |
| `Atomtype W not found` | Force-field `#include` lines missing or in the wrong order in `system.top` |
| `Invalid group Solvent` | `index.ndx` not passed (`-n`) or groups misnamed |
| LINCS warnings / crash in first ps | `dt` too large too early — ramp 5 → 10 → 20 fs; check EM converged |
| Water freezes into a lattice | Add ~10% `WF` antifreeze beads |
| Runs are painfully slow | You're running on `/mnt/c/`. Move to `~/cgmd/` |
| Barostat + gen_vel warning | Expected; `-maxwarn 1` on the production grompp only |

---

## 11. Corrections this file makes to `README.md`

The README is sound on strategy. Four things need updating from what was verified here:

1. **§3 Step 1** — `-el 0.5 -eu 0.9` → **`-el 0 -eu 0.85`** (current tutorial values).
2. **§2 tools table** — DSSP/`mkdssp` is a trap on modern Ubuntu. Use `mdtraj` + bare
   `-dssp`.
3. **§3 Step 1** — must state explicitly: **one martinize2 run per chain**, never on the
   docked complex. This is the difference between a valid and an invalid study.
4. **§3 Step 0** — must state that the HDOCK poses have a 15-residue gap in zein that
   needs repairing before use.

---

## 12. References

- Martini 3 protein models tutorial — [Tutorial I.I: setup and structure bias models](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html) · [I.III: IDRs](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut3.html)
- Martini 3 force field download — [particle definitions (`martini_v300.zip`)](https://cgmartini.nl/docs/downloads/force-field-parameters/martini3/particle-definitions.html)
- Reference production `.mdp` — [Martini MD parameters](https://cgmartini.nl/docs/downloads/example-input-files/md-parameters.html)
- martinize2 / vermouth — [GitHub](https://github.com/marrink-lab/vermouth-martinize) · [docs](https://vermouth-martinize.readthedocs.io/en/latest/) · [eLife 2024](https://doi.org/10.7554/eLife.90627.2)
- GROMACS installation — [current install guide](https://manual.gromacs.org/documentation/current/install-guide/index.html)
- Martini3-IDP — [Nat. Commun. 2025](https://www.nature.com/articles/s41467-025-58199-2) · [parameters repo](https://github.com/Martini-Force-Field-Initiative/Martini3-IDP-parameters)
- Rescaling protein–protein interactions for flexible proteins — [Nat. Commun. 2024](https://www.nature.com/articles/s41467-024-50647-9)
