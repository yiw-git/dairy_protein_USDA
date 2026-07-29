#!/bin/bash
# ============================================================================
# LightDock docking pipeline: zein-alpha (receptor) vs whey proteins (ligand)
# Runs TWO docking jobs: zein vs alpha-whey, and zein vs beta-whey.
#
# WHERE TO RUN THIS: inside your WSL Ubuntu terminal (the bash prompt showing
# "(lightdock) nnjj1@YiWangDell:...$"), NOT inside PyMOL's console.
#
# PREREQUISITE: the "lightdock" conda environment must already exist with
# LightDock installed and working (you already confirmed this with
# `lightdock3.py --help`).
#
# HOW TO RUN:
#   1. Save this file (already done if you're reading it from your workspace
#      folder in Windows Explorer / Cowork).
#   2. In WSL: cd to wherever you saved it, e.g.
#      cd /mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/scripts_and_methods/Molecular_docking
#   3. Make it executable and run it:
#      chmod +x run_zein_whey_docking.sh
#      ./run_zein_whey_docking.sh
#
# This will take a while (100 GSO steps per swarm, ~40-55 swarms per run,
# x2 runs). Grab a coffee.
# ============================================================================

set -e  # stop immediately if any step fails, so you don't waste time on bad output

# --- 1. Activate the conda environment -------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate lightdock

# --- 2. Paths ----------------------------------------------------------------
# Your Windows folder is auto-mounted in WSL under /mnt/c
SRC="/mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/scripts_and_methods/Molecular_docking/protein_structures"

# Work in your Linux home directory instead of /mnt/c directly -- much faster
# disk I/O for the thousands of small files LightDock generates.
WORK="$HOME/lightdock_runs/zein_whey"
mkdir -p "$WORK"
cd "$WORK"

echo ">>> Copying structures and restraint list into $WORK"
cp "$SRC/zein_model.pdb" .
cp "$SRC/whey/alpha_whey_clean.pdb" .
cp "$SRC/whey/beta_whey_clean.pdb" .
cp "$SRC/zein_restricted.txt" .

# --- 3. Convert zein_restricted.txt into LightDock's restraint format -------
# zein_restricted.txt looks like:  4,20,22,25,...,231:A
# LightDock needs one line per residue:  R A.<RESNAME>.<RESNUM>
# ("R" = receptor restraint, since zein is the receptor in this run)
echo ">>> Building restraints.list from zein_restricted.txt"
python3 << 'PYEOF'
with open("zein_restricted.txt") as f:
    line = f.read().strip()
resnums_str, chain = line.split(":")
resnums = [int(x) for x in resnums_str.split(",")]

resname_map = {}
with open("zein_model.pdb") as f:
    for l in f:
        if l.startswith("ATOM"):
            resname = l[17:20].strip()
            resnum = int(l[22:26].strip())
            resname_map[resnum] = resname

with open("restraints.list", "w") as out:
    for n in resnums:
        out.write(f"R {chain}.{resname_map[n]}.{n}\n")

print(f"    restraints.list written with {len(resnums)} receptor restraints")
PYEOF

# Number of CPU cores to use (adjust if needed -- 8 matches what PyMOL
# detected on your machine earlier)
CORES=8

run_docking () {
  local run_name=$1
  local ligand_file=$2

  echo ""
  echo "============================================================"
  echo ">>> RUN: zein vs $run_name"
  echo "============================================================"

  mkdir -p "run_${run_name}"
  cd "run_${run_name}"
  cp ../zein_model.pdb ../"$ligand_file" ../restraints.list .

  # --- Setup: generates swarms, restrained to the zein surface residues ---
  lightdock3_setup.py zein_model.pdb "$ligand_file" --noxt --noh --now -anm -rst restraints.list

  # --- Simulation: 100 GSO steps across all swarms ---
  lightdock3.py setup.json 100 -s fastdfire -c $CORES

  # --- Generate docked structures + cluster per swarm ---
  s=$(ls -d ./swarm_* | wc -l)
  swarms=$((s-1))
  echo ">>> Generating and clustering structures for $s swarms..."
  for i in $(seq 0 $swarms); do
    (cd swarm_$i && lgd_generate_conformations.py ../zein_model.pdb "../$ligand_file" gso_100.out 200 > /dev/null 2>&1)
    (cd swarm_$i && lgd_cluster_bsas.py gso_100.out > /dev/null 2>&1)
  done

  # --- Rank by score, then filter by restraint satisfaction (>=40%) ---
  lgd_rank.py $s 100
  lgd_filter_restraints.py --cutoff 5.0 --fnat 0.4 rank_by_scoring.list restraints.list A B

  echo ">>> DONE: zein vs $run_name"
  echo "    Top poses ranked in: run_${run_name}/filtered/rank_filtered.list"
  cd ..
}

# --- 4. Run both dockings -----------------------------------------------------
run_docking "alpha_whey" "alpha_whey_clean.pdb"
run_docking "beta_whey"  "beta_whey_clean.pdb"

echo ""
echo "============================================================"
echo "ALL DONE."
echo "Results:"
echo "  $WORK/run_alpha_whey/filtered/rank_filtered.list"
echo "  $WORK/run_beta_whey/filtered/rank_filtered.list"
echo "Open the top-ranked swarm_X/lightdock_*.pdb files in PyMOL to inspect."
echo "============================================================"
