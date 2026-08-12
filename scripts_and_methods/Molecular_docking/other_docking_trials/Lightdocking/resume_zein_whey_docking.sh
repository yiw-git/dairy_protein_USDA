#!/bin/bash
# ============================================================================
# RESUME script for the zein-whey LightDock pipeline.
#
# Use this instead of re-running run_zein_whey_docking.sh from scratch.
# It figures out what's already done and only does what's left:
#   - run_alpha_whey: simulation already finished (gso_100.out present for
#     swarm_0), but the generate/cluster/rank/filter post-processing never
#     ran. This script checks EVERY swarm for a complete gso_100.out (in
#     case some swarms were still mid-simulation when the process died),
#     resumes only the incomplete ones, then does post-processing.
#   - run_beta_whey: never started. This script runs it from scratch.
#
# Also: the generate+cluster loop is now PARALLELIZED across your 8 cores
# (via xargs -P), instead of doing all ~800 swarms one at a time. That step
# alone should be much faster than before.
#
# RUN THIS IN YOUR WSL TERMINAL:
#   cd /mnt/c/Users/nnjj1/UMD-work/dairy_protein_USDA/scripts_and_methods/Molecular_docking
#   chmod +x resume_zein_whey_docking.sh
#   ./resume_zein_whey_docking.sh
#
# TIP: if you're worried about the terminal closing again mid-run, launch it
# with `nohup ./resume_zein_whey_docking.sh > resume.log 2>&1 &` instead,
# then check progress anytime with `tail -f resume.log`. That way it keeps
# running even if you close the terminal window.
# ============================================================================

set -e

source ~/miniconda3/etc/profile.d/conda.sh
conda activate lightdock

WORK="$HOME/lightdock_runs/zein_whey"
cd "$WORK"

CORES=8
GLOWWORMS=200   # default used during setup; 200 lines of data + 1 header = 201

finish_run () {
  local run_dir=$1
  local ligand_file=$2

  echo ""
  echo "============================================================"
  echo ">>> PROCESSING: $run_dir"
  echo "============================================================"

  cd "$run_dir"

  # --- Check every swarm for a COMPLETE gso_100.out ---
  echo ">>> Checking swarm completion..."
  incomplete=()
  for d in swarm_*/; do
    swarm_id=$(basename "$d" | sed 's/swarm_//')
    n=$(wc -l < "$d/gso_100.out" 2>/dev/null || echo 0)
    expected=$((GLOWWORMS + 1))
    if [ "$n" -lt "$expected" ]; then
      incomplete+=("$swarm_id")
    fi
  done

  if [ ${#incomplete[@]} -gt 0 ]; then
    echo ">>> ${#incomplete[@]} swarm(s) incomplete, resuming simulation for those only..."
    lightdock3.py setup.json 100 -s fastdfire -c $CORES -l "${incomplete[@]}"
  else
    echo ">>> All swarms already fully simulated, skipping to post-processing."
  fi

  # --- Generate conformations + cluster, in parallel across $CORES ---
  echo ">>> Generating conformations and clustering (parallel, $CORES at a time)..."
  echo "    This can take a while with ~800 swarms -- be patient."
  ls -d swarm_*/ | sed 's#/##' | xargs -P $CORES -I{} bash -c \
    "cd {} && lgd_generate_conformations.py ../zein_model.pdb ../$ligand_file gso_100.out $GLOWWORMS > /dev/null 2>&1; lgd_cluster_bsas.py gso_100.out > /dev/null 2>&1"

  # --- Rank and filter ---
  echo ">>> Ranking and filtering by restraint satisfaction..."
  s=$(ls -d ./swarm_* | wc -l)
  lgd_rank.py $s 100
  if [ -d "filtered" ]; then
    echo "    'filtered' folder already exists, skipping filter step (delete it manually to redo)."
  else
    lgd_filter_restraints.py --cutoff 5.0 --fnat 0.4 rank_by_scoring.list restraints.list A B
  fi

  echo ">>> DONE: $run_dir"
  echo "    Results: $run_dir/filtered/rank_filtered.list"
  cd ..
}

# --- Finish the alpha-whey run ---
finish_run "run_alpha_whey" "alpha_whey_clean.pdb"

# --- Beta-whey run: set up + simulate from scratch if it doesn't exist yet ---
if [ ! -d "run_beta_whey" ]; then
  echo ""
  echo "============================================================"
  echo ">>> STARTING: zein vs beta_whey (from scratch)"
  echo "============================================================"
  mkdir -p run_beta_whey
  cd run_beta_whey
  cp ../zein_model.pdb ../beta_whey_clean.pdb ../restraints.list .
  lightdock3_setup.py zein_model.pdb beta_whey_clean.pdb --noxt --noh --now -anm -rst restraints.list
  lightdock3.py setup.json 100 -s fastdfire -c $CORES
  cd ..
fi

finish_run "run_beta_whey" "beta_whey_clean.pdb"

echo ""
echo "============================================================"
echo "ALL DONE."
echo "  $WORK/run_alpha_whey/filtered/rank_filtered.list"
echo "  $WORK/run_beta_whey/filtered/rank_filtered.list"
echo "============================================================"
