#!/bin/bash -i
#SBATCH --job-name=sli_hmo
#SBATCH --output=logs_hmo/sim_%A_%a.out
#SBATCH --error=logs_hmo/sim_%A_%a.err
#SBATCH --array=1-80
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=48:00:00
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
cd "$PROJECT_DIR"
# Site-specific: location of the glpk shared libraries on the UvA omics cluster.
# Replace with the equivalent path on your system, or remove if glpk is already on the loader path.
export LD_LIBRARY_PATH=/zfs/omics/software/lib:$LD_LIBRARY_PATH

# Tasks 1-40 are the prebiotic arm, 41-80 the synbiotic arm
if [ ${SLURM_ARRAY_TASK_ID} -le 40 ]; then
  IDX=${SLURM_ARRAY_TASK_ID}
  ARM="hmo"
else
  IDX=$((SLURM_ARRAY_TASK_ID - 40))
  ARM="hmo_pro"
fi
ID=$(sed -n "${IDX}p" id_list_preterm.txt)
Rscript simulate_arena_hmo.R $ID $ARM
