#!/bin/bash -i
#SBATCH --job-name=sli_pro_f
#SBATCH --output=logs_pro_formula/sim_%A_%a.out
#SBATCH --error=logs_pro_formula/sim_%A_%a.err
#SBATCH --array=1-40
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=48:00:00

PROJECT_DIR="${PROJECT_DIR:-$PWD}"
cd "$PROJECT_DIR"

# Site-specific: location of the glpk shared libraries on the UvA omics cluster.
# Replace with the equivalent path on your system, or remove if glpk is already on the loader path.
export LD_LIBRARY_PATH=/zfs/omics/software/lib:$LD_LIBRARY_PATH

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" id_list_preterm.txt)

Rscript simulate_arena_formula_pro.R $ID

