#!/bin/bash -i
#SBATCH --job-name=sli_supp
#SBATCH --output=logs_supplement/supp_%A_%a.out
#SBATCH --error=logs_supplement/supp_%A_%a.err
#SBATCH --array=1-40
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
cd "$PROJECT_DIR"
# Site-specific: location of the glpk shared libraries on the UvA omics cluster.
# Replace with the equivalent path on your system, or remove if glpk is already on the loader path.
export LD_LIBRARY_PATH=/zfs/omics/software/lib:$LD_LIBRARY_PATH
ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" id_list_preterm.txt)
Rscript predict_supplements.R $ID