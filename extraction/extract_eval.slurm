#!/bin/bash -i
#SBATCH --job-name=sli_extract
#SBATCH --output=logs_extract/ext_%A_%a.out
#SBATCH --error=logs_extract/ext_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=0:30:00
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
cd "$PROJECT_DIR"
# Site-specific: location of the glpk shared libraries on the UvA omics cluster.
# Replace with the equivalent path on your system, or remove if glpk is already on the loader path.
export LD_LIBRARY_PATH=/zfs/omics/software/lib:$LD_LIBRARY_PATH

# The condition is passed on the command line so that one script covers all eight batches:
#   sbatch --array=1-51 extract_eval.slurm no_treat_breast
#   sbatch --array=1-40 extract_eval.slurm pro_breast
COND=$1
if [ -z "$COND" ]; then
  echo "Usage: sbatch --array=<range> extract_eval.slurm <condition>"; exit 1
fi

# The untreated runs cover all 51 subjects, the intervention runs only the 40 preterm infants
case "$COND" in
  no_treat_*) LIST=id_list.txt ;;
  *)          LIST=id_list_preterm.txt ;;
esac

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)
Rscript extract_eval.R $ID $COND
