# Simulation

These jobs run on a cluster and are not reproducible on a laptop. Each SLURM array task builds the
arena for one infant under one medium and one intervention arm and writes a BacArena Eval object.

| Job script | R script | Arm |
|---|---|---|
| simulate_array_breast.slurm | simulate_arena_breast.R | Baseline, breast milk |
| simulate_array_formula.slurm | simulate_arena_formula.R | Baseline, formula |
| simulate_array_breast_pro.slurm | simulate_arena_breast_pro.R | Probiotic arm, breast milk |
| simulate_array_formula_pro.slurm | simulate_arena_formula_pro.R | Probiotic arm, formula |
| simulate_array_supp.slurm | simulate_arena_supp.R | Personalised nutrient supplementation |
| simulate_array_hmo.slurm | simulate_arena_hmo.R | HMO arm, formula background |
| simulate_array_syn.slurm | simulate_arena_syn.R | HMO and probiotic arm, formula background |
| predict_supplements.slurm | predict_supplements.R | Reduced-cost prediction of the limiting substrates |

prepare_hmo_formula_diet.R has no job script. It transfers the 13 HMO exchange reactions of the
breast milk medium onto the gap-filled formula medium and was run interactively. Its output is
data/media/new_formula_hmo_diet.txt.

Order of execution: prepare_hmo_formula_diet.R; the two baseline arms followed by
extraction/extract_eval.R on their output; predict_supplements.slurm, which needs the producers
identified from the baseline simulations; then the probiotic, supplement, HMO and synbiotic arms in
any order, each followed by extract_eval.R.

Before submitting, create the log directories named in the SBATCH output lines, because SLURM writes
the log at job start and fails if the directory is missing: logs_breast, logs_formula,
logs_pro_breast, logs_pro_formula, logs_supp_sim, logs_hmo, logs_syn, logs_supplement, logs_extract.

LD_LIBRARY_PATH is set to the location of the glpk shared libraries on the cluster these jobs were
run on. Replace it with the equivalent path on your system, or remove the line if glpk is already on
the loader path. R was available on PATH rather than loaded through a module; R 4.3.2 with the
packages listed in the root README is required.
