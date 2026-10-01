# Derived objects

These six files are produced by earlier stages of the pipeline and are not raw data. They are
included so that the analysis script runs without first re-running the constraint step or the
cluster jobs.

| File | How it is produced | Used from |
|---|---|---|
| `breast_all_GEMs_needed.rds` | `GEMs_constraint()` in Part 1, the 191 models constrained to the breast milk medium | Part 1 onwards |
| `formula_all_GEMs_needed.rds` | The same, formula medium | Part 1 onwards |
| `breast_PRO_GEMs_needed.rds` | The seven probiotic strains constrained to the breast milk medium | Parts 6 and 7 |
| `formula_PRO_GEMs_needed.rds` | The seven probiotic strains constrained to the formula medium | Parts 6 and 7 |
| `breast_supplement_mets.rds` | Output of `simulation/predict_supplements.R` on the cluster, breast milk condition | Part 8 |
| `formula_supplement_mets.rds` | The same, formula condition | Part 8 |

The first four can be regenerated from the model files and the media with `GEMs_constraint()`; the
call is the line immediately above each `readRDS` in the script. The last two require the cluster
job to be re-run.
