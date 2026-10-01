# Extraction

extract_eval.R reduces each BacArena Eval object to the small per-infant tables the analysis script
reads. This step exists because an Eval object cannot be held in memory as a whole.

Tables produced per infant and condition:

| Table | Contents |
|---|---|
| Abundance | Number of individuals of each model at the final time step |
| Metabolite concentrations | Concentration of each metabolite in the environment at every time step |
| Exchange flux | Flux of every exchange reaction of every individual at every time step |
| _mediac.rds | Name mapping for the medium components |

The exchange flux table is the raw material for the flux attribution in Parts 4, 5, 6 and 10 of the
analysis script. load_extracted() in Part 0 reads the .rds files written here.
