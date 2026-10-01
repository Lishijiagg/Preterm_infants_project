# Data

Every file here is read by `GEMs-BacArena_analysis.R` through `file.path(data_dir, ...)`, where
`data_dir` is set to `data` in the path block of Part 0.

## metadata/

| File | Contents |
|---|---|
| `metadata.tsv` | Subject group (Term, S-Preterm, N-Preterm) and the mapping between sequencing and simulation identifiers |
| `all_select_meta.txt` | Clinical variables of the preterm infants: gestational age, birth weight, postnatal age at sampling, days of antibiotic exposure, blood-culture isolate, and the interval between sampling and diagnosis |

## features/

| File | Contents |
|---|---|
| `asv_count.tsv` | ASV count table after contaminant removal, 545 ASVs |
| `asv_rel.tsv` | Relative abundance table |
| `asv_rel_191_final.tsv` | Relative abundances of the 191 modelled ASVs, the direct input to the simulations |
| `asv_268_with_taxonomy.tsv` | SILVA annotation of the 268 retained ASVs |
| `asv_id_mapping.tsv` | ASV hash to identifier mapping |
| `all_amount_data.txt` | Inoculum sizes after scaling to a community total of 500 |

## media/

| File | Contents |
|---|---|
| `breast_diet.tsv` | Breast milk medium, VMH composition plus 13 HMO exchange reactions |
| `formula_diet.tsv` | Formula medium, VMH composition |
| `new_breast_diet.txt` | Breast milk medium after gap filling, 145 exchange reactions |
| `new_formula_diet.txt` | Formula medium after gap filling, 133 exchange reactions |
| `new_formula_hmo_diet.txt` | Formula medium carrying the 13 HMO exchange reactions, the background of the factorial design |

## mapping/

| File | Contents |
|---|---|
| `ASV_GEM_map_with_genus.txt` | Assignment of each ASV to a strain and a model file, with genus-level annotation. The same information is given in Supplementary Table S3 |

## derived/

Objects produced by earlier stages of the pipeline rather than raw inputs. See `derived/README.md`.

## Not in this repository

**BacArena `Eval` objects.** A single object exceeds 100 GB. The analysis script reads the reduced
tables produced by `extraction/extract_eval.R` instead.

**The 191 `.mat` metabolic models.** These belong to AGORA2 (https://www.vmh.life) and to the
HMO-extended reconstructions of Shaaban et al., and are not redistributed here.
`mapping/ASV_GEM_map_with_genus.txt` is sufficient to reassemble exactly the set used in this study.

## extracted/

Per-infant tables reduced from the BacArena simulation output, ten conditions. These are the direct
input to the analysis script. See extracted/README.txt.
