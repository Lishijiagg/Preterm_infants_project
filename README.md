# Personalized in silico gut microbiome metabolic modelling of preterm infants

Code and data accompanying the manuscript by Li et al.

Personalized genome-scale metabolic communities were built for 51 infants from 16S rRNA gene
sequencing data, simulated under a breast milk and a formula medium with BacArena, and used to test
three interventions in the 40 preterm communities.

## Layout

    GEMs-BacArena_analysis.R    all analyses and figures, in the order of the manuscript
    simulation/                 cluster job scripts and the R scripts they call
    extraction/                 reduction of BacArena output to the per-infant tables
    data/                       input tables and extracted simulation output, see data/README.md
    extras/                     analyses explored but not reported
    sessionInfo.txt             package versions

## Pipeline

The simulations are not reproduced by the analysis script. The pipeline has three stages.

**Simulation.** simulation/ holds the SLURM array scripts and the R scripts they call. Each task
builds the arena for one infant under one medium and one intervention arm and writes a BacArena Eval
object. A single Eval object exceeds 100 GB, so these are not distributed.

**Extraction.** extraction/extract_eval.R reduces each Eval object to the per-infant tables the
analysis uses: end-point abundances, metabolite concentrations, and the exchange flux of every
organism at every time step. These tables are in data/extracted/ and are the direct input to the
analysis.

**Analysis.** GEMs-BacArena_analysis.R runs top to bottom from the repository root and produces every
number, table and figure in the manuscript. Set the four model directories in the path block of
Part 0 first. Parts 0 and 1 must be run before the rest; after that each part depends only on the
parts above it.

## Script structure

| Part | Contents | Manuscript |
|---|---|---|
| 0 | Environment, parameters, helper functions | setup |
| 1 | Input data, GEM library, media, arena construction | Methods |
| 2 | Fidelity of the simulated communities | Results 1 |
| 3 | Blood-culture isolates and gut abundance | Results 1 |
| 4 | Baseline fermentation product production | Results 2 |
| 5 | Nutritional conditions and HMO utilisation | Results 3 |
| 6 | Probiotic intervention | Results 4 |
| 7 | Determinacy of the flux attribution | Results 2, closing paragraph |
| 8 | Personalised nutrient supplementation | Results 5 |
| 9 | Combined probiotic and nutrient intervention | Results 5 |
| 10 | Prebiotic and synbiotic intervention on formula | Results 6 |
| 11 | Robustness of the lactate results | Results 6 |
| 12 | Supplement composition and cross-arm comparison | Results 5, supporting |

Part 7 is the one departure from manuscript order. The sensitivity analyses report on both the
resident models and the seven introduced strains, so they need the probiotic models already
constrained to each medium by Part 6.

## Metabolic models

The 191 genome-scale metabolic models and the seven probiotic models are too large for this
repository and are archived at Zenodo: https://doi.org/10.5281/zenodo.23090424

They originate from AGORA2 (Heinken et al., Nat Biotechnol 2023, https://www.vmh.life) and from the
HMO-extended reconstructions of Shaaban et al. (Commun Med 2024), the latter covering 119 of the 191
strains. data/mapping/ASV_GEM_map_with_genus.txt gives the strain and model file assigned to every
ASV; the same information is in Supplementary Table S3.

## Requirements

R 4.3.2 with BacArena 1.8.2, sybil, glpkAPI, R.matlab, tidyverse, vegan, rstatix and patchwork.
Exact versions are in sessionInfo.txt.

## Licence

Code is released under the MIT Licence. The metabolic models are distributed under the terms of
their original sources and are not covered by this licence.
