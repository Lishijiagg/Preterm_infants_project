# ================================================================================================== #
# Sensitivity analysis 3: stochastic stability of the attribution under repeated simulation
#
# Purpose
#   Analysis 1 bounds what each flux could have been within a single LP. It cannot say whether the
#   realised attribution is stable, because BacArena introduces two further sources of variation that
#   sit outside the LP: the random placement of organisms on the grid, and the order in which
#   organisms are solved within a time step, which determines who reaches a shared substrate first.
#   Both are controlled by the random seed, which was held fixed across the two media so that diet
#   was the only difference between paired simulations. Holding it fixed makes the paired comparison
#   valid but says nothing about how much the result would move under a different draw.
#
#   This script re-simulates a subset of subjects under several seeds and reports how much of the
#   attribution survives. It answers the reviewer question that analysis 1 does not reach: not
#   "could the solver have returned something else" but "does the reported pattern reproduce".
#
#   The same machinery serves stage B of the strain-choice analysis. If sensitivity_02 shows that
#   candidate Bifidobacterium strains differ in their HMO repertoire, set assignment_variants below
#   and the script re-simulates under alternative ASV-to-strain maps instead of, or alongside,
#   alternative seeds.
#
# Cost
#   One simulation per subject per medium per variant. With 10 subjects, 2 media and 5 seeds that is
#   100 simulations. Scale n_subjects and seeds to the wall time available; 5 seeds on 10 subjects is
#   enough to report a concordance figure, and more subjects buy more than more seeds.
#
# THIS SCRIPT NEEDS ONE EDIT BEFORE IT RUNS
#   The simulation call in run_one() is a placeholder. Replace it with the simEnv call and the
#   extraction used in extract_eval.R, so that the sensitivity runs go through exactly the same code
#   path as the main results. Anything else and the comparison is not like for like.
# ================================================================================================== #

library(BacArena)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

dir.create('sensitivity/reseed', showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------------------------------- #
# Configuration
# -------------------------------------------------------------------------------------------------- #
seeds = c(123, 2024, 7, 31415, 99)     # 123 is the seed used for the main results, kept first
n_timesteps = 12

# Which subjects. Term infants are the informative set for the Bifidobacterium question because they
# are the only group carrying the genus; the preterm subjects are the informative set for the
# attribution questions in sections 4 and 6.
subjects_term    = meta$ID[meta$group == 'Term']
subjects_preterm = c(head(meta$ID[meta$group == 'S-Preterm'], 5),
                     head(meta$ID[meta$group == 'N-Preterm'], 5))

subject_set = subjects_preterm        # switch to subjects_term for the strain-choice question

# Alternative ASV-to-strain assignments. Leave as NULL to vary the seed only. To run stage B, supply
# a named list of ASV_GEM_map variants; sensitivity_02 writes the candidate strains that go here.
assignment_variants = NULL

# -------------------------------------------------------------------------------------------------- #
# One simulation. Returns the same four tables that extract_eval.R produces, so that every downstream
# function in major_script.R can be pointed at this output unchanged.
# -------------------------------------------------------------------------------------------------- #
run_one = function(subject, diet_label, seed, map_data = ASV_GEM_map) {

  diet_table = if (diet_label == 'Breast') new_breast else new_formula
  bac_object = if (diet_label == 'Breast') breast_con_GEMs_list else formula_con_GEMs_list

  arena = generate_arena(comp_data = scaled_abundance,
                         metadata  = meta[meta$ID == subject, , drop = FALSE],
                         diet      = diet_table,
                         map_data  = map_data,
                         bac_object = bac_object)[[subject]]

  set.seed(seed)

  # ---------------------------------------------------------------------------------------------- #
  # REPLACE THIS BLOCK with the simEnv call and extraction from extract_eval.R.
  # The arguments that matter for comparability are sec_obj, diff_par, and the number of time steps;
  # they must match the main run exactly.
  # ---------------------------------------------------------------------------------------------- #
  sim = BacArena::simEnv(arena, time = n_timesteps, sec_obj = SEC_OBJ_USED_IN_MAIN_RUN)

  extract_eval_tables(sim, subject_id = subject)    # the extractor from extract_eval.R
}

# -------------------------------------------------------------------------------------------------- #
# Driver. Writes one rds per subject x medium x variant so that an interrupted run resumes rather
# than restarts.
# -------------------------------------------------------------------------------------------------- #
run_grid = function(subjects, diets = c('Breast', 'Formula'),
                    seed_set = seeds, variants = assignment_variants) {

  grid = if (is.null(variants)) {
    expand.grid(subject = subjects, diet = diets, seed = seed_set,
                variant = 'baseline_map', stringsAsFactors = FALSE)
  } else {
    expand.grid(subject = subjects, diet = diets, seed = seed_set[1],
                variant = names(variants), stringsAsFactors = FALSE)
  }

  for (i in seq_len(nrow(grid))) {
    g = grid[i, ]
    tag = sprintf('%s_%s_seed%s_%s', g$subject, g$diet, g$seed, g$variant)
    out_file = file.path('sensitivity/reseed', paste0(tag, '.rds'))
    if (file.exists(out_file)) next

    map_use = if (g$variant == 'baseline_map') ASV_GEM_map else variants[[g$variant]]
    res = run_one(g$subject, g$diet, g$seed, map_data = map_use)
    res = lapply(res, function(tb) {
      if (is.null(tb) || !nrow(tb)) return(tb)
      tb$seed = g$seed; tb$variant = g$variant; tb$Diet = g$diet; tb
    })
    saveRDS(res, out_file)
    message(i, '/', nrow(grid), '  ', tag)
  }
}

run_grid(subject_set)

# ================================================================================================== #
# Analysis of the re-simulation output
# ================================================================================================== #

load_reseed = function(dir_path = 'sensitivity/reseed') {
  f = list.files(dir_path, pattern = '\\.rds$', full.names = TRUE)
  parts = lapply(f, readRDS)
  out = lapply(c('scfa_end', 'hmo_curve', 'growth', 'flux'),
               function(k) bind_rows(lapply(parts, `[[`, k)))
  names(out) = c('scfa_end', 'hmo_curve', 'growth', 'flux')
  out
}

rs = load_reseed()

# -------------------------------------------------------------------------------------------------- #
# Concordance 1: does the same set of taxa come out as the leading contributors?
#
# This is the figure that answers the degeneracy question for section 2. The manuscript reports the
# five leading contributors to each product; the question is how often those five are the same five
# under a different draw. Jaccard similarity against the seed used for the main results is the
# natural summary because the claim in the text is about set membership, not about rank.
# -------------------------------------------------------------------------------------------------- #
top_contributors = function(flux_tbl, n_top = 5) {
  flux_tbl %>%
    filter(rea %in% scfas) %>%
    group_by(seed, variant, Diet, ID, rea, species) %>%
    summarise(net = sum(flux), .groups = 'drop') %>%
    filter(net > 0) %>%
    group_by(seed, variant, Diet, ID, rea) %>%
    slice_max(net, n = n_top, with_ties = FALSE) %>%
    summarise(taxa = list(sort(species)), .groups = 'drop')
}

tops = top_contributors(rs$flux)

reference = tops %>% filter(seed == seeds[1]) %>%
  select(variant, Diet, ID, rea, ref_taxa = taxa)

jaccard = function(a, b) length(intersect(a, b)) / length(union(a, b))

top5_concordance = tops %>%
  filter(seed != seeds[1]) %>%
  left_join(reference, by = c('variant', 'Diet', 'ID', 'rea')) %>%
  filter(!map_lgl(ref_taxa, is.null)) %>%
  mutate(jaccard = map2_dbl(taxa, ref_taxa, jaccard),
         sub = factor(ifelse(rea %in% names(scfa_flux_labels),
                             scfa_flux_labels[rea], rea), levels = scfa_levels)) %>%
  group_by(Diet, sub) %>%
  summarise(median_jaccard = median(jaccard),
            q25 = quantile(jaccard, 0.25), q75 = quantile(jaccard, 0.75),
            pct_identical = 100 * mean(jaccard == 1), .groups = 'drop')

print(top5_concordance, n = 20)

# -------------------------------------------------------------------------------------------------- #
# Concordance 2: does the qualitative class assignment survive?
#
# Section 2 does not claim which individual strain leads, it claims which functional group does.
# That is the weaker and more defensible claim, and this is the number that supports it.
# -------------------------------------------------------------------------------------------------- #
class_shares = rs$flux %>%
  filter(rea %in% scfas) %>%
  group_by(seed, variant, Diet, ID, rea, species) %>%
  summarise(net = sum(flux), .groups = 'drop') %>%
  filter(net > 0) %>%
  left_join(ASV_GEM_map %>% select(strain, Genus), by = c('species' = 'strain')) %>%
  mutate(Class = case_when(
    grepl('Escherichia|Klebsiella|Enterobacter|Citrobacter|Serratia|Proteus', Genus) ~ 'Enterobacteriaceae',
    grepl('Bacteroides|Parabacteroides|Prevotella|Alistipes', Genus)                 ~ 'Bacteroidetes',
    grepl('Bifidobacterium', Genus)                                                  ~ 'Bifidobacterium',
    TRUE                                                                              ~ 'Other')) %>%
  group_by(seed, variant, Diet, ID, rea, Class) %>%
  summarise(net = sum(net), .groups = 'drop') %>%
  group_by(seed, variant, Diet, ID, rea) %>%
  mutate(share = 100 * net / sum(net)) %>%
  ungroup()

class_stability = class_shares %>%
  group_by(variant, Diet, ID, rea, Class) %>%
  summarise(mean_share = mean(share), sd_share = sd(share), .groups = 'drop') %>%
  group_by(Diet, rea, Class) %>%
  summarise(median_mean_share = median(mean_share),
            median_sd_share = median(sd_share, na.rm = TRUE), .groups = 'drop') %>%
  mutate(sub = ifelse(rea %in% names(scfa_flux_labels), scfa_flux_labels[rea], rea))

print(class_stability, n = 40)

# -------------------------------------------------------------------------------------------------- #
# Concordance 3: how much do the end-point concentrations move?
#
# The manuscript reports medians of these. A coefficient of variation across seeds that is small
# relative to the between-group differences supports the group comparisons; one that is comparable
# to them does not, and the affected products should then be reported qualitatively.
# -------------------------------------------------------------------------------------------------- #
conc_stability = rs$scfa_end %>%
  filter(rea %in% scfas) %>%
  group_by(variant, Diet, ID, rea) %>%
  summarise(mean_c = mean(conc), sd_c = sd(conc),
            cv = ifelse(mean_c > 0, sd_c / mean_c, NA_real_), .groups = 'drop') %>%
  group_by(Diet, rea) %>%
  summarise(median_cv = median(cv, na.rm = TRUE),
            q75_cv = quantile(cv, 0.75, na.rm = TRUE),
            n_with_zero_mean = sum(mean_c == 0), .groups = 'drop') %>%
  mutate(sub = ifelse(rea %in% names(scfa_flux_labels), scfa_flux_labels[rea], rea))

print(conc_stability, n = 20)

# -------------------------------------------------------------------------------------------------- #
# Concordance 4: HMO consumption. This is the most robust measurement in the study because it rests
# on mass balance rather than on flux partitioning, so the expectation is a very small CV. Confirming
# that is worth doing, since the manuscript leans on these numbers most heavily.
# -------------------------------------------------------------------------------------------------- #
hmo_stability = rs$hmo_curve %>%
  group_by(variant, Diet, ID, rea, seed) %>%
  summarise(consumed = 100 * (first(conc) - last(conc)) / first(conc), .groups = 'drop') %>%
  group_by(variant, Diet, ID, rea) %>%
  summarise(mean_c = mean(consumed), sd_c = sd(consumed), .groups = 'drop') %>%
  group_by(Diet, rea) %>%
  summarise(median_sd = median(sd_c, na.rm = TRUE),
            max_sd = max(sd_c, na.rm = TRUE), .groups = 'drop') %>%
  mutate(sub = names(hmo_reactions)[match(rea, hmo_reactions)])

print(hmo_stability, n = 20)

# -------------------------------------------------------------------------------------------------- #
# Supplementary table for the manuscript
# -------------------------------------------------------------------------------------------------- #
supp_table_S23 = bind_rows(
  top5_concordance %>%
    transmute(Analysis = 'Top-5 contributor set', Medium = Diet, Compound = as.character(sub),
              Statistic = 'Median Jaccard vs main run', Value = round(median_jaccard, 3)),
  conc_stability %>%
    transmute(Analysis = 'End-point concentration', Medium = Diet, Compound = sub,
              Statistic = 'Median CV across seeds', Value = round(median_cv, 3)),
  hmo_stability %>%
    transmute(Analysis = 'HMO consumption', Medium = Diet, Compound = sub,
              Statistic = 'Median SD across seeds (pp)', Value = round(median_sd, 3))
)

write.table(supp_table_S23, file = 'Supplementary_Table_S23.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)
