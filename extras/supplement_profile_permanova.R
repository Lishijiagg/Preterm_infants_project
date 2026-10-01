# ================================================================================================== #
# Nutritional requirement profiles: PERMANOVA between the two preterm groups
#
# What this is for
#   The STATISTICAL ANALYSIS section states that the predicted supplementation profiles were compared
#   between the two preterm groups by PERMANOVA, but no such result appears in the Results. This
#   script produces it.
#
#   The expected outcome is a null result, consistent with the absence of any systematic difference
#   between the two preterm groups elsewhere in the study, and it is worth reporting for that reason
#   rather than in spite of it.
#
# Why the test is on substrates rather than on nutritional categories
#   An earlier version aggregated the substrates into nutritional categories, as the Methods section
#   currently describes. Category boundaries are a judgement call that has to be defended to a
#   reviewer, and a first pass left 113 of 205 substrates unclassified, so the counts would have
#   reflected the size of the residual category more than the composition of the profiles. Testing
#   the substrate composition directly needs no such scheme and uses the full resolution of the data.
#   The Methods sentence changes accordingly; the wording is at the bottom of this file.
#
# Prerequisites
#   breast_supplements, formula_supplements, meta, preterm_ids, preterm_levels
#   Packages: vegan, dplyr, tidyr, tibble
# ================================================================================================== #

library(vegan)
library(dplyr)
library(tidyr)
library(tibble)

dir.create('sensitivity', showWarnings = FALSE)

permutations = 999
seed = 123

# -------------------------------------------------------------------------------------------------- #
# One row per infant per substrate. round == 1 keeps the substrates identified in the first pass of
# the reduced-cost procedure, matching how the profiles were defined for the intervention itself.
# -------------------------------------------------------------------------------------------------- #
supp_long = rbind(breast_supplements  %>% mutate(Diet = 'Breast'),
                  formula_supplements %>% mutate(Diet = 'Formula')) %>%
  filter(round == 1, ID %in% preterm_ids) %>%
  distinct(Diet, ID, group, react) %>%
  as_tibble()

cat('\n--- Profiles available ---\n')
print(supp_long %>% count(Diet, group, name = 'substrate_entries'))
print(supp_long %>% distinct(Diet, ID) %>% count(Diet, name = 'infants'))

# -------------------------------------------------------------------------------------------------- #
# The two media give very different profiles, a median of 5 substrates under breast milk against 82
# under formula, so they are tested separately rather than pooled.
# -------------------------------------------------------------------------------------------------- #
substrate_matrix = function(diet_label) {
  wide = supp_long %>%
    filter(Diet == diet_label) %>%
    mutate(present = 1L) %>%
    pivot_wider(names_from = react, values_from = present, values_fill = 0L)
  list(mat  = wide %>% select(-ID, -group) %>% as.data.frame(),
       meta = wide %>% select(ID, group))
}

run_permanova = function(diet_label) {
  p = substrate_matrix(diet_label)
  m = p$mat
  g = factor(p$meta$group, levels = preterm_levels)

  # An infant with an empty profile gives an all-zero row, for which Bray-Curtis is undefined.
  keep = rowSums(m) > 0
  if (any(!keep)) {
    message(diet_label, ': ', sum(!keep), ' infants with an empty profile were dropped')
  }
  m = m[keep, , drop = FALSE]
  g = g[keep]

  d = vegdist(m, method = 'bray')

  set.seed(seed)
  ad = adonis2(d ~ g, permutations = permutations)
  set.seed(seed)
  bd = permutest(betadisper(d, g), permutations = permutations)

  list(diet = diet_label,
       n = length(g),
       n_by_group = table(g),
       n_substrates = ncol(m),
       # How many infants have a profile shared with no other. A low count bounds what the test can
       # resolve regardless of the p value, and belongs in the Results if it is low.
       n_distinct_profiles = nrow(unique(m)),
       median_size = median(rowSums(m)),
       adonis = ad,
       betadisper = bd)
}

res_breast  = run_permanova('Breast')
res_formula = run_permanova('Formula')

for (r in list(res_breast, res_formula)) {
  cat('\n================ ', r$diet, ' ================\n')
  cat('infants:', r$n,
      ' substrates:', r$n_substrates,
      ' distinct profiles:', r$n_distinct_profiles,
      ' median profile size:', r$median_size, '\n')
  print(r$n_by_group)
  cat('\nPERMANOVA\n');  print(r$adonis)
  cat('\nbetadisper\n'); print(r$betadisper)
}

# -------------------------------------------------------------------------------------------------- #
# The numbers the Results sentence needs
# -------------------------------------------------------------------------------------------------- #
permanova_summary = tibble(
  Medium = c('Breast', 'Formula'),
  n = c(res_breast$n, res_formula$n),
  n_substrates = c(res_breast$n_substrates, res_formula$n_substrates),
  n_distinct_profiles = c(res_breast$n_distinct_profiles, res_formula$n_distinct_profiles),
  R2 = c(res_breast$adonis$R2[1], res_formula$adonis$R2[1]),
  p = c(res_breast$adonis$`Pr(>F)`[1], res_formula$adonis$`Pr(>F)`[1]),
  betadisper_p = c(res_breast$betadisper$tab$`Pr(>F)`[1],
                   res_formula$betadisper$tab$`Pr(>F)`[1]))

cat('\n--- For the Results sentence ---\n')
print(as.data.frame(permanova_summary))

write.table(permanova_summary, file = 'sensitivity/supplement_profile_permanova.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

# ================================================================================================== #
# Methods wording to replace the category-based sentence
#
# Search for:
#   For nutritional supplementation requirements, identified nutrients were aggregated into count
#   profiles according to their nutritional categories, and differences in nutritional requirement
#   profiles between the two preterm infant groups were assessed using PERMANOVA based on
#   Bray-Curtis dissimilarities.
#
# Replace with:
#   Differences in the composition of the predicted supplementation profiles between the two preterm
#   groups were assessed using PERMANOVA on Bray-Curtis dissimilarities of the substrate
#   presence-absence matrix, with 999 permutations and separately for each medium, and homogeneity of
#   within-group dispersion was evaluated using betadisper.
# ================================================================================================== #
