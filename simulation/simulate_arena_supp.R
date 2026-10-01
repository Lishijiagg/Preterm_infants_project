.libPaths("/home/sli/personal/sli/R_libs/4.3")
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript simulate_arena_supp.R <ID> <Breast|Formula>")
}
Ind = args[1]
Diet = args[2]
set.seed(sum(utf8ToInt(Ind)) + 1)
library(BacArena)
library(stringr)

setwd("~/personal/sli/New_BacArena_project/")
supp_conc = 0.05      # mM, concentration at which each predicted supplement is supplied
diet_lc = tolower(Diet)

meta = read.table("metadata.tsv", header = TRUE, row.names = 1, sep = '\t')
all_amount_data = read.table('all_amount_data.txt', sep = '\t', check.names = FALSE)
ASV_GEM_map = read.table('ASV_GEM_map_with_genus.txt', sep = '\t',
                         header = TRUE, row.names = 1)
GEMs_data_list = readRDS('all_GEMs_raw.rds')
mets_used = readRDS(paste0(diet_lc, '_supplement_mets.rds'))

# Personalised medium: the troubleshot diet plus this individual's predicted supplements
base_medium = read.table(paste0("new_", diet_lc, "_diet.txt"), sep = '\t')
colnames(base_medium)[1:2] = c("EXCHANGE.REACTION", "MOLECULAR.CONCENTRATION.IN.TOTAL")
add_mets = setdiff(mets_used[[Ind]], base_medium$EXCHANGE.REACTION)
Ind_medium = rbind(
  base_medium[, 1:2],
  data.frame(EXCHANGE.REACTION = add_mets,
             MOLECULAR.CONCENTRATION.IN.TOTAL = rep(supp_conc, length(add_mets)),
             stringsAsFactors = FALSE))
cat(Ind, " ", Diet, ": ", length(mets_used[[Ind]]), " predicted, ",
    length(add_mets), " newly added to the medium\n", sep = '')

# Constrain every model with this individual's medium
GEMs_constraint = function(GEMs_list, diet) {
  cons = list()
  for (i in names(GEMs_list)) {
    diet_ex = as.character(diet[,1])
    diet_val = as.numeric(diet[,2])
    reaction_v = grep(pattern = 'EX_', GEMs_list[[i]]@react_id, value = TRUE, fixed = TRUE)
    diet_pos = which(diet_ex %in% reaction_v, arr.ind = TRUE)
    diet_ex = diet_ex[diet_pos]
    diet_val = -1 * diet_val[diet_pos]
    diet_val[abs(diet_val) > 1000] = -1000
    m = changeBounds(model = GEMs_list[[i]], lb = rep(0, length(reaction_v)),
                     ub = rep(1000, length(reaction_v)), react = reaction_v)
    m = changeBounds(model = m, lb = diet_val,
                     ub = rep(1000, length(diet_val)), react = diet_ex)
    cons[[i]] = m
  }
  cons
}
cons_GEMs_list = GEMs_constraint(GEMs_data_list, Ind_medium)

temp_meta = subset(meta, ID == Ind)
temp_GEMs_amount = all_amount_data[rownames(temp_meta)]
temp_GEMs_amount = temp_GEMs_amount[temp_GEMs_amount != 0, , drop = FALSE]

add_bacs = list()
for (i in names(cons_GEMs_list)) {
  add_bacs[[i]] = Bac(model = cons_GEMs_list[[i]], limit_growth = TRUE, deathrate = 0.1,
                      minweight = runif(1, min = 0.2, max = 0.4),
                      maxweight = runif(1, min = 0.8, max = 1),
                      growtype = "exponential")
}

ind_arena = Arena(n = 100, m = 100)
for (asv in rownames(temp_GEMs_amount)) {
  temp_GEM_input = ASV_GEM_map[ASV_GEM_map$ASV_ID == asv, ]$renamed_gem
  if (length(temp_GEM_input) != 1 || is.null(add_bacs[[temp_GEM_input]])) {
    stop("GEM not found for ", asv)
  }
  ind_arena = addOrg(object = ind_arena,
                     specI = add_bacs[[temp_GEM_input]],
                     amount = temp_GEMs_amount[asv, ])
}

Ind_diet_ex = as.character(Ind_medium[, 1])
Ind_diet_val = as.numeric(Ind_medium[, 2])
for (met in 1:length(Ind_diet_ex)) {
  ind_arena = addSubs(object = ind_arena,
                      smax = ifelse(Ind_diet_val[met] > 1e-3,
                                    Ind_diet_val[met]/20, Ind_diet_val[met]),
                      mediac = Ind_diet_ex[met],
                      unit = 'mM', addAnyway = TRUE)
}

# Run the simulation for 12 time steps; substrate depletion and spatial saturation of the arena
# are both reached by step 10-11.
cat("Running supplement simulation for ", Ind, " (", Diet, ")\n", sep = '')
ind_eval = simEnv(ind_arena, time = 12, diffusion = TRUE)
out_dir = paste0("Ind_eval_results_supp_", diet_lc)
dir.create(out_dir, showWarnings = FALSE)
saveRDS(ind_eval, file.path(out_dir, paste0(Ind, "_eval_result_supp_", diet_lc, ".rds")))
cat("Finished ", Ind, " (", Diet, ")\n", sep = '')