.libPaths("/home/sli/personal/sli/R_libs/4.3")
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript predict_supplements.R <ID>")
}
Ind = args[1]
library(BacArena)
library(stringr)

setwd("~/personal/sli/New_BacArena_project/")
scfas = c("EX_ac(e)","EX_ppa(e)","EX_succ(e)","EX_lac_D(e)","EX_lac_L(e)","EX_for(e)","EX_but(e)")
growth_gain_cutoff = 1.5

# The unconstrained models are used here: the function applies the medium itself, so a model that
# has already been through GEMs_constraint would be constrained twice and give a different answer
GEMs_data_list = readRDS('all_GEMs_raw.rds')

beneficial_met_pred = function(model, medium_comp, eval_obj,
                               redcost_cutoff = -1e-6, max_rounds = 3) {
  diet_ex = as.character(medium_comp[,1])
  diet_val = as.numeric(medium_comp[,2])
  reactions_v = grep(pattern = "EX_", model@react_id, value = TRUE, fixed = TRUE)
  diet_pos = which(diet_ex %in% reactions_v, arr.ind = TRUE)
  diet_ex = diet_ex[diet_pos]
  diet_val = -1 * diet_val[diet_pos]
  
  m = changeBounds(model = model, lb = rep(0, length(reactions_v)),
                   ub = rep(1000, length(reactions_v)), react = reactions_v)
  m = changeBounds(model = m, lb = diet_val,
                   ub = rep(1000, length(diet_val)), react = diet_ex)
  
  found = data.frame(react = character(0), round = integer(0), obj = numeric(0),
                     stringsAsFactors = FALSE)
  for (r in seq_len(max_rounds)) {
    o = optimizeProb(m, poCmd = list("getRedCosts"))
    rc = as.numeric(postProc(o)@pa[[1]])
    ex_idx = grep("EX_", m@react_id, fixed = TRUE)
    lim = m@react_id[ex_idx][rc[ex_idx] < redcost_cutoff]
    lim = setdiff(lim, found$react)
    if (length(lim) == 0) break
    found = rbind(found, data.frame(react = lim, round = r, obj = o@lp_obj,
                                    stringsAsFactors = FALSE))
    m = changeBounds(m, lb = rep(-1000, length(lim)),
                     ub = rep(1000, length(lim)), react = lim)
  }
  
  found[!found$react %in% names(getVarSubs(object = eval_obj, show_products = TRUE)), ]
}

run_one = function(diet_type) {
  eval_file = file.path(paste0("Ind_eval_results_no_treat_", tolower(diet_type)),
                        paste0(Ind, "_eval_result_", tolower(diet_type), ".rds"))
  eval_obj = readRDS(eval_file)
  medium_comp = read.table(paste0("new_", tolower(diet_type), "_diet.txt"), sep = '\t')
  
  mod_id_to_key = setNames(names(GEMs_data_list),
                           vapply(GEMs_data_list, function(m) m@mod_id, character(1)))
  scfa_flux = plotReaActivity(simlist = eval_obj, reactions = scfas, ret_data = TRUE)
  producer_ids = unique(subset(scfa_flux, mflux > 0)$spec)
  producer_keys = mod_id_to_key[producer_ids]
  producer_keys = producer_keys[!is.na(producer_keys)]
  
  cat(Ind, " ", diet_type, ": ", length(producer_keys), " SCFA producers\n", sep = '')
  
  res = lapply(setNames(producer_keys, producer_keys), function(k) {
    out = beneficial_met_pred(model = GEMs_data_list[[k]], medium_comp = medium_comp,
                              eval_obj = eval_obj)
    if (nrow(out) > 0) out$model = k
    out
  })
  
  res_df = do.call(rbind, res[vapply(res, nrow, integer(1)) > 0])
  if (is.null(res_df)) {
    res_df = data.frame(react = character(0), round = integer(0), obj = numeric(0),
                        model = character(0), stringsAsFactors = FALSE)
  }
  rownames(res_df) = NULL
  
  cat(Ind, " ", diet_type, ": ", length(unique(res_df$react)),
      " limiting metabolites across ", length(unique(res_df$model)), " producers\n", sep = '')
  
  out_dir = paste0("supplement_pred_", tolower(diet_type))
  dir.create(out_dir, showWarnings = FALSE)
  saveRDS(res_df, file.path(out_dir, paste0(Ind, "_needed_mets.rds")))
}

run_one('Breast')
run_one('Formula')
cat("Finished ", Ind, "\n", sep = '')
