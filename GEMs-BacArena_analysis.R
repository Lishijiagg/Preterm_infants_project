# ================================================================================================================== #
# Personalized in silico gut microbiome metabolic modelling of preterm infants
#
# Analysis script accompanying the manuscript. Everything reported in the Results and the
# Supplementary Material is produced here, in the order the manuscript presents it.
#
#   Part 0   Environment, parameters and helper functions
#   Part 1   Input data, GEM library, in silico media and arena construction
#   Part 2   Fidelity of the simulated communities                        Results 1
#   Part 3   Blood-culture isolates and gut abundance                     Results 1
#   Part 4   Baseline fermentation product production                     Results 2
#   Part 5   Nutritional conditions and HMO utilisation                   Results 3
#   Part 6   Probiotic intervention                                       Results 4
#   Part 7   Determinacy of the flux attribution                          Results 2, closing paragraph
#   Part 8   Personalised nutrient supplementation                        Results 5
#   Part 9   Combined probiotic and nutrient intervention                 Results 5
#   Part 10  Prebiotic and synbiotic intervention on formula              Results 6
#   Part 11  Robustness of the lactate results                            Results 6
#   Part 12  Supplement composition and cross-arm comparison              Results 5, supporting
#
# On the position of Part 7
#   The sensitivity analyses report on both the resident models and the seven introduced strains, so
#   they require the probiotic models already constrained to each medium by Part 6. They are therefore
#   placed after Part 6 rather than next to the Results section they support. Every other part follows
#   the order of the manuscript.
#
# Running this script
#   Parts 0 and 1 must be run first. After that each part depends only on the parts above it.
#   The simulations themselves are not run here. They are produced on the cluster by the scripts in
#   simulation/ and reduced to the per-infant tables this script reads by extraction/. See README.md.
#
# Session
#   R 4.3.2, BacArena 1.8.2, sybil with glpkAPI (simplex). Package versions in sessionInfo.txt.
# ================================================================================================================== #

# ================================================================================================================== #
# Part 0  Environment, parameters and helper functions
#
#   Manuscript: setup
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Environment. Both the personal and the system library are on the search path and each holds its own copy of
# several packages, so tidyverse is attached first to fix the ggplot2 version before any dependent package can pull
# in a different one. Packages re-exported by tidyverse are not attached separately, and plyr is not attached at
# all because it masks several dplyr verbs; where a plyr function is required it is called as plyr::function().
# ---------------------------------------------------------------------------------------------------------------- #
# Run this script from the root of the repository. Uncomment and edit if you work from elsewhere.
# setwd("/path/to/Preterm_infants_project")
# .libPaths(new = "/path/to/your/R/library")
library(tidyverse)
library(BacArena)
library(R.matlab)
library(vegan)
library(ape)
library(patchwork)
library(ggnewscale)
library(ggtext)
library(ggforce)
library(ggpubr)
library(ggsignif)
library(ggsci)
library(ggvenn)
library(RColorBrewer)
library(introdataviz)
library(rstatix)
library(car)
library(ggh4x)
library(ggraph)
library(igraph)



# ---------------------------------------------------------------------------------------------------------------- #
# Global parameters. Every constant used in more than one place is defined here so that a change propagates through
# the whole analysis; nothing below redefines any of them.
# ---------------------------------------------------------------------------------------------------------------- #
group_levels = c('S-Preterm','N-Preterm','Term')
preterm_levels = c('S-Preterm','N-Preterm')
group_colors = c("S-Preterm" = "red", "N-Preterm" = "orange", "Term" = "blue")
diet_colors = c("Breast" = "deepskyblue1", "Formula" = "#00b32d")
status_colors = c("Non-treated" = "#00AFBB", "Probiotics" = "#E7B800", "Supplemented" = "#E7B800")
class_colors = c("Enterobacteriaceae" = "#ff7f00",
                 "Bacteroidetes" = "#1f78b4",
                 "Bifidobacterium" = "#e31a1c",
                 "Other" = "grey70")

arena_n = 100                 # arena side length in grid cells
sim_time = 12                 # simulation steps; substrate depletion and spatial saturation occur by step 10-11
diet_scale = 20               # diet fluxes are divided by this to obtain arena substrate concentrations
death_rate = 0.1
pro_amount = 10               # initial abundance of each introduced probiotic strain
supp_conc = 0.05              # mM, concentration at which each predicted supplement is supplied
utiliser_cutoff = 5           # percent of an HMO consumed above which an individual counts as a utiliser
resp_cutoff = 1e-5            # increase in log(mM+1) above which an individual counts as a responder
supplement_round = 1          # reduced-cost round retained for the supplement definition
exclude_mets = c("EX_o2(e)")  # not administrable as a dietary supplement

# Figure output directory. Created here so that the ggsave calls below work in a fresh clone of the
# repository, where the directory is not under version control.
dir.create('Figures', showWarnings = FALSE, recursive = TRUE)

# SCFA exchange reactions and the standardisation of the names BacArena returns for them. Each acid comes back under
# several labels because the source models differ in how they spell metNames.
scfas = c("EX_ac(e)","EX_ppa(e)","EX_succ(e)","EX_lac_D(e)","EX_lac_L(e)","EX_for(e)","EX_but(e)")
# Only acetate, propionate and butyrate are short-chain fatty acids under every definition; formate (C1) is
# included among them by some. Succinate and the two lactate isomers are not fatty acids at all, but they are
# retained because they are the central intermediates of the cross-feeding network measured here — succinate is
# the precursor of propionate via the succinate pathway and lactate the precursor of butyrate. The collective term
# used throughout is therefore "fermentation products" rather than "SCFA".
scfa_levels = c("Acetate","Propionate","Butyrate","Formate","Succinate","(R)-Lactate","(S)-Lactate")
scfa_class = c("Acetate" = "SCFA", "Propionate" = "SCFA", "Butyrate" = "SCFA", "Formate" = "SCFA",
               "Succinate" = "Organic acid", "(R)-Lactate" = "Organic acid", "(S)-Lactate" = "Organic acid")
scfa_class_levels = c("SCFA", "Organic acid")
name_map = c('acetate' = 'Acetate','Acetate' = 'Acetate',
             'butyrate' = 'Butyrate','Butyrate' = 'Butyrate','but' = 'Butyrate',
             'propionate' = 'Propionate','Propionate' = 'Propionate',
             'succinate' = 'Succinate','Succinate' = 'Succinate',
             'formate' = 'Formate','Formate' = 'Formate',
             '(S)-lactate' = '(S)-Lactate','L-lactate' = '(S)-Lactate',
             '(R)-lactate' = '(R)-Lactate','D-lactate' = '(R)-Lactate')
scfa_colors = c("EX_ac(e)" = "#6633FF", "EX_but(e)" = "#00CC99", "EX_for(e)" = "#FF3333",
                "EX_lac_D(e)" = "#CC00CC", "EX_lac_L(e)" = "#FF6699",
                "EX_succ(e)" = "#FFCC00", "EX_ppa(e)" = "#FF9570")
scfa_flux_labels = c("EX_ac(e)" = "Acetate", "EX_but(e)" = "Butyrate", "EX_for(e)" = "Formate",
                     "EX_lac_D(e)" = "(R)-Lactate", "EX_lac_L(e)" = "(S)-Lactate",
                     "EX_succ(e)" = "Succinate", "EX_ppa(e)" = "Propionate")

# HMO exchange reactions, with the same two-label problem. Eight of the thirteen oligosaccharide exchanges carried
# by the breast medium are quantified: the six most frequently reported HMOs plus two further fucosylated species
# for which the model library provides comparable coverage (22 and 24 of the 191 modelled strains). The remaining
# five are represented by three or fewer models, or, in the case of free fucose, are degradation products rather
# than intact oligosaccharides. Ordered here by structural class: fucosylated, neutral core, then sialylated.
hmo_reactions = c("2'-FL"  = "EX_2fuclac(e)",     "3'-FL"  = "EX_3fuclac(e)",
                  "DFL"    = "EX_dfuclac(e)",     "LDFT"   = "EX_lacdfucttr(e)",
                  "LNT"    = "EX_lacnttr(e)",     "LNnT"   = "EX_lacnnttr(e)",
                  "6'-SL"  = "EX_6slac(e)",       "3'-SL"  = "EX_3slac(e)")
hmo_levels = names(hmo_reactions)
hmo_name_map = c('2-Fucosyllactose'    = "2'-FL",  '2fuclac'    = "2'-FL",
                 '3-Fucosyllactose'    = "3'-FL",  '3fuclac'    = "3'-FL",
                 'Difucosyllactose'    = "DFL",    'dfuclac'    = "DFL",
                 'Lactodifucotetraose' = "LDFT",   'lacdfucttr' = "LDFT",
                 'Lacto-N-tetraose'    = "LNT",    'lacnttr'    = "LNT",
                 'Lacto-N-neotetraose' = "LNnT",   'lacnnttr'   = "LNnT",
                 '6-Sialyllactose'     = "6'-SL",  '6slac'      = "6'-SL",
                 '3-Sialyllactose'     = "3'-SL",  '3slac'      = "3'-SL")

# ---------------------------------------------------------------------------------------------------------------- #
# Paths. The input tables are shipped with this repository; the metabolic models are not, because they belong to
# AGORA2 (https://www.vmh.life) and to the HMO-extended reconstructions of Shaaban et al. Point the four model
# directories at a local copy of those libraries. data/mapping/ASV_GEM_map_with_genus.txt gives the model assigned
# to each ASV, which is sufficient to reassemble exactly the set used here.
# ---------------------------------------------------------------------------------------------------------------- #
data_dir  = "data"
gem_dir   = "<path to the 191 renamed .mat models used in the simulations>"
hmo_dir   = "<path to the Shaaban HMO-extended .mat models>"
agora_dir = "<path to the AGORA2 .mat models>"
old_dir   = "<path to the seven probiotic .mat models>"

# Probiotic strains and their short labels
probiotics = c(
  PRO_REU  = file.path(old_dir, "PRO_1_Lactobacillus_reuteri_SD2112_ATCC_55730.mat"),
  PRO_BB12 = file.path(old_dir, "PRO_2_Bifidobacterium_animalis_lactis_BB_12.mat"),
  PRO_LGG  = file.path(old_dir, "PRO_3_Lactobacillus_rhamnosus_GG_ATCC_53103.mat"),
  PRO_FER  = file.path(old_dir, "PRO_4_Lactobacillus_fermentum_IFO_3956.mat"),
  PRO_CBU  = file.path(old_dir, "PRO_5_Clostridium_butyricum_DSM_10702.mat"),
  PRO_BINF = file.path(hmo_dir, "Bifidobacterium_longum_infantis_ATCC_15697.mat"),
  PRO_BBIF = file.path(hmo_dir, "Bifidobacterium_bifidum_BGN4.mat"))
pro_labels = c(PRO_REU = "L.reuteri", PRO_BB12 = "B.lactis", PRO_LGG = "L.rhamnosus",
               PRO_FER = "L.fermentum", PRO_CBU = "C.butyricum",
               PRO_BINF = "B.infantis", PRO_BBIF = "B.bifidum")
pro_names = names(probiotics)   # model identifiers of the introduced strains


# ---------------------------------------------------------------------------------------------------------------- #
# Helper functions. Collected here rather than beside their first use so that the analysis sections below contain
# only the steps specific to each question.
# ---------------------------------------------------------------------------------------------------------------- #
# A single Eval object occupies well over a hundred gigabytes once loaded, so the simulation results are never held
# in memory here. extract_eval.R reduces each object on the cluster to four small tables (end-point SCFA, the HMO
# time course, the growth curve of every organism, and the SCFA exchange flux of every organism); the loader below
# reads those tables and binds them across subjects. Everything downstream works from these tables.
load_extracted = function(cond, expected = NULL) {
  dir_path = paste0("extracted_", cond)
  f = list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  f = f[basename(f) != "_mediac.rds"]
  parts = lapply(f, readRDS)
  out = lapply(c('scfa_end','hmo_curve','growth','flux'),
               function(k) do.call(rbind, lapply(parts, `[[`, k)))
  names(out) = c('scfa_end','hmo_curve','growth','flux')
  out$mediac = readRDS(file.path(dir_path, "_mediac.rds"))
  # A silently incomplete batch would propagate into every downstream table, so the count is checked here
  found = length(f)
  if (!is.null(expected) && found != expected) {
    warning(cond, ": ", found, " subjects extracted, ", expected, " expected", call. = FALSE)
  } else {
    message(cond, ": ", found, " subjects")
  }
  out
}

# Define a function for loading GEM files as the BacArena required format
readMATmod_fix <- function(file) {
  
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("Need R.matlab: install.packages('R.matlab')")
  }
  if (!requireNamespace("sybil", quietly = TRUE)) {
    stop("Need sybil: BacArena depends on it")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Need Matrix")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Need stringr")
  }
  
  print(system.time(data <- R.matlab::readMat(file)))
  dat.mat <- data[[1]]
  mod.var <- dimnames(dat.mat)[[1]]
  
  # ─── 1) model name ─────────────────────────────────────────
  if ("modelID" %in% mod.var) {
    mod.id <- as.character(dat.mat[[which(mod.var == "modelID")]])
  } else {
    mod.id <- names(data)[1]
  }
  if ("modelName" %in% mod.var) {
    mod.name <- as.character(dat.mat[[which(mod.var == "modelName")]])
  } else {
    mod.name <- mod.id
  }
  if ("description" %in% mod.var) {
    mod.desc <- as.character(dat.mat[[which(mod.var == "description")]])
  } else {
    mod.desc <- mod.id
  }
  
  # ─── 2) stoich matrix ──────────────────────────────────────
  mod.S <- Matrix::Matrix(dat.mat[[which(mod.var == "S")]], sparse = TRUE)
  
  # ─── 3) rxn ────────────────────────────────────────────────
  mod.react_id   <- unlist(dat.mat[[which(mod.var == "rxns")]])
  
  if ("rxnNames" %in% mod.var) {
    mod.react_name <- unlist(dat.mat[[which(mod.var == "rxnNames")]])
  } else {
    mod.react_name <- mod.react_id  # default: use the reaction id as the name
  }
  
  if ("rev" %in% mod.var) {
    mod.react_rev <- as.vector(dat.mat[[which(mod.var == "rev")]]) == TRUE
  } else {
    mod.react_rev <- as.vector(dat.mat[[which(mod.var == "lb")]]) < 0
  }
  
  # ─── 4) met ────────────────────────────────────────────────
  mod.met_id   <- unlist(dat.mat[[which(mod.var == "mets")]])
  if ("metNames" %in% mod.var) {
    mod.met_name <- unlist(dat.mat[[which(mod.var == "metNames")]])
  } else {
    mod.met_name <- mod.met_id
  }
  
  # ─── 5) genes ─────────────────
  if ("genes" %in% mod.var) {
    mod.gene_id <- unname(unlist(dat.mat[[which(mod.var == "genes")]]))
  } else {
    mod.gene_id <- character(0)
  }
  
  if ("rxnGeneMat" %in% mod.var) {
    mod.GeneMat <- dat.mat[[which(mod.var == "rxnGeneMat")]]
    mod.genes <- apply(mod.GeneMat, 1, function(row) {
      x <- unname(mod.gene_id[which(row != 0)])
      if (length(x) == 0) "" else x
    })
    mod.gpr <- sapply(sapply(dat.mat[[which(mod.var == "grRules")]], unlist),
                      function(entry) {
                        if (length(entry) == 0) "" else unname(entry)
                      })
  } else {
    mod.GeneMat <- Matrix::Matrix(0, 
                                  nrow = length(mod.react_id),
                                  ncol = length(mod.gene_id))
    mod.genes <- list()
    
    if ("rules" %in% mod.var) {
      rules.list <- unlist(dat.mat[[which(mod.var == "rules")]], recursive = FALSE)
      if (length(rules.list) != length(mod.react_id)) {
        stop("Length of rules not same as length of reactions.")
      }
      for (i in seq_along(mod.react_id)) {
        rule.tmp <- rules.list[[i]]
        if (length(rule.tmp) == 0) {
          mod.genes[[i]] <- ""
          next
        }
        j <- as.numeric(unlist(stringr::str_extract_all(
          rule.tmp, "(?<=x\\()[0-9]+?(?=\\))")))
        if (length(j) > 0) {
          mod.GeneMat[i, j] <- 1
          mod.genes[[i]] <- mod.gene_id[j]
        } else {
          mod.genes[[i]] <- ""
        }
      }
    } else {
      for (i in seq_along(mod.react_id)) {
        mod.genes[[i]] <- ""
      }
    }
    
    mod.gpr <- rep("", length(mod.react_id))
  }
  
  if ("rules" %in% mod.var) {
    mod.gprRules <- sapply(sapply(dat.mat[[which(mod.var == "rules")]], unlist),
                           function(entry) {
                             if (length(entry) == 0) ""
                             else {
                               numbers <- as.numeric(unlist(
                                 stringr::str_extract_all(entry, "[0-9]+")))
                               if (length(numbers) == 0) return("")
                               dict <- as.character(numbers - min(numbers) + 1)
                               names(dict) <- as.character(numbers)
                               gsub("\\(([0-9]+)\\)", "\\[\\1\\]",
                                    stringr::str_replace_all(entry, dict))
                             }
                           })
  } else {
    mod.gprRules <- sapply(seq_along(mod.gpr), function(i) {
      if (mod.gpr[i] == "") "" 
      else {
        genes <- unlist(mod.genes[i])
        dict <- paste0("x[", seq_along(genes), "]")
        names(dict) <- genes
        gsub("and", "&",
             gsub("or", "|",
                  stringr::str_replace_all(mod.gpr[i], dict)))
      }
    })
  }
  
  # ─── 6) bounds ─────────────────────────────────────────────
  mod.lb <- as.vector(dat.mat[[which(mod.var == "lb")]])
  mod.ub <- as.vector(dat.mat[[which(mod.var == "ub")]])
  
  # ─── 7) compartments ───────────────────────────────────────
  met_comp <- stringr::str_extract_all(mod.met_id, "(?<=\\[)[a-z](?=\\])")
  if (all(sapply(met_comp, length) == 0)) {
    met_comp <- stringr::str_extract_all(mod.met_id, "(?<=_)[a-z][0-9]?(?=$)")
  }
  mod.mod_compart <- unique(unlist(met_comp))
  mod.met_comp <- match(met_comp, mod.mod_compart)
  
  # ─── 8) subsystems ─────────────────────────────────────────
  if ("subSystems" %in% mod.var) {
    sub <- sapply(dat.mat[[which(mod.var == "subSystems")]], unlist)
    sub.unique <- unique(sub)
    mod.subSys <- Matrix::Matrix(FALSE, 
                                 nrow = length(sub),
                                 ncol = length(sub.unique), sparse = TRUE)
    for (i in seq_along(sub)) {
      j <- match(sub[i], sub.unique)
      if (!is.na(j)) mod.subSys[i, j] <- TRUE
    }
    colnames(mod.subSys) <- sub.unique
  } else {
    mod.subSys <- Matrix::Matrix(FALSE, nrow = length(mod.react_id),
                                 ncol = 1, sparse = TRUE)
  }
  
  model = sybil::modelorg(id = mod.id, name = mod.name)
  model@mod_desc = mod.desc
  model@S = mod.S
  model@lowbnd = mod.lb
  model@uppbnd = mod.ub
  model@met_id = mod.met_id
  model@met_name = mod.met_name
  model@met_num = length(mod.met_id)
  model@react_id = mod.react_id
  model@react_name = mod.react_name
  model@react_num = length(mod.react_id)
  model@react_rev = mod.react_rev
  model@genes = mod.genes
  model@gprRules = mod.gprRules
  model@gpr = mod.gpr
  model@allGenes = mod.gene_id
  model@mod_compart = mod.mod_compart
  model@met_comp = mod.met_comp
  model@subSys = mod.subSys
  
  obj.idx <- which(dat.mat[[which(mod.var == "c")]] != 0)
  if (length(obj.idx) > 0) {
    model <- sybil::changeObjFunc(model, react = obj.idx)
    print(sybil::optimizeProb(model))
  }
  
  return(model)
}

# Define a function GEMs_constraint
GEMs_constraint = function(GEMs_list,diet) {
  cons_GEMs_list = list()
  for (i in names(GEMs_list)) {
    # get nutrients information of nutrients
    diet_ex = as.character(diet[,1])
    diet_val = as.numeric(diet[,2])
    
    # find exchange reactions, their values, and their position in gem
    reaction_v = grep(pattern = 'EX_',GEMs_list[[i]]@react_id,value = TRUE, fixed = TRUE)
    reaction_p = grep(pattern = 'EX_',GEMs_list[[i]]@react_id,value = TRUE, fixed = TRUE)
    
    # find the diet's exchange reactions in the gem and values
    diet_pos = which(diet_ex %in% reaction_v, arr.ind = TRUE)
    diet_ex = diet_ex[diet_pos]
    
    # set nutrient concentration as lower bounds
    diet_val = -1*diet_val[diet_pos]
    diet_val[abs(diet_val) > 1000] = -1000
    
    # remove original diet (default in GEMs)
    print(paste("The original optimized obj_value of ",i," is ",optimizeProb(GEMs_list[[i]])@lp_obj,sep = ""))
    model_new = changeBounds(model = GEMs_list[[i]],lb = rep(x = 0,length(reaction_v)),
                             ub = rep(x = 1000,length(reaction_v)),
                             react = reaction_v)
    
    # applying defined medium
    model_new = changeBounds(model = model_new,lb = diet_val,
                             ub = rep(x = 1000, length(diet_val)),
                             react = diet_ex)
    print(paste("The optimized obj_value of modified ",i," is ",optimizeProb(model_new)@lp_obj,sep = ""))
    cons_GEMs_list[[i]] = model_new
  }
  return(cons_GEMs_list)
}


# Generate the initial arena object for each individual.
# The empty arena consist of 10000 grids with each one can be occupied by one GEM (only one microorganism).
# Add troubleshoot nutrient to the arena as the in silico medium.
# Add GEMs with their scaled abundance to the arena.
# ---------------------------------------------------------------------------------------------------------------- #
generate_arena = function(comp_data,metadata,diet,map_data,bac_object){
  sample_id = levels(factor(metadata$ID))
  Ind_arena = list()
  for (Ind in sample_id) {
    temp_meta = subset(metadata,ID == Ind)
    temp_GEMs_amount = comp_data[rownames(temp_meta)]
    temp_GEMs_amount = temp_GEMs_amount[temp_GEMs_amount != 0,,drop = FALSE]
    Ind_arena[[Ind]] = Arena(n = 100,m=100)
    r = 0
    for (asv in rownames(temp_GEMs_amount)) {
      r = r+1
      temp_GEM_input = map_data[map_data$ASV_ID == asv,]$renamed_gem
      Ind_arena[[Ind]] = addOrg(object = Ind_arena[[Ind]],
                                specI = bac_object[[temp_GEM_input]],
                                amount = temp_GEMs_amount[asv,])
    }
    new_diet_ex = as.character(diet[,1])
    new_diet_val = as.numeric(diet[,2])
    for (met in 1:length(new_diet_ex)) {
      Ind_arena[[Ind]] = addSubs(object = Ind_arena[[Ind]],
                                 smax = ifelse(new_diet_val[met]>1e-3,new_diet_val[met]/20,new_diet_val[met]),
                                 mediac = new_diet_ex[met],
                                 unit = 'mM',addAnyway = TRUE)
    }
  }
  
  return(Ind_arena)
}


heat_plot = function(genus_data,metadata,data_type) {
  genus_data_long = genus_data %>%
    pivot_longer(-Genus,names_to = "ID",values_to = "Abundance") %>%
    left_join(metadata %>% select(ID,group),by = "ID") %>%
    mutate(Id = factor(ID,levels = unique(ID[order(group)])))
  
  heatmap_plot = ggplot(genus_data_long,aes(x = Genus, y = ID, fill = Abundance))+
    geom_tile()+
    scale_fill_gradient(low = "white",high = "deepskyblue1")+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8))+
    labs(title = paste(data_type," genera relative abundance composition",sep = ""))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          title = element_text(size = 8),
          legend.title = element_text(size = 10))+
    labs(fill = "Relative Abundance")
  
  group_bar = genus_data_long %>% distinct(ID,group)
  group_bar_plot = ggplot(group_bar,aes(x = 1,y = ID,fill = group))+
    geom_tile()+
    scale_fill_manual(values = c("S-Preterm" = "red",
                                 "N-Preterm" = "orange",
                                 "Term" = "blue"))+
    theme_void()+
    theme(legend.position = 'none',
          axis.text = element_blank(),legend.title = element_text(size = 10))+
    labs(fill = "Group")
  
  output_plot = group_bar_plot + heatmap_plot + 
    plot_layout(guides = "collect",widths = c(0.5,20))+
    plot_annotation(tag_levels = NULL)
  return(output_plot)
}


# ---------------------------------------------------------------------------------------------------------------- #
aggregate_to_genus <- function(comp_data, map) {
  asv2genus <- setNames(map$Genus, map$ASV_ID)
  genus_vec <- asv2genus[rownames(comp_data)]
  if (any(is.na(genus_vec))) {
    missing <- rownames(comp_data)[is.na(genus_vec)]
    keep <- !is.na(genus_vec)
    comp_data <- comp_data[keep, , drop = FALSE]
    genus_vec <- genus_vec[keep]
  }
  
  genus_abund <- rowsum(comp_data, group = genus_vec)
  genus_rel <- sweep(genus_abund, 2, colSums(genus_abund), "/")
  genus_rel_df <- as.data.frame(genus_rel)
  genus_rel_df$Genus <- rownames(genus_rel_df)
  
  return(genus_rel_df)
}

alpha_diversity = function(comp_data, sample_ids) {
  x = as.data.frame(t(comp_data))
  data.frame(SampleID = sample_ids,
             Shannon = as.numeric(sprintf("%0.4f", vegan::diversity(x, index = 'shannon', base = 2))),
             Simpson = as.numeric(sprintf("%0.4f", vegan::diversity(x, index = 'simpson')))) %>%
    left_join(meta %>% mutate(SampleID = ID) %>% select(SampleID, group), by = 'SampleID') %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    pivot_longer(cols = c(Shannon, Simpson), names_to = 'Index', values_to = 'Value')
}

extract_scfa = function(extracted, metadata, diet_type) {
  extracted$scfa_end %>%
    mutate(sub = ifelse(sub %in% names(name_map), name_map[sub], sub)) %>%
    group_by(ID, sub) %>%
    summarise(value = sum(value), .groups = 'drop') %>%
    complete(ID, sub, fill = list(value = 0)) %>%
    left_join(metadata %>% select(ID, group), by = 'ID') %>%
    mutate(group = factor(group, levels = group_levels),
           sub = factor(sub, levels = scfa_levels),
           Class = factor(scfa_class[as.character(sub)], levels = scfa_class_levels),
           value = log(value + 1),
           Diet = diet_type)
}

scfa_group_box = function(scfa_df, data_type, adjust = TRUE) {
  stat_res = scfa_df %>%
    group_by(sub) %>%
    rstatix::wilcox_test(value ~ group, p.adjust.method = 'none') %>%
    ungroup()
  
  if (adjust) {
    stat_res = stat_res %>%
      rstatix::adjust_pvalue(method = 'BH') %>%
      rstatix::add_significance('p.adj')
    label_col = 'p.adj.signif'
  } else {
    stat_res = stat_res %>% rstatix::add_significance('p')
    label_col = 'p.signif'
  }
  
  # Bracket heights computed within each product so that the free y scales are preserved
  panel_max = scfa_df %>%
    group_by(sub) %>%
    summarise(top = max(value, na.rm = TRUE), .groups = 'drop')
  
  # The panels are nested inside the compound class, so the bracket table must carry both facet variables
  stat_res = stat_res %>%
    left_join(panel_max, by = 'sub') %>%
    group_by(sub) %>%
    mutate(y.position = top * (1.05 + 0.12 * row_number())) %>%
    ungroup() %>%
    mutate(sub = factor(sub, levels = scfa_levels),
           Class = factor(scfa_class[as.character(sub)], levels = scfa_class_levels))
  
  ggplot(scfa_df, aes(x = group, y = value))+
    geom_jitter(width = 0.15, size = 0.4, alpha = 1, colour = 'grey25')+
    geom_boxplot(aes(fill = group), width = 0.6, alpha = 1, outlier.shape = NA)+
    # Nesting the panels inside the compound class keeps the distinction between the true short-chain fatty
    # acids and the related organic acids visible on the figure itself
    ggh4x::facet_nested_wrap(~ Class + sub, scales = 'free_y', nrow = 1)+
    scale_fill_manual(values = group_colors, name = "Group")+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))+
    stat_pvalue_manual(stat_res, label = label_col,
                       tip.length = 0.01, size = 6, bracket.size = 0.3)+
    theme_bw()+
    labs(title = paste(data_type, " predicted fermentation products", sep = ''),
         x = "", y = "Fermentation product, log(mM+1)")+
    theme(strip.text = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          panel.grid.minor = element_blank())
}

# Substrates limiting a strain's growth on a given medium. A negative reduced cost marks an exchange whose
# relaxation would raise growth, whether or not the metabolite is already in the medium; relieving those and
# re-solving exposes the constraints that only become binding once the first ones are lifted. Reproduced here for
# reference; predict_supplements.R runs the equivalent on the cluster.
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

# ================================================================================================================== #
# Part 1  Input data, GEM library, in silico media and arena construction
#
#   Manuscript: Methods
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Read OTU table and metadata of the preterm infants
# ---------------------------------------------------------------------------------------------------------------- #
count_data = read.table(file.path(data_dir, 'features/asv_count.tsv'),sep = '\t',check.names = F)

raw_data = read.table(file = file.path(data_dir, "features/asv_rel_191_final.tsv"),
                      row.names = 1,header = TRUE,comment.char = '',
                      check.names = FALSE,sep = '\t')

meta = read.table(file.path(data_dir, "metadata/metadata.tsv"),header = TRUE,row.names = 1,sep = '\t')

ASV_GEM_map = read.table(file.path(data_dir, 'mapping/ASV_GEM_map_with_genus.txt'),sep = '\t')

all_amount_data = raw_data %>% 
  `*`(500) %>%  # Scale the relative abundance by multiplying by 500
  round(digits = 0) %>%   # Round to integers
  filter(rowSums(.) > 1)
#write.table(all_amount_data,file = 'all_amount_data.txt',row.names = TRUE,col.names = TRUE,sep = '\t')
all_amount_data = read.table(file.path(data_dir, 'features/all_amount_data.txt'),sep = '\t',check.names = F)

# The 191 renamed models, loaded once and cached; predict_supplements.R reads the cached copy
GEMs_files = list.files(path = gem_dir,pattern = "\\.mat$",full.names = TRUE)
# Create an empty list to store GEMs files
GEMs_data_list = list()
# Iterate each file in GEMs_files and load them.
for (mat_file in GEMs_files) {
  GEMs_data <- readMATmod_fix(mat_file)
  bname <- basename(mat_file)
  species_name <- sub("^ASV\\d+_", "", bname)
  species_name <- sub("\\.mat$", "", species_name)
  
  GEMs_data@mod_id   <- species_name
  GEMs_data@mod_name <- species_name
  GEMs_data@mod_desc <- species_name
  
  GEMs_data_list[[bname]] <- GEMs_data
}

#saveRDS(GEMs_data_list,file = "all_GEMs_raw.rds")


# Load the in silico nutrient data and troubleshoot it to ensure that all GEMs can be simulated to grow on it.
# The original nutrient data contains all compounds of mixed 50g breast and formula milk, and was obtained from 
# the VMH database.
# The nutrient data must contain at least 2 columns, which are the substrate exchange reaction IDs (EXCHANGE.REACTION) 
# and the initial concentrations of each substrate (MOLECULAR.CONCENTRATION.IN.TOTAL).
# ---------------------------------------------------------------------------------------------------------------- #
# Formula diet
formula_medium = read.table(file = file.path(data_dir, 'media/formula_diet.tsv'),sep = '\t')

# Breast milk diet
breast_medium = read.table(file.path(data_dir, 'media/breast_diet.tsv'),sep = '\t')

# Define a function for adding the found compounds until there is growth
GM_compounds_added_until_growth <- function(model, cutoff_for_growth) {
  vector.exchanges <- vector() # initial vector to save the exchange reactions
  mpd <- model # model
  vec.reactions <- grep(pattern = "EX_", x = mpd@react_id, fixed = T,
                        value = T) # report all exchange reactions
  vec.reactions.positions <- grep(pattern = "EX_",
                                  x = mpd@react_id, fixed = T,
                                  value = F)  # report positions of all exchangre reactions
  # optimization results - as seen in example of sybil documentation
  x2 <- optimizeProb(mpd, poCmd = list("getRedCosts"))
  vd <- vector()  # interim vector to save the found exchangw reactions
  while (x2@lp_obj <= cutoff_for_growth) { # 1e-6 cut-off value
    x3 <- postProc(x2) # as seen in example of sybil documentation
    x4 <- x3@pa[[1]] # as seen in example of sybil documentation
    # make a dataframe to save the results - formatting
    df <- matrix(nrow = length(vec.reactions), ncol = 2)
    df[, 1] <- vec.reactions
    df1 <- as.data.frame(df)
    # only the reduced costs of the exchange reactions
    df1$V2 <- as.numeric(x4[vec.reactions.positions])
    df2 <- df1[abs((df1[, 2])) > (1e-6), ] # default cutoff
    df2 <- as.data.frame(df2)
    # alphabetically ordered
    dfmin <- sort(df2$V1[which(df2$V2 == min(df2$V2))])[[1]]
    #optimization - preparing the results for the next round of the loop
    #add -1 (assumption)
    mpd <- changeBounds(model = mpd, lb = -1, ub = 1000,
                        react = as.character(dfmin))
    x2 <- optimizeProb(mpd, poCmd = list("getRedCosts"))
    vd <- unique(c(vd, as.character(dfmin)))
  }
  #cat("These compounds have been added until growth achieved: \n")
  #print(vd)
  #initial vector to save the exchange reactions
  ve <- unique(c(vector.exchanges, vd))
  return(ve)
}

# Define a function to test what nutrients are essential for cell growth, and 
# add those nutrients to the last row of the human milk file.
addSubs_to_def_medium <- function(model,medium_comp) {
  print(model@mod_name)
  #get nutrients infomation of nutrints
  diet_ex <- as.character(medium_comp[,1])
  diet_val <- as.numeric(medium_comp[,2])
  #find exchanges reactions, their values, and their position in model
  reactions_v = grep(pattern = "EX_",model@react_id,value = TRUE,fixed = TRUE)
  reactions_p = grep(pattern = "EX_",model@react_id,value = FALSE,fixed = TRUE)
  #find the diet's exchange reactions in the model and their values
  diet_pos = which(diet_ex %in% reactions_v, arr.ind = TRUE)
  diet_ex = diet_ex[diet_pos]
  diet_val = -1*diet_val[diet_pos]
  diet_val[diet_val < -1000] = -1000
  #remove original diet (default in GEMs)
  model_new = changeBounds(model = model,
                           lb = rep(x = 0,length(reactions_v)),
                           ub = rep(x = 1000,length(reactions_v)),
                           react = reactions_v)
  #applying defined medium
  model_new = changeBounds(model = model_new,
                           lb = diet_val,
                           ub = rep(x = 1000,length(diet_val)),
                           react = diet_ex)
  #optimization with the new setting
  opt_results = optimizeProb(model_new)
  if (opt_results@lp_obj > 1e-6) {
    cat("This GEM can grow on the defined medium: \n")
    print(model_new@mod_name)
  } else {
    model_new_df_g <- GM_compounds_added_until_growth(model = model_new,cutoff_for_growth = 1e-6)
    cat("For model",model@mod_name," -- need substrates for growth: \n")
    print(model_new_df_g)
    #add needed substrates into the defined medium data as the last row
    for (i in model_new_df_g) {
      new_row <- data.frame(
        EXCHANGE.REACTION = as.character(i),
        MOLECULAR.CONCENTRATION.IN.TOTAL = 0.01,
        UNIT = 'mM')
      medium_comp = rbind(medium_comp,new_row)
    }
  }
  return(medium_comp)
}

new_formula = formula_medium
for (model in names(GEMs_data_list)) {
  new_formula = addSubs_to_def_medium(model = GEMs_data_list[[model]],
                                      medium_comp = new_formula)
}
#write.table(new_formula,file = 'new_formula_diet.txt',row.names = TRUE,col.names = TRUE,sep = '\t')
new_formula = read.table(file.path(data_dir, "media/new_formula_diet.txt"),sep = '\t')

new_breast = breast_medium
for (model in names(GEMs_data_list)) {
  new_breast = addSubs_to_def_medium(model = GEMs_data_list[[model]],
                                     medium_comp = new_breast)
}
#write.table(new_breast,file = 'new_breast_diet.txt',row.names = TRUE,col.names = TRUE,sep = '\t')
new_breast = read.table(file.path(data_dir, 'media/new_breast_diet.txt'),sep = '\t')

# The HMO-supplemented formula medium, built by simulation/prepare_hmo_formula_diet.R, which transfers the
# thirteen HMO exchange reactions of the breast medium onto the gap-filled formula medium. It is read here
# rather than rebuilt so that the sensitivity analysis in Part 7 and the factorial design in Part 10 use the
# same medium.
new_formula_hmo = read.table(file.path(data_dir, 'media/new_formula_hmo_diet.txt'), sep = '\t')


# ---------------------------------------------------------------------------------------------------------------- #
# According to the troubleshoot nutrient data, constraint exchange reaction lowerbounds to mimic the simulated 
# environment for each GEM .
# The available substrate exchange reaction lowerbounds were set to -1*initial concentration, and the substrate 
# non-exit in the in silico medium (nutrient) were set to 0.

# ---------------------------------------------------------------------------------------------------------------- #
# Constrain every model with each diet, build the Bac objects and the per-individual arenas. The arenas are
# reproduced here for reference; the cluster scripts perform the equivalent steps one subject at a time.
# ---------------------------------------------------------------------------------------------------------------- #
formula_con_GEMs_list = GEMs_constraint(GEMs_list = GEMs_data_list,diet = new_formula)
#saveRDS(formula_con_GEMs_list,file = 'formula_all_GEMs_needed.rds')
formula_con_GEMs_list = readRDS(file.path(data_dir, "derived/formula_all_GEMs_needed.rds"))

breast_con_GEMs_list = GEMs_constraint(GEMs_list = GEMs_data_list,diet = new_breast)
#saveRDS(breast_con_GEMs_list,file = 'breast_all_GEMs_needed.rds')
breast_con_GEMs_list = readRDS(file.path(data_dir, "derived/breast_all_GEMs_needed.rds"))


formula_add_bacs = list()
for (i in names(formula_con_GEMs_list)) {
  formula_add_bacs[[i]] <- Bac(model = formula_con_GEMs_list[[i]],
                               limit_growth = TRUE,  # Set limited growth for each GEM
                               deathrate = 0.1,  # Set 10% deathrate for each GEM during simulations
                               minweight = runif(1,min = 0.2,max = 0.4), # Set the minimum initial biomass (unit pg)
                               maxweight = runif(1,min = 0.8,max = 1.0), # Set the maximum initial biomass (unit pg)
                               growtype = "exponential") # The growth type should be exponential
}

breast_add_bacs <- list()
for (i in names(breast_con_GEMs_list)) {
  breast_add_bacs[[i]] <- Bac(model = breast_con_GEMs_list[[i]],
                              limit_growth = TRUE,  # Set limited growth for each GEM
                              deathrate = 0.1,  # Set 10% deathrate for each GEM during simulations
                              minweight = runif(1,min = 0.2,max = 0.4), # Set the minimum initial biomass (unit pg)
                              maxweight = runif(1,min = 0.8,max = 1.0), # Set the maximum initial biomass (unit pg)
                              growtype = "exponential") # The growth type should be exponential
}



formula_Ind_arena = generate_arena(comp_data = all_amount_data, metadata = meta,
                                   diet = new_formula, map_data = ASV_GEM_map,
                                   bac_object = formula_add_bacs)
breast_Ind_arena = generate_arena(comp_data = all_amount_data, metadata = meta,
                                  diet = new_breast, map_data = ASV_GEM_map,
                                  bac_object = breast_add_bacs)

preterm_ids = meta$ID[meta$group %in% preterm_levels]
writeLines(preterm_ids, 'id_list_preterm.txt')

# ---------------------------------------------------------------------------------------------------------------- #
# Taxonomic features of the measured communities. These are derived from the sequencing data alone and are referred
# to by several later sections, so they are computed once here rather than beside their first use. Two abundance
# scales are used throughout: the measured relative abundance from the sequencing table, and the scaled inoculum
# that actually enters a simulated community, which is the relative abundance multiplied by 500 and rounded. A taxon
# below roughly 0.2% in a given infant is present in the sequencing data but absent from that infant's simulation,
# so the two scales are reported separately wherever the distinction matters.
# ---------------------------------------------------------------------------------------------------------------- #
asv_rel = read.table(file.path(data_dir, "features/asv_rel.tsv"),
                     header = TRUE, row.names = 1, sep = "\t",
                     skip = 1, comment.char = "", check.names = FALSE)
asv_rel = asv_rel[rowSums(asv_rel) > 0, ]

# Hash-to-ASV-ID correspondence for the 268 retained features; asv_rel is keyed on the hashes and the mapping
# tables on the renamed IDs, so both are needed to look up a measured relative abundance by genus
asv_id_mapping = read.table(file.path(data_dir, "features/asv_id_mapping.tsv"),
                            header = TRUE, sep = "\t", check.names = FALSE)

# SILVA annotation of the 268 retained features, including the taxa that could not be matched to a model. It is
# needed wherever a genus has to be located in the sequencing data rather than in the simulated community.
asv_taxonomy = read.table(file.path(data_dir, "features/asv_268_with_taxonomy.tsv"),
                          header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Relative abundance of a set of genera per infant, keyed by simulation ID
genus_rel_abund = function(genera, map_data = ASV_GEM_map, id_map = asv_id_mapping) {
  # SILVA uses bracketed names for some clades, so the genus match is widened by pattern
  matched = grep(paste(genera, collapse = "|"), unique(map_data$Genus), value = TRUE)
  asv = map_data$ASV_ID[map_data$Genus %in% matched]
  hash = id_map$asv_hash[id_map$new_id %in% asv]
  setNames(colSums(asv_rel[rownames(asv_rel) %in% hash, rownames(meta), drop = FALSE]),
           meta$ID)
}

# Scaled inoculum of a set of genera per infant, keyed by simulation ID
genus_inoculum = function(genera, map_data = ASV_GEM_map) {
  matched = grep(paste(genera, collapse = "|"), unique(map_data$Genus), value = TRUE)
  asv = map_data$ASV_ID[map_data$Genus %in% matched]
  setNames(colSums(all_amount_data[rownames(all_amount_data) %in% asv,
                                   rownames(meta), drop = FALSE]),
           meta$ID)
}

# Genera present in a given infant's simulated community, from a set of candidates
genera_in_community = function(genera, map_data = ASV_GEM_map) {
  matched = grep(paste(genera, collapse = "|"), unique(map_data$Genus), value = TRUE)
  asv = map_data$ASV_ID[map_data$Genus %in% matched]
  vapply(rownames(meta), function(s) {
    present = asv[all_amount_data[asv, s] > 0]
    g = unique(map_data$Genus[map_data$ASV_ID %in% present])
    if (length(g) == 0) "none" else paste(sort(g), collapse = "; ")
  }, character(1), USE.NAMES = FALSE)
}

bif_genera = "Bifidobacterium"
entero_genera = c("Escherichia","Enterobacter","Klebsiella","Citrobacter","Serratia","Yersinia")
but_genera = c("Roseburia","Faecalibacterium","Eubacterium","Anaerostipes","Subdoligranulum")
bac_genera = c("Bacteroides","Parabacteroides","Prevotella")

# ---- Per-infant taxonomic summary --------------------------------------------------------------------------------

taxa_features = data.frame(
  ID = meta$ID,
  group = factor(meta$group, levels = group_levels),
  total_inoculum = colSums(all_amount_data[, rownames(meta), drop = FALSE]),
  bif_rel = 100 * genus_rel_abund(bif_genera),
  bif_inoculum = genus_inoculum(bif_genera),
  entero_inoculum = genus_inoculum(entero_genera),
  but_rel = 100 * genus_rel_abund(but_genera),
  but_inoculum = genus_inoculum(but_genera),
  but_genera_present = genera_in_community(but_genera),
  bac_rel = 100 * genus_rel_abund(bac_genera),
  row.names = NULL) %>%
  mutate(entero_pct = 100 * entero_inoculum / total_inoculum)

# Kept for compatibility with code written before taxa_features existed
bif_abund = setNames(taxa_features$bif_inoculum, taxa_features$ID)

# ---- Enterobacteriaceae ------------------------------------------------------------------------------------------
entero_summary = taxa_features %>%
  group_by(group) %>%
  summarise(n = n(),
            median_pct = median(entero_pct),
            q25 = quantile(entero_pct, .25),
            q75 = quantile(entero_pct, .75),
            median_bif = median(bif_inoculum),
            n_bif_zero = sum(bif_inoculum == 0), .groups = 'drop')

entero_kruskal = taxa_features %>% rstatix::kruskal_test(entero_pct ~ group)
entero_stat = taxa_features %>%
  rstatix::wilcox_test(entero_pct ~ group) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

# ---- Butyrate producers ------------------------------------------------------------------------------------------
# The universal zero butyrate prediction rests on these taxa being absent from almost every simulated community.
# Sequencing detection and simulated presence differ in the Term group, where the taxa are widely detected but at
# abundances too low to survive the scaling, so both are reported.
but_producer_summary = taxa_features %>%
  group_by(group) %>%
  summarise(n = n(),
            median_rel = median(but_rel),
            n_rel_zero = sum(but_rel == 0),
            median_inoculum = median(but_inoculum),
            max_inoculum = max(but_inoculum),
            n_inoculum_zero = sum(but_inoculum == 0), .groups = 'drop')

supp_table_butyrate = taxa_features %>%
  select(ID, Group = group,
         `Measured relative abundance (%)` = but_rel,
         `Scaled inoculum` = but_inoculum,
         `Genera in simulated community` = but_genera_present) %>%
  mutate(`Measured relative abundance (%)` = round(`Measured relative abundance (%)`, 4)) %>%
  arrange(Group, desc(`Measured relative abundance (%)`))

write.table(supp_table_butyrate, "Supplementary_Table_S9.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)


# ---------------------------------------------------------------------------------------------------------------- #
# Read the baseline simulation results
# ---------------------------------------------------------------------------------------------------------------- #
breast_ext  = load_extracted('no_treat_breast', expected = 51)
formula_ext = load_extracted('no_treat_formula', expected = 51)

# ================================================================================================================== #
# Part 2  Fidelity of the simulated communities
#
#   Manuscript: Results 1; Figure 1; Figures S1, S2
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Growth curve of every organism in each community, taken from the extracted tables
# ---------------------------------------------------------------------------------------------------------------- #
growth_curve_list = function(growth_df) {
  split(growth_df, growth_df$ID)
}
growth_curve_plots = function(growth_df) {
  lapply(split(growth_df, growth_df$ID), function(d) {
    ggplot(d, aes(x = time, y = value, color = species))+
      geom_line(linewidth = 1.1)+
      labs(title = paste(unique(d$ID), "gut microbiota growth curves", sep = " "),
           x = "Time (1 hour/sim_step)", y = "Biomass")+
      theme(text = element_text(size = 17))
  })
}

breast_Growth_curves  = growth_curve_list(breast_ext$growth)
formula_Growth_curves = growth_curve_list(formula_ext$growth)
breast_Growth_curves_plot  = growth_curve_plots(breast_ext$growth)
formula_Growth_curves_plot = growth_curve_plots(formula_ext$growth)

# ---------------------------------------------------------------------------------------------------------------- #
# Obtain species composition data for all individuals after 12 hours of simulation.
# Obtain the abundance composition of GEMs in the in silico microbiota of each individual at the last simulation time step.
# ---------------------------------------------------------------------------------------------------------------- #
breast_sim_df_list = list()
for (Ind in names(breast_Growth_curves)) {
  temp_data = subset(breast_Growth_curves[[Ind]],time == 12)
  colnames(temp_data)[3] = Ind
  breast_sim_df_list[[Ind]] = temp_data %>% select(1,3)
}

formula_sim_df_list = list()
for (Ind in names(formula_Growth_curves)) {
  temp_data = subset(formula_Growth_curves[[Ind]],time == 12)
  colnames(temp_data)[3] = Ind
  formula_sim_df_list[[Ind]] = temp_data %>% select(1,3)
}


# ---------------------------------------------------------------------------------------------------------------- #
# Combine each individual's simulated species composition at the last simulated time step to generate the 
# new simulated species composition table.
# ---------------------------------------------------------------------------------------------------------------- #
breast_sim_table = breast_sim_df_list[[1]]
for (i in 2:length(breast_sim_df_list)) {
  breast_sim_table = full_join(breast_sim_table,breast_sim_df_list[[i]],by = "species")
}
breast_sim_table[is.na(breast_sim_table)] = 0
breast_sim_table = breast_sim_table %>% mutate(
  species = species %>%
    str_replace_all(",","") %>%
    str_replace_all("\\.","") %>%
    str_replace_all("_"," ") %>%
    str_replace_all("-"," ")
)

formula_sim_table = formula_sim_df_list[[1]]
for (i in 2:length(formula_sim_df_list)) {
  formula_sim_table = full_join(formula_sim_table,formula_sim_df_list[[i]],by = "species")
}
formula_sim_table[is.na(formula_sim_table)] = 0
formula_sim_table = formula_sim_table %>% mutate(
  species = species %>%
    str_replace_all(",","") %>%
    str_replace_all("\\.","") %>%
    str_replace_all("_"," ") %>%
    str_replace_all("-"," ")
)


# ---------------------------------------------------------------------------------------------------------------- #
# Manually convert the species column in sim_table to the species annotation column in the common OTU table.
# Set the row name of sim_table to the OTU ID corresponding to GEMs, and rename it to sim_data for subsequent analysis.
# ---------------------------------------------------------------------------------------------------------------- #
breast_sim_Data = breast_sim_table %>% 
  left_join(ASV_GEM_map %>% select(strain, ASV_ID,Genus),
            by = c("species" = "strain")) %>% column_to_rownames(var = 'ASV_ID')
breast_sim_Genus = breast_sim_Data$Genus
breast_sim_data = breast_sim_Data[,-c(1,ncol(breast_sim_Data))]

formula_sim_Data = formula_sim_table %>% 
  left_join(ASV_GEM_map %>% select(strain, ASV_ID,Genus),
            by = c("species" = "strain")) %>% column_to_rownames(var = 'ASV_ID')
formula_sim_Genus = formula_sim_Data$Genus
formula_sim_data = formula_sim_Data[,-c(1,ncol(formula_sim_Data))]

# ---------------------------------------------------------------------------------------------------------------- #
# Shared index vectors for this part. The measured table is keyed on sequencing IDs and the simulated tables on
# simulation IDs, so both orderings are kept; common_asv restricts every comparison to the features present in all
# three datasets, and sample_levels orders the subjects by group for the figures.
# ---------------------------------------------------------------------------------------------------------------- #
selected_all_data = subset(count_data, rownames(count_data) %in% rownames(breast_sim_data))
common_asv = Reduce(intersect, list(rownames(selected_all_data),
                                    rownames(breast_sim_data),
                                    rownames(formula_sim_data)))
sample_levels = meta$ID[order(meta$group)]
sample_levels = sample_levels[sample_levels %in% colnames(breast_sim_data)]
seq_levels = rownames(meta)[match(sample_levels, meta$ID)]   # matching sequencing IDs of the measured table
n = length(sample_levels)



# ============================================================================
# Beta Diversity Analysis of Gut Microbiota in Three Groups of Preterm Infants (NMDS Version)
# Original Sequencing Microbiota vs. BacArena Simulated Microbiota
#
# Workflow: Bray–Curtis distance → NMDS ordination → PERMANOVA → betadisper → ggplot
# ============================================================================

raw_comm <- as.data.frame(t(selected_all_data))

raw_meta_aligned <- meta[match(rownames(raw_comm), rownames(meta)), ]

raw_bray_dist <- vegdist(raw_comm, method = "bray")

set.seed(123)
raw_nmds <- metaMDS(raw_comm, distance = "bray",
                    autotransform = FALSE,
                    try = 50, trymax = 200,
                    trace = FALSE)

set.seed(123)
raw_adonis <- as.data.frame(
  adonis2(raw_bray_dist ~ group, data = raw_meta_aligned, permutations = 999)
)
raw_p  <- raw_adonis$`Pr(>F)`[1]
raw_r2 <- raw_adonis$R2[1]

raw_bd_p <- anova(betadisper(raw_bray_dist, raw_meta_aligned$group))$`Pr(>F)`[1]

meta_join_raw <- meta
meta_join_raw$SampleID <- rownames(meta)

raw_result <- as.data.frame(raw_nmds$points) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
  dplyr::left_join(meta_join_raw[, c("SampleID", "group")], by = "SampleID") %>%
  dplyr::mutate(NMDS1 = as.numeric(NMDS1), NMDS2 = as.numeric(NMDS2))

beta_ori <- ggplot(raw_result, aes(NMDS1, NMDS2, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(alpha = 0.5) +
  scale_color_manual(values = group_colors, name = "Group") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 10)) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("PERMANOVA: italic(P) == ", signif(raw_p, 3), "*', ' ~ R^2 == ", signif(raw_r2, 3)),
           parse = TRUE, hjust = -0.1, vjust = 1.05, size = 5)+
  theme(legend.position = 'none')

d_ori = vegdist(t(selected_all_data[common_asv, seq_levels]), method = 'bray')
bd_ori = betadisper(d_ori, meta$group[match(sample_levels, meta$ID)])
anova(bd_ori)


# breast-feeding results
breast_sim_comm <- as.data.frame(t(breast_sim_data))

breast_sim_meta_aligned <- meta[match(rownames(breast_sim_comm), meta$ID), ]

breast_sim_bray_dist <- vegdist(breast_sim_comm, method = "bray")

set.seed(123)
breast_sim_nmds <- metaMDS(breast_sim_comm, distance = "bray",
                           autotransform = FALSE,
                           try = 50, trymax = 200,
                           trace = FALSE)

set.seed(123)
breast_sim_adonis <- as.data.frame(
  adonis2(breast_sim_bray_dist ~ group, data = breast_sim_meta_aligned, permutations = 999)
)
breast_sim_p  <- breast_sim_adonis$`Pr(>F)`[1]
breast_sim_r2 <- breast_sim_adonis$R2[1]

meta_join_sim <- meta
meta_join_sim$SampleID <- meta$ID

breast_sim_result <- as.data.frame(breast_sim_nmds$points) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
  dplyr::left_join(meta_join_sim[, c("SampleID", "group")], by = "SampleID") %>%
  dplyr::mutate(NMDS1 = as.numeric(NMDS1), NMDS2 = as.numeric(NMDS2))

breast_beta_sim <- ggplot(breast_sim_result, aes(NMDS1, NMDS2, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(alpha = 0.5) +
  scale_color_manual(values = group_colors, name = "Group") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 10)) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("PERMANOVA: italic(P) == ", signif(breast_sim_p, 3), 
                          "*', ' ~ R^2 == ", signif(breast_sim_r2, 3)),
           parse = TRUE, hjust = -0.1, vjust = 1.05, size = 4)

d_breast = vegdist(t(breast_sim_data[common_asv, sample_levels]), method = 'bray')
bd_breast = betadisper(d_breast, meta$group[match(sample_levels, meta$ID)])
anova(bd_breast)

# formula-feeding results
formula_sim_comm <- as.data.frame(t(formula_sim_data))

formula_sim_meta_aligned <- meta[match(rownames(formula_sim_comm), meta$ID), ]

formula_sim_bray_dist <- vegdist(formula_sim_comm, method = "bray")

set.seed(123)
formula_sim_nmds <- metaMDS(formula_sim_comm, distance = "bray",
                            autotransform = FALSE,
                            try = 50, trymax = 200,
                            trace = FALSE)

set.seed(123)
formula_sim_adonis <- as.data.frame(
  adonis2(formula_sim_bray_dist ~ group, data = formula_sim_meta_aligned, permutations = 999)
)
formula_sim_p  <- formula_sim_adonis$`Pr(>F)`[1]
formula_sim_r2 <- formula_sim_adonis$R2[1]

formula_sim_result <- as.data.frame(formula_sim_nmds$points) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
  dplyr::left_join(meta_join_sim[, c("SampleID", "group")], by = "SampleID") %>%
  dplyr::mutate(NMDS1 = as.numeric(NMDS1), NMDS2 = as.numeric(NMDS2))

formula_beta_sim <- ggplot(formula_sim_result, aes(NMDS1, NMDS2, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(alpha = 0.5) +
  scale_color_manual(values = group_colors, name = "Group") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 10)) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("PERMANOVA: italic(P) == ", signif(formula_sim_p, 3), 
                          "*', ' ~ R^2 == ", signif(formula_sim_r2, 3)),
           parse = TRUE, hjust = -0.1, vjust = 1.05, size = 4)

d_formula = vegdist(t(formula_sim_data[common_asv, sample_levels]), method = 'bray')
bd_formula = betadisper(d_formula, meta$group[match(sample_levels, meta$ID)])
anova(bd_formula)

# ---------------------------------------------------------------------------------------------------------------- #
# Genus-level composition of the measured and simulated communities
# ---------------------------------------------------------------------------------------------------------------- #

measured_genus  = aggregate_to_genus(selected_all_data, ASV_GEM_map)
id_map = setNames(meta$ID, rownames(meta))
colnames(measured_genus) = ifelse(colnames(measured_genus) %in% names(id_map),
                                  id_map[colnames(measured_genus)], colnames(measured_genus))
breast_genus    = aggregate_to_genus(breast_sim_data,  ASV_GEM_map)
formula_genus   = aggregate_to_genus(formula_sim_data, ASV_GEM_map)

measured_heatmap = heat_plot(genus_data = measured_genus,metadata = meta,data_type = 'Measured')+
  labs(title = NULL)+theme(legend.position = 'none')
breast_sim_heatmap = heat_plot(genus_data = breast_genus,metadata = meta,data_type = "Simulated (breastfed)")+
  labs(title = NULL)+theme(legend.position = 'none')
formula_sim_heatmap = heat_plot(genus_data = formula_genus,metadata = meta,data_type = "Simulated (formulafed)")+
  labs(title = NULL)+theme(legend.position = 'none')


# ---------------------------------------------------------------------------------------------------------------- #
# Genus-level composition of the measured and simulated communities, averaged within each subject group. Genera are
# ordered by their overall mean abundance and share a single palette across all panels so that the same genus keeps
# the same colour and stacking position in every group; genera outside the most abundant set are pooled as Others.
# ---------------------------------------------------------------------------------------------------------------- #
n_top_genus = 20

# Genera ranked by mean abundance across the whole cohort, computed once so that all panels share the same set
genus_rank = measured_genus %>%
  pivot_longer(-Genus, names_to = 'ID', values_to = 'Abundance') %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = 'drop') %>%
  arrange(desc(mean_abund))
top_genera = head(genus_rank$Genus, n_top_genus)
genus_levels = c(top_genera, 'Others')

genus_palette = setNames(
  c(colorRampPalette(brewer.pal(12, 'Paired'))(n_top_genus), 'grey70'),
  genus_levels)

genus_stack_data = function(genus_data, metadata, data_type) {
  genus_data %>%
    pivot_longer(-Genus, names_to = 'ID', values_to = 'Abundance') %>%
    mutate(Genus = ifelse(Genus %in% top_genera, Genus, 'Others')) %>%
    left_join(metadata %>% select(ID, group), by = 'ID') %>%
    group_by(group, ID, Genus) %>%
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    group_by(group, Genus) %>%
    summarise(Abundance = mean(Abundance), .groups = 'drop') %>%
    group_by(group) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    mutate(Genus = factor(Genus, levels = genus_levels),
           group = factor(group, levels = group_levels),
           Dataset = data_type)
}

genus_stack_all = rbind(
  genus_stack_data(measured_genus, meta, 'Measured'),
  genus_stack_data(breast_genus,   meta, 'Simulated, breast milk'),
  genus_stack_data(formula_genus,  meta, 'Simulated, formula')) %>%
  mutate(Dataset = factor(Dataset, levels = c('Measured',
                                              'Simulated, breast milk',
                                              'Simulated, formula')))

genus_stack_plot = ggplot(genus_stack_all, aes(x = group, y = Abundance, fill = Genus))+
  geom_bar(stat = 'identity', width = 0.7)+
  facet_wrap(~ Dataset, nrow = 1)+
  scale_fill_manual(values = genus_palette, name = "Genus")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0))+
  theme_bw()+
  labs(x = "", y = "Relative abundance")+
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 8, face = 'italic'),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid = element_blank())

cor_breast = sapply(1:n, function(i) cor(as.numeric(selected_all_data[common_asv, seq_levels[i]]),
                                         as.numeric(breast_sim_data[common_asv, sample_levels[i]]), method = 'spearman'))
cor_formula = sapply(1:n, function(i) cor(as.numeric(selected_all_data[common_asv, seq_levels[i]]),
                                          as.numeric(formula_sim_data[common_asv, sample_levels[i]]), method = 'spearman'))
cor_breast[is.na(cor_breast)] = 0
cor_formula[is.na(cor_formula)] = 0

fidelity_df = data.frame(Sample = sample_levels,
                         group = meta$group[match(sample_levels, meta$ID)],
                         Breast = cor_breast,
                         Formula = cor_formula,
                         Diff = cor_breast - cor_formula)
fidelity_test = wilcox.test(cor_breast, cor_formula, paired = TRUE)

# Vertices of the equilateral triangle and evenly spaced nodes along each edge
apex_top = c(x = 0.5, y = sqrt(3)/2); apex_left = c(x = 0, y = 0); apex_right = c(x = 1, y = 0)
tri_center = c(x = 0.5, y = sqrt(3)/6)
edge_points = function(p1, p2, n, margin = 0.06) {
  t = seq(margin, 1 - margin, length.out = n)
  data.frame(x = as.numeric(p1['x']) + t*(as.numeric(p2['x']) - as.numeric(p1['x'])),
             y = as.numeric(p1['y']) + t*(as.numeric(p2['y']) - as.numeric(p1['y'])))
}
bottom_pts = edge_points(apex_left, apex_right, n)   # measured
left_pts = edge_points(apex_left, apex_top, n)       # breast
right_pts = edge_points(apex_right, apex_top, n)     # formula

node_group = meta$group[match(sample_levels, meta$ID)]

node_df = rbind(data.frame(bottom_pts, edge = 'Measured', sample = sample_levels, group = node_group),
                data.frame(left_pts,   edge = 'Breast',   sample = sample_levels, group = node_group),
                data.frame(right_pts,  edge = 'Formula',  sample = sample_levels, group = node_group))
node_df$group = factor(node_df$group, levels = c('S-Preterm','N-Preterm','Term'))
node_df$edge  = factor(node_df$edge,  levels = c('Measured','Breast','Formula'))

# Sample IDs: vertical along the bottom edge, horizontal outside the two slanted edges
label_df = rbind(data.frame(bottom_pts, sample = sample_levels, nx = 0, ny = -1, angle = 90, hj = 1, off = 0.0126),
                 data.frame(left_pts, sample = sample_levels, nx = -1, ny = 0, angle = 0, hj = 1, off = 0.011),
                 data.frame(right_pts, sample = sample_levels, nx = 1, ny = 0, angle = 0, hj = 0, off = 0.011))
label_df$xlab = label_df$x + label_df$off * label_df$nx
label_df$ylab = label_df$y + label_df$off * label_df$ny

# Quadratic Bezier curves (start - control point pulled towards the centroid - end)
make_bezier = function(p_from, p_to, r_val, id, ctrl_pull = 0.5) {
  fx = as.numeric(p_from$x); fy = as.numeric(p_from$y)
  tx = as.numeric(p_to$x); ty = as.numeric(p_to$y)
  ctrlx = (fx+tx)/2 + ctrl_pull*(as.numeric(tri_center['x']) - (fx+tx)/2)
  ctrly = (fy+ty)/2 + ctrl_pull*(as.numeric(tri_center['y']) - (fy+ty)/2)
  data.frame(x = c(fx, ctrlx, tx), y = c(fy, ctrly, ty), group = id, R = r_val)
}
bezier_df = do.call(rbind, lapply(1:n, function(i) {
  rbind(make_bezier(bottom_pts[i,], left_pts[i,], cor_breast[i], paste0('B_', i)),
        make_bezier(bottom_pts[i,], right_pts[i,], cor_formula[i], paste0('F_', i)))
}))

tri_chord_plot = ggplot()+
  geom_bezier(data = bezier_df, aes(x = x, y = y, group = group, color = R),
              linewidth = 0.6, alpha = 0.75)+
  scale_color_gradient(low = "grey80", high = "deepskyblue1",
                       limits = range(c(cor_breast, cor_formula)), name = "Spearman R")+
  ggnewscale::new_scale_color()+
  geom_point(data = node_df, aes(x = x, y = y, shape = edge, color = group), size = 2.5)+
  scale_shape_manual(values = c("Measured" = 16, "Breast" = 17, "Formula" = 15),
                     labels = c("Measured" = "Measured",
                                "Breast"   = "Simulated, breast milk",
                                "Formula"  = "Simulated, formula"),
                     name = "Sample type")+
  scale_color_manual(values = c("S-Preterm" = "red", "N-Preterm" = "orange", "Term" = "blue"),
                     name = "Group")+
  guides(shape = guide_legend(override.aes = list(color = "black")),
         color = guide_legend(override.aes = list(shape = 16)))+
  geom_text(data = label_df, aes(x = xlab, y = ylab, label = sample, angle = angle, hjust = hj),
            size = 2.4, color = 'grey20')+
  coord_fixed(clip = 'off')+
  theme_void()+
  theme(legend.position = 'right', plot.margin = margin(0,5,0,5))+
  labs(title = NULL)#labs(title =NULL#labs(title = "Measured vs Simulated Community Correlation")


# ---------------------------------------------------------------------------------------------------------------- #
# Alpha diversity of the measured microbiome composition and of the simulated compositions under the two feeding
# conditions, compared between the three subject groups (S-Preterm, N-Preterm and Term). Pairwise Wilcoxon rank-sum
# tests were performed within each diversity index, with Benjamini-Hochberg correction applied across all pairwise
# comparisons of a given dataset.
# ---------------------------------------------------------------------------------------------------------------- #


alpha_plot = function(alpha_long, data_type, dodge_width = 0.75, color_alpha = 1) {
  stat_res = alpha_long %>%
    group_by(Index) %>%
    wilcox_test(Value ~ group, p.adjust.method = 'none') %>%
    ungroup() %>%
    adjust_pvalue(method = 'BH') %>%
    add_significance('p.adj') %>%
    add_xy_position(x = 'Index', dodge = dodge_width, step.increase = 0.12)
  
  ggplot(alpha_long, aes(x = Index, y = Value, fill = group))+
    geom_boxplot(width = 0.6, alpha = color_alpha, outlier.size = 1,
                 position = position_dodge(width = dodge_width))+
    scale_fill_manual(values = group_colors, name = "Group")+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    stat_pvalue_manual(stat_res, label = 'p.adj.signif',
                       tip.length = 0.01, size = 5, bracket.size = 0.4)+
    theme_bw()+
    labs(title = paste(data_type, " alpha diversity", sep = ''), y = "Value", x = "")+
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          title = element_text(size = 14),
          axis.title.y = element_text(size = 12))
}

alpha_stat_table = function(alpha_long, data_type) {
  alpha_long %>%
    group_by(Index) %>%
    wilcox_test(Value ~ group, p.adjust.method = 'none') %>%
    ungroup() %>%
    adjust_pvalue(method = 'BH') %>%
    add_significance('p.adj') %>%
    mutate(Dataset = data_type)
}
ori_alpha = alpha_diversity(selected_all_data[common_asv, seq_levels], sample_levels)
breast_alpha = alpha_diversity(breast_sim_data[common_asv, sample_levels], sample_levels)
formula_alpha = alpha_diversity(formula_sim_data[common_asv, sample_levels], sample_levels)

measured_alpha_plot = alpha_plot(ori_alpha, 'Measured')+labs(title = NULL)
breast_alpha_plot = alpha_plot(breast_alpha, 'Simulated (breastfed)')+labs(title = NULL)
formula_alpha_plot = alpha_plot(formula_alpha, 'Simulated (formulafed)')+labs(title = NULL)

alpha_stat_df = rbind(alpha_stat_table(ori_alpha, 'Measured'),
                      alpha_stat_table(breast_alpha, 'Breast'),
                      alpha_stat_table(formula_alpha, 'Formula'))



# ---------------------------------------------------------------------------------------------------------------- #
# To evaluate the realism of the simulation results, the subject group in which the highest mean relative abundance
# occurs was determined for each genus, and compared between the measured microbiome and the microbiomes simulated
# under the two feeding conditions. Genera showing an inconsistent pattern across the three datasets are highlighted
# with a grey background. Asterisks denote genera that differ significantly between the three subject groups within
# that dataset (Kruskal-Wallis test, BH-corrected, p.adj < 0.05).
# ---------------------------------------------------------------------------------------------------------------- #
higher_in_group = function(genus_data, metadata, data_type) {
  group_levels = c('S-Preterm', 'N-Preterm', 'Term')
  sample_cols = setdiff(colnames(genus_data), 'Genus')
  group_means = sapply(group_levels, function(g) {
    ids = metadata$ID[metadata$group == g]
    rowMeans(genus_data[, ids[ids %in% sample_cols], drop = FALSE])
  })
  data.frame(Genus = genus_data$Genus, group_means, check.names = FALSE) %>%
    mutate(higher_in = ifelse(rowSums(group_means) == 0, 'Absent',
                              group_levels[max.col(group_means, ties.method = 'first')]),
           Status = data_type)
}

genus_kw_test = function(genus_data, metadata, data_type) {
  sample_cols = setdiff(colnames(genus_data), 'Genus')
  grp = factor(metadata$group[match(sample_cols, metadata$ID)])
  pvals = apply(genus_data[, sample_cols], 1, function(v) {
    if (length(unique(as.numeric(v))) < 2) return(NA_real_)
    kruskal.test(as.numeric(v) ~ grp)$p.value
  })
  data.frame(Genus = genus_data$Genus, p = pvals,
             p.adj = p.adjust(pvals, method = 'BH'), Status = data_type)
}

measured_higher = higher_in_group(measured_genus, meta, 'Measured')
breast_higher = higher_in_group(breast_genus, meta, 'Breast')
formula_higher = higher_in_group(formula_genus, meta, 'Formula')

abun_comp_df = rbind(measured_higher, breast_higher, formula_higher) %>%
  mutate(Status = factor(Status, levels = c('Formula', 'Breast', 'Measured')))

species_consistency = abun_comp_df %>%
  filter(higher_in != 'Absent') %>%
  group_by(Genus) %>%
  summarise(n_groups = n_distinct(higher_in), .groups = 'drop') %>%
  mutate(background = ifelse(n_groups > 1, 'darkgray', 'white'))

plot_data = expand.grid(Genus = unique(abun_comp_df$Genus),
                        Status = levels(abun_comp_df$Status)) %>%
  left_join(abun_comp_df, by = c('Genus', 'Status')) %>%
  left_join(species_consistency, by = 'Genus') %>%
  mutate(background = ifelse(is.na(background), 'white', background),
         Status = factor(Status, levels = c('Formula', 'Breast', 'Measured')))

genus_kw_df = rbind(genus_kw_test(measured_genus, meta, 'Measured'),
                    genus_kw_test(breast_genus,   meta, 'Breast'),
                    genus_kw_test(formula_genus,  meta, 'Formula'))

star_df = genus_kw_df %>%
  mutate(star = ifelse(!is.na(p.adj) & p.adj < 0.05, '*', ''),
         Status = factor(Status, levels = c('Formula', 'Breast', 'Measured'))) %>%
  filter(star != '', Genus %in% plot_data$Genus)

status_labels = c("Measured"  = "Measured",
                  "Breast"    = "Simulated, breast milk",
                  "Formula"   = "Simulated, formula")

abun_comp_plot = ggplot(plot_data, aes(x = Genus, y = Status))+
  geom_tile(aes(fill = background), color = 'black', width = 1, height = 1)+
  scale_fill_identity()+
  geom_point(data = subset(plot_data, !is.na(higher_in)),
             aes(color = higher_in), size = 4.5)+
  geom_text(data = star_df, aes(x = Genus, y = Status, label = star),
            inherit.aes = FALSE, size = 6, color = 'white',
            fontface = 'bold', vjust = 0.72)+
  scale_color_manual(values = c("S-Preterm" = "red", "N-Preterm" = "orange",
                                "Term" = "blue", "Absent" = "black"),
                     labels = c("S-Preterm" = 'Highest in S-Preterm',
                                "N-Preterm" = 'Highest in N-Preterm',
                                "Term" = 'Highest in Term',
                                "Absent" = 'Absent in all groups'))+
  scale_y_discrete(labels = status_labels)+
  coord_fixed()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic', size = 8),
        axis.text.y = element_text(face = 'bold', size = 8),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title = element_blank())


# ---------------------------------------------------------------------------------------------------------------- #
# Overlap of the genome-scale metabolic models used to construct the in silico communities of the three subject
# groups. A model is counted as used in a group when the corresponding ASV is present in at least one member of
# that group. All three groups draw on the same library of 191 models and the same ASV-to-GEM mapping, so the
# comparison reflects differences in community membership rather than differences in model availability.
# ---------------------------------------------------------------------------------------------------------------- #

gem_used_by_group = lapply(group_levels, function(g) {
  seq_ids = rownames(meta)[meta$group == g]
  sub_data = all_amount_data[, seq_ids, drop = FALSE]
  asv_present = rownames(sub_data)[rowSums(sub_data) > 0]
  ASV_GEM_map$strain[ASV_GEM_map$ASV_ID %in% asv_present]
})
names(gem_used_by_group) = c("S-Preterm","N-Preterm","Term")

GEMs_venn_plot = ggvenn(
  data = gem_used_by_group,
  fill_color = c("red","orange","blue"),
  stroke_size = 0.8,
  set_name_size = 4,
  text_size = 4.5)

# Number of models used per group and per subject
gem_usage_summary = data.frame(
  group = factor(group_levels, levels = group_levels),
  n_subjects = sapply(group_levels, function(g) sum(meta$group == g)),
  n_models_used = sapply(gem_used_by_group, length),
  row.names = NULL)

gem_per_subject = data.frame(
  ID = meta$ID,
  group = factor(meta$group, levels = group_levels),
  n_models = colSums(all_amount_data[, rownames(meta), drop = FALSE] > 0),
  row.names = NULL)

gem_per_subject_stat = gem_per_subject %>%
  rstatix::wilcox_test(n_models ~ group, p.adjust.method = 'none') %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  rstatix::add_xy_position(x = 'group', step.increase = 0.10)

gem_per_subject_plot = ggplot(gem_per_subject, aes(x = group, y = n_models, fill = group))+
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1)+
  scale_fill_manual(values = group_colors, name = "Group")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  stat_pvalue_manual(gem_per_subject_stat %>% filter(p.adj < 0.05),
                     label = 'p.adj.signif', tip.length = 0.01, size = 5)+
  theme_bw()+
  labs(title = "Number of metabolic models per individual community",
       x = "", y = "Models used")+
  theme(axis.text = element_text(size = 12), title = element_text(size = 12),
        legend.position = 'none')

# ================================================================================================================== #
# Part 3  Blood-culture isolates and the gut abundance of the same taxa
#
#   Manuscript: Results 1, opening paragraph
#   Moved here from the end of the original script. Depends only on asv_rel and
#   asv_taxonomy from Part 1, so it runs at this position unchanged.
# ================================================================================================================== #

# ---- Blood-culture isolates and the gut abundance of the same taxa -----------------------------------------------
# Each septic infant is compared with its matched control for the abundance of the taxon group its own pathogen
# belongs to. Klebsiella is absent from the SILVA annotation of this dataset, which does not resolve it from
# Escherichia-Shigella in the V4 region, so the comparison is made at the level of the Enterobacteriaceae.
preterm_meta = read.table(file.path(data_dir, "metadata/all_select_meta.txt"),
                          header = TRUE, sep = "\t", quote = '"', check.names = FALSE)

entero_genera_silva = c("Escherichia-Shigella","Enterobacter","Citrobacter",
                        "Serratia","Yersinia","Pantoea")

# Relative abundance of a genus set, keyed on the sequencing IDs used by the abundance table
genus_rel_by_seqid = function(genera) {
  hash = asv_taxonomy$asv_hash[asv_taxonomy$genus %in% genera]
  100 * colSums(asv_rel[rownames(asv_rel) %in% hash, , drop = FALSE])
}

rel_staph  = genus_rel_by_seqid("Staphylococcus")
rel_entco  = genus_rel_by_seqid("Enterococcus")
rel_entero = genus_rel_by_seqid(entero_genera_silva)

path_pairs = preterm_meta %>%
  mutate(seqid = sub("_S\\d+_L001$", "", rownames(preterm_meta))) %>%
  select(seqid, Group, matched, Pathogen, DoL, Age) %>%
  mutate(path_group = case_when(
    Group == 'Control'                 ~ NA_character_,
    grepl("Staphylococcus", Pathogen)  ~ "Staphylococcus",
    grepl("Enterococcus", Pathogen)    ~ "Enterococcus",
    TRUE                               ~ "Enterobacteriaceae")) %>%
  group_by(matched) %>%
  mutate(path_group = path_group[Group == 'Sepsis']) %>%
  ungroup() %>%
  mutate(rel = case_when(
    path_group == "Staphylococcus"     ~ rel_staph[seqid],
    path_group == "Enterococcus"       ~ rel_entco[seqid],
    path_group == "Enterobacteriaceae" ~ rel_entero[seqid]))

path_wide = path_pairs %>%
  select(matched, Group, path_group, rel) %>%
  pivot_wider(names_from = Group, values_from = rel) %>%
  mutate(diff = Sepsis - Control) %>%
  arrange(matched)

# Days between stool sampling and the sepsis diagnosis, to establish that the samples precede the event
path_timing = preterm_meta %>%
  filter(Group == 'Sepsis') %>%
  summarise(n = n(),
            median_days = median(DoL - Age),
            min_days = min(DoL - Age),
            max_days = max(DoL - Age))

path_stat_overall = wilcox.test(path_wide$Sepsis, path_wide$Control, paired = TRUE)

path_stat_bygroup = path_wide %>%
  group_by(path_group) %>%
  summarise(n = n(),
            n_higher = sum(diff > 0),
            median_sepsis = median(Sepsis),
            median_control = median(Control),
            p = wilcox.test(Sepsis, Control, paired = TRUE)$p.value, .groups = 'drop')

path_pair_plot = path_pairs %>%
  filter(!is.na(rel)) %>%
  mutate(Group = factor(Group, levels = c('Control','Sepsis'))) %>%
  ggplot(aes(x = Group, y = rel + 0.01, group = matched))+
  geom_line(colour = 'grey60', linewidth = 0.4)+
  geom_point(aes(colour = path_group), size = 2)+
  scale_y_log10()+
  scale_colour_brewer(palette = 'Dark2', name = "Pathogen group")+
  theme_bw()+
  labs(x = "", y = "Relative abundance of the pathogen's taxon group (%)")+
  theme(axis.text = element_text(size = 10), legend.position = 'bottom')

supp_table_pathogen = path_wide %>%
  left_join(preterm_meta %>% filter(Group == 'Sepsis') %>% select(matched, Pathogen),
            by = 'matched') %>%
  mutate(across(c(Sepsis, Control, diff), ~ signif(.x, 3))) %>%
  select(Pair = matched, Pathogen, `Taxon group` = path_group,
         `Sepsis infant (%)` = Sepsis, `Matched control (%)` = Control,
         Difference = diff)

# ================================================================================================================== #
# Part 4  Baseline fermentation product production
#
#   Manuscript: Results 2; Figure 2; Tables S8 to S10
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Predicted SCFA concentrations at the end of the simulation, per individual and per feeding condition. Several of
# the acids are zero in most individuals, so the number of subjects with a non-zero prediction is recorded alongside
# the median: a paired or grouped test on such an acid rests on far fewer observations than the nominal sample size
# and its p value is not treated as evidence.
# ---------------------------------------------------------------------------------------------------------------- #
breast_scfa = extract_scfa(breast_ext, meta, 'Breast')
formula_scfa = extract_scfa(formula_ext, meta, 'Formula')

all_scfa_conc = rbind(breast_scfa, formula_scfa) %>%
  mutate(Diet = factor(Diet, levels = c('Breast','Formula')),
         sub = factor(sub, levels = scfa_levels))
all_scfa_conc_wider = all_scfa_conc %>% pivot_wider(names_from = 'sub', values_from = 'value')
all_scfa_conc_wider[is.na(all_scfa_conc_wider)] = 0

# How much of each acid is actually informative
scfa_detection = all_scfa_conc %>%
  group_by(Diet, sub, group) %>%
  summarise(n = n(),
            n_nonzero = sum(value > 0),
            pct_nonzero = 100 * n_nonzero / n,
            median = median(value),
            iqr = IQR(value), .groups = 'drop')

scfa_bygroup = all_scfa_conc %>%
  group_by(Diet, sub, group) %>%
  summarise(median = median(value), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = median)

# ---------------------------------------------------------------------------------------------------------------- #
# Per-individual profiles
# ---------------------------------------------------------------------------------------------------------------- #
scfa_plot = function(scfa_df, data_type) {
  scfa_df = scfa_df %>% mutate(ID = factor(ID, levels = unique(ID[order(group)])))
  ggplot(scfa_df, aes(x = value, y = ID, fill = group))+
    geom_bar(stat = 'identity')+
    # Nesting the panels inside the compound class keeps the distinction between the true short-chain
    # fatty acids and the related organic acids visible on the figure itself
    ggh4x::facet_nested(. ~ Class + sub, scales = 'free_x', space = 'fixed')+
    scale_fill_manual(values = group_colors, name = "Group")+
    theme_bw()+
    labs(title = paste(data_type, " predicted fermentation products", sep = ''),
         x = "Fermentation product, log(mM+1)", y = "")+
    theme(strip.text.x = element_text(size = 8),
          strip.text.y = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_markdown(size = 7),
          axis.text.x = element_text(angle = 60, vjust = 0.6, size = 8),
          panel.spacing = unit(0.2, "lines"),
          panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6))
}

breast_scfa_plot = scfa_plot(breast_scfa, 'Simulated (breastfed)')
formula_scfa_plot = scfa_plot(formula_scfa, 'Simulated (formulafed)')

comb_scfa_plot = breast_scfa_plot + formula_scfa_plot +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A',tag_prefix = "(",tag_suffix = ")")
ggsave(plot = comb_scfa_plot,filename = 'Figures/comb_scfa_plot.jpeg',
       dpi = 1600,height = 20,width = 36,units = 'cm')

scfa_group_median = all_scfa_conc %>%
  group_by(Diet, sub, group) %>%
  summarise(median = signif(median(value), 4), .groups = 'drop')


# ---------------------------------------------------------------------------------------------------------------- #
# Differences between the three subject groups, tested separately within each feeding condition with pairwise
# Wilcoxon tests and Benjamini-Hochberg correction across all pairwise comparisons of that condition.
# ---------------------------------------------------------------------------------------------------------------- #
scfa_group_stat = function(scfa_df, diet_type) {
  scfa_df %>%
    group_by(sub) %>%
    rstatix::wilcox_test(value ~ group, p.adjust.method = 'none') %>%
    ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance('p.adj') %>%
    mutate(Diet = diet_type)
}
breast_scfa_stat = scfa_group_stat(breast_scfa, 'Breast')
formula_scfa_stat = scfa_group_stat(formula_scfa, 'Formula')
scfa_group_stat_all = rbind(breast_scfa_stat, formula_scfa_stat)

breast_scfa_box = scfa_group_box(breast_scfa, 'Simulated (breastfed)')
formula_scfa_box = scfa_group_box(formula_scfa, 'Simulated (formulafed)')

fig_scfa_group = breast_scfa_box/formula_scfa_box+
  plot_layout(guides = 'collect')&theme(legend.position = 'bottom')
ggsave(plot = fig_scfa_group,filename = 'Figures/fig_scfa_group.jpeg',
       dpi = 1600,width = 26,height = 26,units = 'cm')

supp_table_S8 = rbind(breast_scfa_stat, formula_scfa_stat) %>%
  select(Diet, Metabolite = sub, Group1 = group1, Group2 = group2,
         n1, n2, statistic, p, p.adj, p.adj.signif) %>%
  left_join(scfa_group_median %>% rename(Group1 = group, Median1 = median),
            by = c('Diet', 'Metabolite' = 'sub', 'Group1')) %>%
  left_join(scfa_group_median %>% rename(Group2 = group, Median2 = median),
            by = c('Diet', 'Metabolite' = 'sub', 'Group2')) %>%
  select(Diet, Metabolite,
         Group1, n1, Median1,
         Group2, n2, Median2,
         statistic, p, p.adj, p.adj.signif) %>%
  arrange(Diet, Metabolite)

write.table(supp_table_S8, "Supplementary_Table_S8.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------------------------------------------- #
# Effect of the feeding condition, tested within subject. The comparison is run on the whole cohort and again on the
# preterm infants alone, since the two feeding conditions differ mainly in their oligosaccharide content and the
# preterm communities are the ones that largely cannot use it.
# ---------------------------------------------------------------------------------------------------------------- #
scfa_diet_stat = all_scfa_conc %>%
  arrange(sub, Diet, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(value ~ Diet, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

scfa_diet_stat_preterm = all_scfa_conc %>%
  filter(group %in% preterm_levels) %>%
  arrange(sub, Diet, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(value ~ Diet, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

# The effective sample size of a paired test excludes pairs with no difference, which for the sparse acids leaves
# only a handful of informative pairs; this is reported so that a small p value on such an acid is not over-read
scfa_diet_paired = all_scfa_conc %>%
  select(ID, group, Diet, sub, value) %>%
  pivot_wider(names_from = Diet, values_from = value) %>%
  mutate(Diff = Breast - Formula)

scfa_diet_summary = scfa_diet_paired %>%
  group_by(sub) %>%
  summarise(n = n(),
            n_tied = sum(Diff == 0),
            n_effective = n - n_tied,
            n_breast_higher = sum(Diff > 0),
            median_breast = median(Breast),
            median_formula = median(Formula),
            median_diff = median(Diff), .groups = 'drop') %>%
  left_join(scfa_diet_stat %>% select(sub, p, p.adj, p.adj.signif), by = 'sub')

# ---------------------------------------------------------------------------------------------------------------- #
# Whether the spread of the predictions differs between the two feeding conditions. The Brown-Forsythe form of
# Levene's test is used because the concentrations are strongly right-skewed and carry many zeros, which makes the
# variance ratio test unusable.
# ---------------------------------------------------------------------------------------------------------------- #
scfa_spread_stat = all_scfa_conc %>%
  group_by(sub, group) %>%
  group_modify(~ {
    lt = car::leveneTest(value ~ Diet, data = .x, center = median)
    tibble(p_levene = lt$`Pr(>F)`[1],
           iqr_breast = IQR(.x$value[.x$Diet == 'Breast']),
           iqr_formula = IQR(.x$value[.x$Diet == 'Formula']))
  }) %>%
  ungroup() %>%
  rstatix::adjust_pvalue(p.col = 'p_levene', method = 'BH') %>%
  arrange(p_levene)


# ---------------------------------------------------------------------------------------------------------------- #
# Which taxa produce each acid. Net flux per organism is summed over the simulation, so that a taxon appearing
# transiently on the positive side is not counted as a producer; only net secretors are retained. This turns the
# compositional differences between the groups into a statement about which pathways are running.
# ---------------------------------------------------------------------------------------------------------------- #
scfa_producers = function(extracted, metadata, diet_type) {
  extracted$flux %>%
    filter(rea %in% scfas) %>%
    group_by(ID, spec, rea) %>%
    summarise(net = sum(mflux, na.rm = TRUE), .groups = 'drop') %>%
    filter(net > 0) %>%
    left_join(metadata %>% select(ID, group), by = 'ID') %>%
    mutate(group = factor(group, levels = group_levels),
           sub = ifelse(rea %in% names(scfa_flux_labels), scfa_flux_labels[rea], rea),
           Diet = diet_type)
}

producer_long = rbind(scfa_producers(breast_ext, meta, 'Breast'),
                      scfa_producers(formula_ext, meta, 'Formula')) %>%
  mutate(Diet = factor(Diet, levels = c('Breast','Formula')),
         sub = factor(sub, levels = scfa_levels))



# ---------------------------------------------------------------------------------------------------------------- #
# The five largest contributors to each fermentation product in each group, ranked by the mean net flux across the
# subjects whose community contains the taxon. Ranking by the mean rather than by the group total keeps the
# comparison independent of group size, which differs between the preterm groups and the term group.
# ---------------------------------------------------------------------------------------------------------------- #
producer_rank = producer_long %>%
  group_by(Diet, sub, group, spec) %>%
  summarise(n_subjects = n_distinct(ID),
            mean_net = mean(net),
            total_net = sum(net), .groups = 'drop') %>%
  arrange(Diet, sub, group, desc(mean_net))

producer_top5 = producer_rank %>%
  group_by(Diet, sub, group) %>%
  slice_max(mean_net, n = 5) %>%
  ungroup()

# ---------------------------------------------------------------------------------------------------------------- #
# Composition of the producer pool, summarised at the level of the taxonomic groups that carry the pathways of
# interest. The individual strains that make up the five largest contributors change with the nutrient environment,
# so the shares are computed for both media and reported side by side; the strain-level detail is in the
# supplementary figures. Taxa are grouped by their role in the pathways examined here rather than by rank.
# ---------------------------------------------------------------------------------------------------------------- #
producer_class_data = producer_top5 %>%
  mutate(sub = factor(sub, levels = scfa_levels),
         Class = factor(scfa_class[as.character(sub)], levels = scfa_class_levels),
         taxon_group = factor(case_when(
           grepl("Bacteroides|Parabacteroides|Prevotella", spec) ~ "Bacteroidetes",
           grepl("Bifidobacterium", spec)                        ~ "Bifidobacterium",
           grepl("Escherichia|Enterobacter|Citrobacter|Yersinia|Klebsiella|Serratia|Pantoea", spec)
           ~ "Enterobacteriaceae",
           TRUE ~ "Other"),
           levels = c("Enterobacteriaceae","Bacteroidetes","Bifidobacterium","Other"))) %>%
  group_by(Diet, sub, Class, group, taxon_group) %>%
  summarise(total = sum(mean_net), .groups = 'drop') %>%
  group_by(Diet, sub, group) %>%
  mutate(pct = 100 * total / sum(total)) %>%
  ungroup()

# The same shares in wide form, for the values quoted in the text
producer_class_wide = producer_class_data %>%
  select(-total, -Class) %>%
  mutate(pct = round(pct)) %>%
  pivot_wider(names_from = taxon_group, values_from = pct, values_fill = 0)

class_colors = c("Enterobacteriaceae" = "#ff7f00",
                 "Bacteroidetes" = "#1f78b4",
                 "Bifidobacterium" = "#e31a1c",
                 "Other" = "grey70")

producer_flux_plot = ggplot(producer_class_data,
                            aes(x = group, y = pct, fill = taxon_group))+
  geom_bar(stat = 'identity', width = 0.7)+
  # Nesting the panels inside the compound class keeps the distinction between the true short-chain fatty acids
  # and the related organic acids visible on the figure itself, as in the concentration panels
  ggh4x::facet_nested(Diet ~ Class + sub)+
  scale_fill_manual(values = class_colors, name = "")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  labs(x = "", y = "Share of the five largest contributors (%)")+
  theme(axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14,face = 'italic'),
        legend.position = 'bottom')
ggsave(plot = producer_flux_plot,filename = 'Figures/producer_flux_plot.jpeg',
       dpi = 1200,width = 32,height = 16,units = 'cm')


# ---------------------------------------------------------------------------------------------------------------- #
entero_bracket = entero_stat %>%
  filter(p.adj < 0.05) %>%
  rstatix::add_xy_position(x = 'group', step.increase = 0.08)

entero_plot = ggplot(taxa_features, aes(x = group, y = entero_pct, fill = group))+
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 1, alpha = 0.8, colour = 'grey25')+
  scale_fill_manual(values = group_colors, name = "Group")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  stat_pvalue_manual(entero_bracket, label = 'p.adj.signif',
                     tip.length = 0.01, size = 8)+
  theme_bw()+
  labs(x = "",y = expression(italic(Enterobacteriaceae) ~ "relative abundance (%)"))+
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size = 12),legend.position = 'none')
ggsave(plot = entero_plot,filename = 'Figures/entero_plot.jpeg',
       dpi = 1200,height = 10,width = 10,units = 'cm')


bif_rel_stat = taxa_features %>%
  rstatix::wilcox_test(bif_rel ~ group) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

bif_rel_bracket = bif_rel_stat %>%
  filter(p.adj < 0.05) %>%
  rstatix::add_xy_position(x = 'group', step.increase = 0.08)

bif_rel_plot = ggplot(taxa_features, aes(x = group, y = bif_rel, fill = group))+
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 1, alpha = 0.8, colour = 'grey25')+
  scale_fill_manual(values = group_colors, name = "Group")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  stat_pvalue_manual(bif_rel_bracket, label = 'p.adj.signif',
                     tip.length = 0.01, size = 8)+
  theme_bw()+
  labs(x = "",y = expression(italic(Bifidobacterium) ~ "relative abundance (%)"))+
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size = 12),legend.position = 'none')
ggsave(plot = bif_rel_plot,filename = 'Figures/bif_ref_plot.jpeg',
       dpi = 1200,height = 10,width = 10,units = 'cm')



taxa_features %>%
  group_by(group) %>%
  summarise(n = n(),
            median_rel = round(median(bac_rel), 3),
            q25 = round(quantile(bac_rel, .25), 3),
            q75 = round(quantile(bac_rel, .75), 3),
            n_zero = sum(bac_rel == 0), .groups = 'drop')

bac_rel_stat = taxa_features %>%
  rstatix::wilcox_test(bac_rel ~ group) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')
bac_rel_stat %>% select(group1, group2, statistic, p, p.adj, p.adj.signif)

taxa_features %>% rstatix::kruskal_test(bac_rel ~ group)

bac_rel_bracket = bac_rel_stat %>%
  filter(p.adj < 0.05) %>%
  rstatix::add_xy_position(x = 'group', step.increase = 0.08)

bac_rel_plot = ggplot(taxa_features, aes(x = group, y = bac_rel, fill = group))+
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 1, alpha = 0.8, colour = 'grey25')+
  scale_fill_manual(values = group_colors, name = "Group")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  stat_pvalue_manual(bac_rel_bracket, label = 'p.adj.signif',
                     tip.length = 0.01, size = 8)+
  theme_bw()+
  labs(x = "",y = expression(italic(Bacteroidetes) ~ "relative abundance (%)"))+
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size = 12),legend.position = 'none')
ggsave(plot = bac_rel_plot,filename = 'Figures/bac_rel_plot.jpeg',
       dpi = 1200,width = 10,height = 10,units = 'cm')


# ---------------------------------------------------------------------------------------------------------------- #
# Net secretion of each fermentation product by taxon, summarised within a subject group. A taxon's net flux is the
# sum of its exchange flux over the whole simulation, averaged across the subjects whose community contains it;
# only taxa with a positive net flux above the display threshold are shown, so each panel lists the producers of
# that product in that group.
# ---------------------------------------------------------------------------------------------------------------- #
flux_display_threshold = 1e4

net_flux_by_taxon = function(extracted, metadata, diet_type) {
  extracted$flux %>%
    group_by(ID, spec, rea) %>%
    summarise(net = sum(mflux, na.rm = TRUE), .groups = 'drop') %>%
    left_join(metadata %>% select(ID, group), by = 'ID') %>%
    mutate(group = factor(group, levels = group_levels),
           sub = factor(ifelse(rea %in% names(scfa_flux_labels),
                               scfa_flux_labels[rea], rea), levels = scfa_levels),
           Diet = diet_type)
}

producer_flux_panels = function(net_flux_df, target_group, threshold = flux_display_threshold) {
  d = net_flux_df %>%
    filter(group == target_group) %>%
    group_by(sub, spec) %>%
    summarise(mean_net = mean(net, na.rm = TRUE),
              n_subjects = n_distinct(ID), .groups = 'drop') %>%
    filter(mean_net > threshold) %>%
    mutate(taxon = spec %>%
             gsub("_(sp|subsp)_", " sp. ", .) %>%
             gsub("_", " ", .))
  
  lapply(levels(d$sub), function(m) {
    dm = d %>% filter(sub == m)
    if (nrow(dm) == 0) return(NULL)
    ggplot(dm, aes(x = reorder(taxon, mean_net), y = mean_net))+
      geom_bar(stat = 'identity', width = 0.7, fill = '#1b9e77')+
      scale_y_continuous(labels = scales::scientific,
                         expand = expansion(mult = c(0, 0.05)))+
      coord_flip()+
      labs(title = paste0(m, " production"), x = "", y = "Mean net flux")+
      theme_minimal()+
      theme(axis.text.y = element_text(size = 7, face = 'italic'),
            axis.text.x = element_text(size = 7),
            plot.title = element_text(size = 12),
            panel.grid.major.y = element_blank())
  }) %>% setNames(levels(d$sub)) %>% Filter(Negate(is.null), .)
}

breast_net_flux = net_flux_by_taxon(breast_ext, meta, 'Breast')
formula_net_flux = net_flux_by_taxon(formula_ext, meta, 'Formula')

sp_flux_panels = producer_flux_panels(breast_net_flux, 'S-Preterm')
np_flux_panels = producer_flux_panels(breast_net_flux, 'N-Preterm')
tm_flux_panels = producer_flux_panels(breast_net_flux, 'Term')

sp_flux_panels_f = producer_flux_panels(formula_net_flux, 'S-Preterm')
np_flux_panels_f = producer_flux_panels(formula_net_flux, 'N-Preterm')
tm_flux_panels_f = producer_flux_panels(formula_net_flux, 'Term')

acetate_fluxes = (sp_flux_panels$Acetate | np_flux_panels$Acetate | tm_flux_panels$Acetate) /
  (sp_flux_panels_f$Acetate | np_flux_panels_f$Acetate | tm_flux_panels_f$Acetate)+
  plot_annotation(tag_levels = 'A',tag_prefix = '(',tag_suffix = ')')
ggsave(plot = acetate_fluxes,filename = 'Figures/acetate_fluxes.jpeg',dpi = 1600,
       height = 30,width = 36,units = 'cm')

propionate_fluxes = (sp_flux_panels$Propionate | np_flux_panels$Propionate | tm_flux_panels$Propionate) /
  (sp_flux_panels_f$Propionate | np_flux_panels_f$Propionate | tm_flux_panels_f$Propionate)+
  plot_annotation(tag_levels = 'A',tag_prefix = '(',tag_suffix = ')')
ggsave(plot = propionate_fluxes,filename = 'Figures/propionate_fluxes.jpeg',dpi = 1600,
       height = 30,width = 36, units = 'cm')

butyrate_fluxes = (sp_flux_panels$Butyrate | np_flux_panels$Butyrate | tm_flux_panels$Butyrate) /
  (sp_flux_panels_f$Butyrate | np_flux_panels_f$Butyrate | tm_flux_panels_f$Butyrate)+
  plot_annotation(tag_levels = 'A',tag_prefix = '(',tag_suffix = ')')
ggsave(plot = butyrate_fluxes,filename = 'Figures/butyrate_fluxes.jpeg',dpi = 1600,
       height = 10,width = 36, units = 'cm')

formate_fluxes = (sp_flux_panels$Formate | np_flux_panels$Formate | tm_flux_panels$Formate) /
  (sp_flux_panels_f$Formate | np_flux_panels_f$Formate | tm_flux_panels_f$Formate)+
  plot_annotation(tag_levels = 'A',tag_prefix = '(',tag_suffix = ')')
ggsave(plot = formate_fluxes,filename = 'Figures/formate_fluxes.jpeg',dpi = 1600,
       height = 36,width = 36, units = 'cm')

succinate_fluxes = (sp_flux_panels$Succinate | np_flux_panels$Succinate | tm_flux_panels$Succinate) /
  (sp_flux_panels_f$Succinate | np_flux_panels_f$Succinate | tm_flux_panels_f$Succinate)+
  plot_annotation(tag_levels = 'A',tag_prefix = '(',tag_suffix = ')')
ggsave(plot = succinate_fluxes,filename = 'Figures/succinate_fluxes.jpeg',dpi = 1600,
       height = 30,width = 36, units = 'cm')

R_Lactate_fluxes = (sp_flux_panels$`(R)-Lactate` | np_flux_panels$`(R)-Lactate` | tm_flux_panels$`(R)-Lactate`) /
  (sp_flux_panels_f$`(R)-Lactate` | np_flux_panels_f$`(R)-Lactate` | tm_flux_panels_f$`(R)-Lactate`)+
  plot_annotation(tag_levels = 'A',tag_prefix = '(',tag_suffix = ')')
ggsave(plot = R_Lactate_fluxes,filename = 'Figures/R_Lactate_fluxes.jpeg',dpi = 1600,
       height = 16,width = 36, units = 'cm')

S_Lactate_fluxes = (sp_flux_panels$`(S)-Lactate` | np_flux_panels$`(S)-Lactate` | tm_flux_panels$`(S)-Lactate`) /
  (sp_flux_panels_f$`(S)-Lactate` | np_flux_panels_f$`(S)-Lactate` | tm_flux_panels_f$`(S)-Lactate`)+
  plot_annotation(tag_levels = 'A',tag_prefix = '(',tag_suffix = ')')
ggsave(plot = S_Lactate_fluxes,filename = 'Figures/S_Lactate_fluxes.jpeg',dpi = 1600,
       height = 16,width = 36, units = 'cm')

supp_table_S10 = producer_top5 %>%
  mutate(taxon_group = case_when(
    grepl("Bacteroides|Parabacteroides|Prevotella", spec) ~ "Bacteroidetes",
    grepl("Bifidobacterium", spec)                        ~ "Bifidobacterium",
    grepl("Escherichia|Enterobacter|Citrobacter|Yersinia|Klebsiella|Serratia|Pantoea", spec)
    ~ "Enterobacteriaceae",
    TRUE ~ "Other"),
    Taxon = spec %>% gsub("_(sp|subsp)_", " sp. ", .) %>% gsub("_", " ", .)) %>%
  group_by(Diet, sub, group) %>%
  mutate(taxon_share = 100 * mean_net / sum(mean_net)) %>%
  arrange(desc(mean_net), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  group_by(Diet, sub, group, taxon_group) %>%
  mutate(class_share = sum(taxon_share)) %>%
  ungroup() %>%
  mutate(`Taxon share (%)` = round(taxon_share, 1),
         `Class share (%)` = round(class_share, 1),
         `Mean net flux` = signif(mean_net, 4)) %>%
  select(Medium = Diet,
         Metabolite = sub,
         Group = group,
         Rank,
         Taxon,
         `Taxonomic group` = taxon_group,
         `Subjects with net secretion` = n_subjects,
         `Mean net flux`,
         `Taxon share (%)`,
         `Class share (%)`) %>%
  arrange(Medium, Metabolite, Group, Rank)

write.table(supp_table_S10, "Supplementary_Table_S10.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Paired comparison restricted to the preterm infants -----------------------------------------------------------
scfa_diet_summary_preterm = scfa_diet_paired %>%
  filter(group != 'Term') %>%
  group_by(sub) %>%
  summarise(n = n(),
            n_tied = sum(Diff == 0),
            n_effective = n - n_tied,
            n_breast_higher = sum(Diff > 0),
            median_breast = signif(median(Breast), 3),
            median_formula = signif(median(Formula), 3),
            median_diff = signif(median(Diff), 3), .groups = 'drop') %>%
  left_join(scfa_diet_stat_preterm %>% select(sub, p, p.adj, p.adj.signif),
            by = 'sub')

supp_table_diet_paired = rbind(
  scfa_diet_summary %>% mutate(Cohort = 'All infants'),
  scfa_diet_summary_preterm %>% mutate(Cohort = 'Preterm only')) %>%
  mutate(across(c(median_breast, median_formula, median_diff), ~ signif(.x, 3)),
         across(c(p, p.adj), ~ signif(.x, 3))) %>%
  select(Cohort, Metabolite = sub,
         Pairs = n, `Tied pairs` = n_tied, `Effective pairs` = n_effective,
         `Pairs higher under breast milk` = n_breast_higher,
         `Median, breast milk` = median_breast,
         `Median, formula` = median_formula,
         `Median difference` = median_diff,
         p, p.adj, p.adj.signif) %>%
  arrange(Cohort, Metabolite)

write.table(supp_table_diet_paired, "Supplementary_Table_S11.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ================================================================================================================== #
# Part 5  Nutritional conditions and HMO utilisation
#
#   Manuscript: Results 3; Figure 3A to 3E
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Side-by-side comparison of the two media. Substrates are arranged in three blocks: those present in both media,
# then those unique to each, so that the qualitative difference between the media is as visible as the quantitative
# one. Concentrations span ten orders of magnitude and include zeros, so the axis is logarithmic with the zeros
# drawn at a floor one decade below the smallest non-zero concentration; breast milk extends to the left of the
# central axis and formula to the right.
# ---------------------------------------------------------------------------------------------------------------- #

# Reaction ID to metabolite name, taken from the arena's own media vector
mediac_names = setNames(names(breast_ext$mediac), breast_ext$mediac)

# Names for the exchange reactions the arena's media vector does not label
extra_names = c("EX_phyQ(e)"   = "Phylloquinone (vitamin K1)",
                "EX_avite1(e)" = "alpha-Tocopherol (vitamin E)")
mediac_names = c(mediac_names, extra_names[!names(extra_names) %in% names(mediac_names)])

med_breast_df  = read.table(file.path(data_dir, 'media/breast_diet.tsv'), sep = '\t', header = TRUE) %>%
  select(reaction = 1, breast = 2)
med_formula_df = read.table(file.path(data_dir, 'media/formula_diet.tsv'), sep = '\t', header = TRUE) %>%
  select(reaction = 1, formula = 2)

media_compare = full_join(med_breast_df, med_formula_df, by = 'reaction') %>%
  mutate(breast  = replace_na(breast, 0),
         formula = replace_na(formula, 0)) %>%
  # substrates listed at zero in both media carry no information and are dropped
  filter(breast > 0 | formula > 0) %>%
  mutate(block = factor(case_when(
    breast > 0 & formula > 0 ~ "Shared",
    breast > 0               ~ "Breast only",
    TRUE                     ~ "Formula only"),
    levels = c("Shared","Breast only","Formula only")),
    label = ifelse(reaction %in% names(mediac_names),
                   mediac_names[reaction], reaction),
    fold = ifelse(breast > 0 & formula > 0, formula / breast, NA)) %>%
  arrange(block, desc(pmax(breast, formula))) %>%
  mutate(label = factor(label, levels = rev(unique(label))))

media_compare = media_compare %>%
  mutate(label = as.character(label),
         label = case_when(
           reaction == "EX_phyQ(e)"   ~ "Vitamin K1",
           reaction == "EX_avite1(e)" ~ "Vitamin E",
           reaction == "EX_3slac(e)"  ~ "3-Sialyllactose",
           reaction == "EX_dca(e)"  ~ "Dicarboxylic Acid",
           reaction == "EX_pro_D[e]"  ~ "D-Proline",
           reaction == "EX_asp_D[e]"  ~ "D-Aspartate",
           reaction == "EX_ala_D[e]"  ~ "D-Alanine",
           reaction == "EX_chol[e]"  ~ "Choline",
           reaction == "EX_chsterol(e)"  ~ "Cholesterol",
           reaction == "EX_ascb_L(e)"  ~ "Vitamin C",
           reaction == "EX_retinol(e)"  ~ "Vitamin A",
           reaction == "EX_retinol(e)"  ~ "Vitamin A",
           reaction == "EX_vitd3(e)"  ~ "Vitamin D3",
           reaction == "EX_hdcea(e)"  ~ "Palmitoyl",
           reaction == "EX_CE2510(e)"  ~ "CE2510",
           reaction == "EX_arachd(e)"  ~ "Arachidonic Acid",
           reaction == "EX_lacnfucpt_i(e)"  ~ "Lacto-N-fucopentaose I",
           reaction == "EX_lacnhx(e)"  ~ "Lacto-N-hexaose",
           reaction == "EX_caro(e)"  ~ "Carotenoid",
           reaction == "EX_octa(e)"  ~ "Octanoate",
           reaction == "EX_adpcbl(e)"  ~ "Adenosylcobalamin",
           TRUE ~ label)) %>%
  arrange(block, desc(pmax(breast, formula))) %>%
  mutate(label = factor(label, levels = rev(unique(label))))

media_compare = media_compare %>%
  mutate(label = as.character(label),
         # drop a trailing parenthesis and its contents, but keep one that opens the name
         label = trimws(sub("(?<=.)\\s*\\([^()]*\\)\\s*$", "", label, perl = TRUE))) %>%
  arrange(block, desc(pmax(breast, formula))) %>%
  mutate(label = factor(label, levels = rev(unique(label))))

write.table(media_compare,file = "Supplementary_Table_S12.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Zeros cannot be placed on a log axis, so they are drawn at a floor below the smallest non-zero value
conc_floor = min(c(media_compare$breast[media_compare$breast > 0],
                   media_compare$formula[media_compare$formula > 0])) / 10

media_long = media_compare %>%
  # the direction of each difference is determined before the two media are folded into one column
  mutate(direction = case_when(breast > formula ~ "Breast",
                               formula > breast ~ "Formula",
                               TRUE             ~ "Equal")) %>%
  pivot_longer(c(breast, formula), names_to = 'medium', values_to = 'conc') %>%
  mutate(medium = factor(medium, levels = c('breast','formula'),
                         labels = c('Breast','Formula')),
         decades = log10(pmax(conc, conc_floor) / conc_floor),
         x = ifelse(medium == 'Breast', -decades, decades),
         # only the arm carrying the higher concentration is coloured
         higher_in = factor(ifelse(as.character(medium) == direction,
                                   direction, "Equal"),
                            levels = c("Breast","Formula","Equal")))

# Axis labels are the concentrations the decades correspond to, mirrored about the centre
decade_max = ceiling(max(media_long$decades))
decade_breaks = seq(-decade_max, decade_max, by = 2)

lollipop_colors = c(diet_colors, "Equal" = "grey65")

media_lollipop = ggplot(media_long, aes(x = x, y = label, colour = higher_in))+
  geom_segment(aes(x = 0, xend = x, yend = label), linewidth = 0.8)+
  geom_point(size = 3.2)+
  geom_vline(xintercept = 0, colour = 'grey30', linewidth = 0.5)+
  ggforce::facet_col(~ block, scales = 'free_y', space = 'free')+
  scale_colour_manual(values = lollipop_colors, name = "Higher in")+
  scale_x_continuous(breaks = decade_breaks,
                     labels = function(b) {
                       v = conc_floor * 10^abs(b)
                       ifelse(abs(b) < 0.5, "0", format(signif(v, 1), scientific = TRUE))
                     })+
  theme_bw()+
  labs(x = "Concentration (mM)", y = "")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
        strip.text = element_text(size = 14),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = 'bottom')

ggsave(filename = 'Figures/Media_comparison_lollipop.jpeg',plot = media_lollipop,
       dpi = 1600,width = 16, height = 40, units = 'cm')

# Check the block sizes before reading the figure
table(media_compare$block)

# ---------------------------------------------------------------------------------------------------------------- #
# Human milk oligosaccharide utilisation under the breast diet. Consumption of each HMO is quantified per individual
# as the percentage decrease between the first and the last simulation step, and an individual is classified as a
# utiliser when consumption exceeds 5%. Utilisation is bimodal, with most individuals consuming nothing while a
# subset consumes nearly all of the substrate, so the number of utilisers is reported alongside the median.
# BacArena returns each HMO under both a long and an abbreviated label; these are merged before the calculation.
# The number of models in the library carrying each exchange reaction is reported so that low utilisation caused by
# limited library coverage can be distinguished from a biological result.
# ---------------------------------------------------------------------------------------------------------------- #

# ---- Per-individual HMO consumption ------------------------------------------------------------------------------
# BacArena returns each HMO under both a long and an abbreviated label, so the two are merged before the
# consumption of each is computed
hmo_consumption = function(extracted) {
  extracted$hmo_curve %>%
    mutate(sub = ifelse(sub %in% names(hmo_name_map), hmo_name_map[sub], sub)) %>%
    group_by(ID, sub, time) %>%
    summarise(value = sum(value), .groups = 'drop') %>%
    group_by(ID, sub) %>%
    summarise(initial = value[time == min(time)],
              final   = value[time == max(time)],
              consumed_pct = 100 * (initial - final) / initial, .groups = 'drop')
}

hmo_use = hmo_consumption(breast_ext) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = group_levels),
         sub   = factor(sub, levels = hmo_levels),
         user  = consumed_pct > utiliser_cutoff,
         bif_amount = bif_abund[ID])

# ---- Summaries ---------------------------------------------------------------------------------------------------
hmo_use_summary = hmo_use %>%
  group_by(sub) %>%
  summarise(n_users = sum(user, na.rm = TRUE),
            pct_users = 100 * n_users / n(),
            median_among_users = median(consumed_pct[user], na.rm = TRUE),
            max_consumed = max(consumed_pct, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(n_users))

hmo_use_bygroup = hmo_use %>%
  group_by(sub, group) %>%
  summarise(n = n(),
            n_users = sum(user, na.rm = TRUE),
            pct_users = 100 * n_users / n,
            median_consumed = median(consumed_pct, na.rm = TRUE),
            median_among_users = median(consumed_pct[user], na.rm = TRUE),
            .groups = 'drop')

# How many of the modelled strains carry each oligosaccharide exchange. Reported alongside the utilisation results
# so that a low utilisation rate can be told apart from a gap in the model library.
hmo_model_coverage = data.frame(
  sub = factor(hmo_levels, levels = hmo_levels),
  reaction = unname(hmo_reactions),
  n_models = sapply(hmo_reactions, function(rx)
    sum(sapply(breast_con_GEMs_list, function(m) rx %in% m@react_id))),
  row.names = NULL)

# ---- Which taxa actually take up the oligosaccharides -------------------------------------------------------------
# Net flux over the whole simulation; negative values are uptake. The division of labour between taxa follows the
# structural class of the oligosaccharide, so the shares are reported per HMO rather than pooled.
hmo_uptake = breast_ext$flux %>%
  filter(rea %in% unname(hmo_reactions)) %>%
  group_by(ID, spec, rea) %>%
  summarise(net = sum(mflux, na.rm = TRUE), .groups = 'drop') %>%
  filter(net < 0) %>%
  mutate(sub = factor(names(hmo_reactions)[match(rea, hmo_reactions)], levels = hmo_levels),
         uptake = -net,
         taxon_group = factor(case_when(
           grepl("Bifidobacterium", spec)                          ~ "Bifidobacterium",
           grepl("Bacteroides|Parabacteroides|Prevotella", spec)   ~ "Bacteroidetes",
           grepl("Lactobacillus|Enterococcus|Streptococcus", spec) ~ "Lactic acid bacteria",
           TRUE ~ "Other"),
           levels = c("Bifidobacterium","Bacteroidetes","Lactic acid bacteria","Other"))) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = group_levels))

# Consumers ranked across the cohort
hmo_consumer_rank = hmo_uptake %>%
  group_by(spec, taxon_group) %>%
  summarise(n_subjects = n_distinct(ID),
            n_hmo = n_distinct(sub),
            total_uptake = sum(uptake), .groups = 'drop') %>%
  arrange(desc(total_uptake))

# Leading consumers of each oligosaccharide within each subject group
hmo_consumer_top = hmo_uptake %>%
  group_by(sub, group, spec, taxon_group) %>%
  summarise(n_subjects = n_distinct(ID),
            mean_uptake = mean(uptake), .groups = 'drop') %>%
  group_by(sub, group) %>%
  slice_max(mean_uptake, n = 3) %>%
  ungroup()

# Share of the total uptake of each oligosaccharide carried by each taxonomic group
hmo_uptake_class = hmo_uptake %>%
  group_by(sub, taxon_group) %>%
  summarise(total = sum(uptake), .groups = 'drop') %>%
  group_by(sub) %>%
  mutate(pct = 100 * total / sum(total)) %>%
  ungroup()

hmo_uptake_class_wide = hmo_uptake_class %>%
  select(-total) %>%
  mutate(pct = round(pct)) %>%
  pivot_wider(names_from = taxon_group, values_from = pct, values_fill = 0)

# ---- Breadth of HMO use per infant and the taxa that predict it ----------------------------------------------------
# Each infant is scored by how many of the eight oligosaccharides its community can use, and this breadth is related
# to the inocula of the two genera that carry the uptake. Both are entered together because they co-occur only
# partly and cover different structural classes.
hmo_breadth = hmo_use %>%
  group_by(ID, group) %>%
  summarise(n_hmo_used = sum(user, na.rm = TRUE), .groups = 'drop') %>%
  left_join(taxa_features %>%
              mutate(bac_inoculum = genus_inoculum(c("Bacteroides","Parabacteroides"))) %>%
              select(ID, bif_inoculum, bac_inoculum), by = 'ID')

hmo_breadth_cor = data.frame(
  taxon = c("Bifidobacterium", "Bacteroidetes"),
  rho = c(cor(hmo_breadth$n_hmo_used, hmo_breadth$bif_inoculum, method = 'spearman'),
          cor(hmo_breadth$n_hmo_used, hmo_breadth$bac_inoculum, method = 'spearman')),
  p = c(cor.test(hmo_breadth$n_hmo_used, hmo_breadth$bif_inoculum,
                 method = 'spearman', exact = FALSE)$p.value,
        cor.test(hmo_breadth$n_hmo_used, hmo_breadth$bac_inoculum,
                 method = 'spearman', exact = FALSE)$p.value),
  row.names = NULL)

# The two genera themselves
bif_bac_cor = cor.test(hmo_breadth$bif_inoculum, hmo_breadth$bac_inoculum,
                       method = 'spearman', exact = FALSE)

# Joint model: do the two genera contribute independently?
hmo_breadth_model = lm(n_hmo_used ~ log1p(bif_inoculum) + log1p(bac_inoculum),
                       data = hmo_breadth)
hmo_breadth_model_summary = summary(hmo_breadth_model)

# ---- Bifidobacterium colonisation --------------------------------------------------------------------------------
bif_bygroup = data.frame(ID = names(bif_abund), bif = as.numeric(bif_abund),
                         stringsAsFactors = FALSE) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  group_by(group) %>%
  summarise(n = n(), median_bif = median(bif),
            q25 = quantile(bif, .25), q75 = quantile(bif, .75),
            n_zero = sum(bif == 0), pct_zero = 100 * n_zero / n, .groups = 'drop')

hmo_bif_link = hmo_use %>%
  group_by(sub, user) %>%
  summarise(n = n(), median_bif = median(bif_amount, na.rm = TRUE), .groups = 'drop')

# ---- Statistics --------------------------------------------------------------------------------------------------
# Consumption between the three subject groups
hmo_group_stat = hmo_use %>%
  group_by(sub) %>%
  rstatix::kruskal_test(consumed_pct ~ group) %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

hmo_pair_stat = hmo_use %>%
  group_by(sub) %>%
  rstatix::wilcox_test(consumed_pct ~ group, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  rstatix::add_xy_position(x = 'sub', dodge = 0.75, step.increase = 0.10)

# Bifidobacterium inoculum between the three subject groups
bif_df = data.frame(bif = as.numeric(bif_abund), ID = names(bif_abund),
                    stringsAsFactors = FALSE) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = group_levels))

bif_group_stat = bif_df %>% rstatix::kruskal_test(bif ~ group)

bif_pair_stat = bif_df %>%
  rstatix::wilcox_test(bif ~ group) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

# Bifidobacterium inoculum of utilisers versus non-utilisers
hmo_bif_stat = hmo_use %>%
  group_by(sub) %>%
  rstatix::wilcox_test(bif_amount ~ user) %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

# ---- Plots -------------------------------------------------------------------------------------------------------
hmo_use_plot = ggplot(hmo_use, aes(x = sub, y = consumed_pct, fill = group))+
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.size = 1,
               position = position_dodge(width = 0.75))+
  scale_fill_manual(values = group_colors, name = "Group")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  stat_pvalue_manual(hmo_pair_stat %>% filter(p.adj < 0.05),
                     label = 'p.adj.signif', tip.length = 0.01,
                     size = 7, bracket.size = 0.4)+
  theme_bw()+
  labs(x = "", y = "HMO consumed (%)")+
  theme(axis.text = element_text(size = 13), axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 17),legend.text = element_text(size = 15),
        legend.position = 'bottom')
ggsave(plot = hmo_use_plot,filename = 'Figures/hmo_use_plot.jpeg',
       dpi = 1200,width = 13,height = 13,units = 'cm')

hmo_user_plot = ggplot(hmo_use_bygroup, aes(x = sub, y = pct_users, fill = group))+
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8),
           width = 0.7, alpha = 0.85)+
  scale_fill_manual(values = group_colors, name = "Group")+
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05)))+
  theme_bw()+
  labs(x = "", y = paste0("Infants with >", utiliser_cutoff, "% consumption (%)"))+
  theme(axis.text = element_text(size = 13), axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 17),legend.text = element_text(size = 15),
        axis.title.y.left = element_text(size = 15),legend.position = 'bottom')
ggsave(plot = hmo_user_plot,filename = 'Figures/hmo_user_plot.jpeg',
       dpi = 1200,width = 13,height = 13,units = 'cm')


# Which taxonomic group carries the uptake of each oligosaccharide
hmo_uptake_plot = ggplot(hmo_uptake_class, aes(x = sub, y = pct, fill = taxon_group))+
  geom_bar(stat = 'identity', width = 0.7)+
  scale_fill_manual(values = c("Bifidobacterium" = "#e31a1c",
                               "Bacteroidetes" = "#1f78b4",
                               "Lactic acid bacteria" = "#33a02c",
                               "Other" = "grey70"), name = "")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  labs(x = "", y = "Share of total HMO uptake (%)")+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.position = 'bottom')+
  guides(color = guide_legend(nrow = 2))
ggsave(plot = hmo_uptake_plot,filename = 'Figures/hmo_uptake_plot.jpeg',
       dpi = 1200,height = 13,width = 15,units = 'cm')

# Proportion of infants whose community contains no Bifidobacterium at all. The inoculum is zero in most preterm
# infants, so a boxplot conveys little; presence or absence of the genus is what the utilisation results turn on.
bif_zero_plot = ggplot(bif_bygroup, aes(x = group, y = pct_zero, fill = group))+
  geom_bar(stat = 'identity', width = 0.6, alpha = 0.85)+
  geom_text(aes(label = paste0(n_zero, "/", n)), vjust = -0.4, size = 3.5)+
  scale_fill_manual(values = group_colors, name = "Group")+
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.08)))+
  theme_bw()+
  labs(x = "", y = "Infants without Bifidobacterium (%)")+
  theme(axis.text = element_text(size = 11), legend.position = 'none')

# Breadth of HMO use against the inoculum of each genus
hmo_breadth_plot = hmo_breadth %>%
  pivot_longer(c(bif_inoculum, bac_inoculum), names_to = 'taxon', values_to = 'inoculum') %>%
  mutate(taxon = factor(taxon, levels = c('bif_inoculum','bac_inoculum'),
                        labels = c('Bifidobacterium','Bacteroidetes'))) %>%
  ggplot(aes(x = inoculum + 1, y = n_hmo_used, colour = group))+
  geom_jitter(width = 0.06, height = 0.18, size = 2, alpha = 0.75)+
  geom_smooth(aes(group = 1), method = 'lm', se = TRUE, colour = 'grey30', linewidth = 0.6)+
  facet_wrap(~ taxon, nrow = 1)+
  scale_x_log10()+
  scale_colour_manual(values = group_colors, name = "Group")+
  theme_bw()+
  labs(x = "Inoculum + 1", y = "Number of HMOs utilised")+
  theme(axis.text = element_text(size = 11), strip.text = element_text(size = 11))


hmo_quadrant = hmo_breadth %>%
  mutate(quadrant = factor(case_when(
    bif_inoculum > 0 & bac_inoculum > 0   ~ "Both",
    bif_inoculum > 0 & bac_inoculum == 0  ~ "Bifidobacterium only",
    bif_inoculum == 0 & bac_inoculum > 0  ~ "Bacteroidetes only",
    TRUE                                  ~ "Neither"),
    levels = c("Neither","Bacteroidetes only","Bifidobacterium only","Both")))

hmo_quadrant_plot = ggplot(hmo_quadrant, aes(x = quadrant, y = n_hmo_used, fill = quadrant))+
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA)+
  geom_jitter(width = 0.15, height = 0.15, size = 1.5, alpha = 0.6, colour = 'grey25')+
  scale_y_continuous(breaks = seq(0, 8, 2))+
  scale_fill_brewer(palette = 'Set2')+
  theme_bw()+
  labs(x = "Colonisation by the two HMO-degrading groups", y = "Number of HMOs utilised")+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank())
ggsave(plot = hmo_quadrant_plot,filename = 'Figures/hmo_quadrant_plot.jpeg',
       dpi = 1200,width = 13,height = 13,units = 'cm')


# Bifidobacterium inoculum of utilisers and non-utilisers, kept for the supplementary material
hmo_bif_plot = ggplot(hmo_use, aes(x = sub, y = bif_amount + 1, fill = user))+
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1,
               position = position_dodge(width = 0.75))+
  scale_y_log10()+
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "deepskyblue1"),
                    labels = c("FALSE" = "Non-utiliser", "TRUE" = "Utiliser"),
                    name = "")+
  theme_bw()+
  labs(x = "", y = "Bifidobacterium inoculum + 1")+
  theme(axis.text = element_text(size = 11), axis.title.y = element_text(size = 11))

# ---- Supplementary tables ----------------------------------------------------------------------------------------
supp_table_hmo_coverage = hmo_model_coverage %>%
  select(HMO = sub, `Exchange reaction` = reaction, `Modelled strains carrying it` = n_models)

supp_table_hmo_tests = hmo_pair_stat %>%
  select(HMO = sub, Group1 = group1, Group2 = group2, n1, n2,
         statistic, p, p.adj, p.adj.signif)


supp_table_hmo_consumers = hmo_uptake %>%
  group_by(spec, taxon_group, sub) %>%
  summarise(n_subjects = n_distinct(ID),
            uptake = sum(uptake), .groups = 'drop') %>%
  group_by(sub) %>%
  mutate(`Share of uptake (%)` = round(100 * uptake / sum(uptake), 1)) %>%
  ungroup() %>%
  mutate(Taxon = spec %>% gsub("_(sp|subsp)_", " sp. ", .) %>% gsub("_", " ", .),
         `Total uptake (a.u.)` = signif(uptake, 3)) %>%
  select(HMO = sub, Taxon, `Taxonomic group` = taxon_group,
         `Subjects with net uptake` = n_subjects,
         `Share of uptake (%)`, `Total uptake (a.u.)`) %>%
  arrange(HMO, desc(`Share of uptake (%)`))

write.table(supp_table_hmo_consumers, "Supplementary_Table_S13.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ================================================================================================================== #
# Part 6  Probiotic intervention
#
#   Manuscript: Results 4; Figure 3F to 3J; Figures S11 to S14
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Probiotic intervention. Seven probiotic strains are introduced into the baseline in silico community of each
# preterm infant at an initial abundance of 10 per strain, following the design used in our previous work. Five
# strains reproduce the combination used there (L. reuteri, B. animalis subsp. lactis BB-12, L. rhamnosus GG,
# L. fermentum, C. butyricum), two of which (L. rhamnosus GG ATCC53103 and B. lactis Bb-12) are also among those
# conditionally recommended by the ESPGHAN 2020 position paper on probiotics for preterm infants. Two further
# bifidobacteria are added on the basis of the HMO utilisation deficit identified here: B. longum subsp. infantis
# ATCC 15697, corresponding to the B. infantis Bb-02 component of the ESPGHAN recommendation, and B. bifidum BGN4,
# which carries the largest number of HMO exchange reactions in the library. The exact commercial strains (Bb-02,
# Bb-12) are unavailable in AGORA2, so same-species reconstructions are substituted. Probiotic models are
# constrained with the same troubleshot media used for the baseline simulations, without additional substrates, so
# that treated and untreated arenas remain directly comparable. Simulations are executed on the cluster through
# simulate_arena_breast_pro.R and simulate_arena_formula_pro.R; the arena construction is reproduced below.
# ---------------------------------------------------------------------------------------------------------------- #
stopifnot(all(file.exists(probiotics)))

PRO_GEMs = list()
for (nm in names(probiotics)) {
  m = readMATmod_fix(probiotics[nm])
  m@mod_id = nm; m@mod_name = nm; m@mod_desc = basename(probiotics[nm])
  PRO_GEMs[[nm]] = m
}

PRO_cons_breast = GEMs_constraint(GEMs_list = PRO_GEMs, diet = new_breast)
#saveRDS(PRO_cons_breast, file = 'breast_PRO_GEMs_needed.rds')
PRO_cons_breast = readRDS(file.path(data_dir, 'derived/breast_PRO_GEMs_needed.rds'))
PRO_cons_formula = GEMs_constraint(GEMs_list = PRO_GEMs, diet = new_formula)
#saveRDS(PRO_cons_formula, file = 'formula_PRO_GEMs_needed.rds')
PRO_cons_formula = readRDS(file.path(data_dir, 'derived/formula_PRO_GEMs_needed.rds'))

# HMO capacity and growth feasibility of each probiotic strain on the two media
hmo_rx = c("EX_2fuclac(e)","EX_3fuclac(e)","EX_lacnttr(e)",
           "EX_lacnnttr(e)","EX_6slac(e)","EX_3slac(e)")
pro_profile = data.frame(
  strain = names(PRO_GEMs),
  species = sub("\\.mat$", "", sub("^PRO_\\d+_", "", basename(probiotics))),
  n_reactions = sapply(PRO_GEMs, function(m) m@react_num),
  n_hmo = sapply(PRO_GEMs, function(m) sum(hmo_rx %in% m@react_id)),
  hmo_list = sapply(PRO_GEMs, function(m)
    paste(gsub("EX_|\\(e\\)", "", hmo_rx[hmo_rx %in% m@react_id]), collapse = "; ")),
  growth_breast = sapply(PRO_cons_breast, function(m) optimizeProb(m)@lp_obj),
  growth_formula = sapply(PRO_cons_formula, function(m) optimizeProb(m)@lp_obj),
  row.names = NULL)

# Subjects receiving the intervention: preterm infants only
writeLines(preterm_ids, 'id_list_preterm.txt')


# ---------------------------------------------------------------------------------------------------------------- #
# Arena construction for the probiotic intervention. Reproduced here for reference; the cluster scripts perform the
# equivalent steps one subject at a time.
# ---------------------------------------------------------------------------------------------------------------- #
generate_arena_pro = function(comp_data, metadata, diet, map_data, cons_list, pro_list,
                              subject_ids, amount = pro_amount) {
  Ind_arena = list()
  for (Ind in subject_ids) {
    set.seed(sum(utf8ToInt(Ind)) + 1)
    temp_meta = subset(metadata, ID == Ind)
    temp_GEMs_amount = comp_data[rownames(temp_meta)]
    temp_GEMs_amount = temp_GEMs_amount[temp_GEMs_amount != 0,,drop = FALSE]
    
    add_bacs = list()
    for (i in c(names(cons_list), names(pro_list))) {
      mod = if (i %in% names(cons_list)) cons_list[[i]] else pro_list[[i]]
      add_bacs[[i]] = Bac(model = mod, limit_growth = TRUE, deathrate = 0.1,
                          minweight = runif(1, min = 0.2, max = 0.4),
                          maxweight = runif(1, min = 0.8, max = 1.0),
                          growtype = "exponential")
    }
    
    Ind_arena[[Ind]] = Arena(n = 100, m = 100)
    for (asv in rownames(temp_GEMs_amount)) {
      temp_GEM_input = map_data[map_data$ASV_ID == asv,]$renamed_gem
      Ind_arena[[Ind]] = addOrg(object = Ind_arena[[Ind]],
                                specI = add_bacs[[temp_GEM_input]],
                                amount = temp_GEMs_amount[asv,])
    }
    for (pro in names(pro_list)) {
      Ind_arena[[Ind]] = addOrg(object = Ind_arena[[Ind]],
                                specI = add_bacs[[pro]], amount = amount)
    }
    new_diet_ex = as.character(diet[,1])
    new_diet_val = as.numeric(diet[,2])
    for (met in 1:length(new_diet_ex)) {
      Ind_arena[[Ind]] = addSubs(object = Ind_arena[[Ind]],
                                 smax = ifelse(new_diet_val[met] > 1e-3,
                                               new_diet_val[met]/20, new_diet_val[met]),
                                 mediac = new_diet_ex[met], unit = 'mM', addAnyway = TRUE)
    }
  }
  return(Ind_arena)
}


# Read the simulated Eval objects of the probiotic-treated preterm infants
breast_pro_ext  = load_extracted('pro_breast', expected = 40)
formula_pro_ext = load_extracted('pro_formula', expected = 40)

# ---------------------------------------------------------------------------------------------------------------- #
# Effect of the probiotic intervention on the preterm infants. Engraftment of the seven introduced strains is
# examined first, so that an absent metabolic response can be attributed either to failed colonisation or to
# colonisation without additional output. HMO utilisation and SCFA production are then compared between the treated
# and untreated simulation of the same subject with paired Wilcoxon tests, and the effect of the intervention is
# contrasted between the two feeding backgrounds to test whether the presence of HMO modifies the response.
# ---------------------------------------------------------------------------------------------------------------- #

# ---- Engraftment of the introduced strains -----------------------------------------------------------------------
engraftment = function(extracted, diet_type, strains) {
  extracted$growth %>%
    filter(species %in% strains) %>%
    group_by(ID, species) %>%
    summarise(initial = value[time == min(time)],
              final = value[time == max(time)],
              fold_change = (final + 1) / (initial + 1), .groups = 'drop') %>%
    mutate(Diet = diet_type)
}

pro_engraft = rbind(engraftment(breast_pro_ext, 'Breast', pro_names),
                    engraftment(formula_pro_ext, 'Formula', pro_names)) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = c('S-Preterm','N-Preterm')),
         Diet = factor(Diet, levels = c('Breast','Formula')),
         species = factor(species, levels = pro_names))

pro_engraft_summary = pro_engraft %>%
  group_by(Diet, species) %>%
  summarise(n = n(), median_final = median(final), median_fold = median(fold_change),
            n_expanded = sum(final > initial), .groups = 'drop')

# Diet dependence of engraftment, tested within strain
engraft_stat = pro_engraft %>%
  arrange(species, Diet, ID) %>%
  group_by(species) %>%
  rstatix::wilcox_test(final ~ Diet, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

engraft_bracket = engraft_stat %>%
  filter(p.adj < 0.05) %>%
  left_join(pro_engraft %>% group_by(species) %>%
              summarise(top = quantile(final, 0.95), .groups = 'drop'), by = 'species') %>%
  mutate(x = match(as.character(species), levels(pro_engraft$species)),
         xmin = x - 0.19, xmax = x + 0.19, y.position = top * 1.10)

pro_engraft_plot = ggplot(pro_engraft, aes(x = species, y = final, fill = Diet))+
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1,
               position = position_dodge(width = 0.75))+
  geom_hline(yintercept = pro_amount, linetype = 'dashed', color = 'grey40')+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  stat_pvalue_manual(engraft_bracket, label = 'p.adj.signif',
                     tip.length = 0.01, size = 4, bracket.size = 0.4)+
  theme_bw()+
  labs(title = "Final abundance of the introduced probiotic strains",
       x = "", y = "Abundance at the final simulation step")+
  theme(axis.text.x = element_text(size = 11, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 12), title = element_text(size = 12))

# ---- SCFA before and after the intervention ----------------------------------------------------------------------
breast_pro_scfa = extract_scfa(breast_pro_ext, meta, 'Breast') %>%
  mutate(Status = 'Prob-treated')
formula_pro_scfa = extract_scfa(formula_pro_ext, meta, 'Formula') %>%
  mutate(Status = 'Prob-treated')

scfa_pro_all = rbind(
  breast_scfa  %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  formula_scfa %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  breast_pro_scfa, formula_pro_scfa) %>%
  mutate(Diet = factor(Diet, levels = c('Breast','Formula')),
         Status = factor(Status, levels = c('Non-treated','Prob-treated')),
         group = factor(group, levels = preterm_levels))

pro_effect_stat = scfa_pro_all %>%
  arrange(Diet, sub, Status, ID) %>%
  group_by(Diet, sub) %>%
  rstatix::wilcox_test(value ~ Status, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  group_by(Diet) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  ungroup()

pro_paired = scfa_pro_all %>%
  select(ID, group, Diet, sub, Status, value) %>%
  pivot_wider(names_from = Status, values_from = value) %>%
  mutate(Diff = `Prob-treated` - `Non-treated`)

pro_effect_summary = pro_paired %>%
  group_by(Diet, sub) %>%
  summarise(n = sum(!is.na(Diff)),
            median_untreated = median(`Non-treated`, na.rm = TRUE),
            median_treated = median(`Prob-treated`, na.rm = TRUE),
            median_diff = median(Diff, na.rm = TRUE),
            n_increased = sum(Diff > 0, na.rm = TRUE), .groups = 'drop') %>%
  left_join(pro_effect_stat %>% select(Diet, sub, p, p.adj, p.adj.signif), by = c('Diet','sub'))
write.table(pro_effect_summary, "Supplementary_Table_S14.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Does the feeding background modify the effect of the intervention?
pro_interaction_stat = pro_paired %>%
  arrange(sub, Diet, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(Diff ~ Diet, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

pro_paired_plot = ggplot(scfa_pro_all, aes(x = Status, y = value))+
  geom_line(aes(group = ID, color = group), alpha = 0.3, linewidth = 0.4)+
  geom_point(aes(color = group), size = 1, alpha = 0.6)+
  geom_boxplot(width = 0.35, alpha = 0.25, outlier.shape = NA, fill = 'grey70')+
  facet_grid(Diet ~ sub, scales = 'free_y')+
  scale_color_manual(values = group_colors, name = "Group")+
  theme_bw()+
  labs(title = "Predicted fermentation products before and after the probiotic intervention",
       x = "", y = "Fermentation product, log(mM+1)")+
  theme(strip.text = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 9), title = element_text(size = 12))

pro_diff_marks = pro_paired %>%
  group_by(sub) %>%
  summarise(y = max(Diff, na.rm = TRUE) * 1.15, .groups = 'drop') %>%
  filter(sub %in% c('Butyrate', 'Formate'))

pro_diff_plot = ggplot(pro_paired, aes(x = sub, y = Diff, fill = Diet))+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1,
               position = position_dodge(width = 0.75))+
  geom_text(data = pro_diff_marks, aes(x = sub, y = y, label = "\u2605"),
            inherit.aes = FALSE, size = 6, colour = 'goldenrod2')+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))+
  theme_bw()+
  labs(title = "Change in predicted fermentation products after the probiotic intervention",
       x = "", y = "Probiotic minus untreated, log(mM+1)")+
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        axis.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        title = element_blank(),
        legend.position = "bottom")
ggsave(plot = pro_diff_plot,filename = 'Figures/pro_diff_plot.jpeg',
       dpi = 1200,width = 13,height = 13,units = 'cm')



# ---- HMO utilisation before and after the intervention -----------------------------------------------------------
hmo_use_pro = hmo_consumption(breast_pro_ext) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = c('S-Preterm','N-Preterm')),
         sub = factor(sub, levels = hmo_levels),
         user = consumed_pct > utiliser_cutoff,
         Status = 'Prob-treated')

hmo_pro_compare = rbind(
  hmo_use %>% filter(ID %in% preterm_ids) %>%
    select(sub, ID, group, consumed_pct, user) %>% mutate(Status = 'Non-treated'),
  hmo_use_pro %>% select(sub, ID, group, consumed_pct, user, Status)) %>%
  mutate(Status = factor(Status, levels = c('Non-treated','Prob-treated')))

hmo_pro_stat = hmo_pro_compare %>%
  arrange(sub, Status, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(consumed_pct ~ Status, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

hmo_pro_summary = hmo_pro_compare %>%
  group_by(sub, Status) %>%
  summarise(n_users = sum(user, na.rm = TRUE), pct_users = 100 * n_users / n(),
            median_consumed = median(consumed_pct, na.rm = TRUE), .groups = 'drop')

hmo_pro_plot = ggplot(hmo_pro_compare, aes(x = sub, y = consumed_pct, fill = Status))+
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1,
               position = position_dodge(width = 0.75))+
  scale_fill_manual(values = c("Non-treated" = "grey70", "Prob-treated" = "deepskyblue1"), name = "")+
  theme_bw()+
  labs(x = "", y = "HMO consumed (%)")+
  theme(axis.text = element_text(size = 13), title = element_blank(),
        axis.title = element_text(size = 15),legend.text = element_text(size = 15),
        legend.position = 'bottom')
ggsave(plot = hmo_pro_plot,filename = 'Figures/hmo_pro_plot.jpeg',
       dpi = 1200,width = 13,height = 13,units = 'cm')

hmo_pro_user_plot = ggplot(hmo_pro_summary, aes(x = sub, y = pct_users, fill = Status))+
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8)+
  scale_fill_manual(values = c("Non-treated" = "grey70", "Prob-treated" = "deepskyblue1"), name = "")+
  theme_bw()+
  labs(x = "", y = paste0("Individuals with >", utiliser_cutoff, "% consumption (%)"))+
  theme(axis.text = element_text(size = 13), title = element_blank(),
        axis.title.y = element_text(size = 15),legend.text = element_text(size = 15),
        legend.position = "bottom")
ggsave(plot = hmo_pro_user_plot,filename = 'Figures/hmo_pro_user_plot.jpeg',
       dpi = 1200,width = 13,height = 13,units = 'cm')

but_producer = breast_pro_ext$flux %>%
  filter(rea == 'EX_but(e)') %>%
  group_by(ID, spec) %>%
  summarise(mflux_mean = mean(mflux, na.rm = TRUE), .groups = 'drop') %>%
  filter(mflux_mean > 1e-6) %>%
  group_by(spec) %>%
  summarise(n_subjects = n(), median_flux = median(mflux_mean), .groups = 'drop') %>%
  arrange(desc(median_flux))
but_producer



# ---------------------------------------------------------------------------------------------------------------- #
# Comparison of the SCFA levels before and after the probiotic treatment, carried out separately under each feeding
# condition. Individuals in which a given acid rises above a threshold of 1e-5 log(mM+1) after the intervention are
# classified as responders for that acid, so that the breadth of the response can be read per subject rather than
# only as a group median.
# ---------------------------------------------------------------------------------------------------------------- #

# breast_pro_scfa and formula_pro_scfa already carry the Status column; extracted once above

all_scfa_pro_treat = rbind(
  breast_scfa  %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  formula_scfa %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  breast_pro_scfa, formula_pro_scfa) %>%
  mutate(Status = factor(Status, levels = c('Non-treated','Prob-treated')),
         Diet = factor(Diet, levels = c('Breast','Formula')),
         group = factor(group, levels = c('S-Preterm','N-Preterm')),
         sub = factor(sub, levels = scfa_levels),
         ID = factor(ID, levels = unique(ID[order(group)])))

# ---- Per-individual comparison, treated overlaid on untreated ----------------------------------------------------
scfa_pro_overlay_plot = function(scfa_df, diet_type) {
  d = scfa_df %>% filter(Diet == diet_type)
  ggplot()+
    geom_bar(data = d %>% filter(Status == 'Prob-treated'),
             aes(x = value, y = ID, fill = interaction(group, Status)),
             stat = 'identity', position = 'identity', alpha = 0.7)+
    geom_bar(data = d %>% filter(Status == 'Non-treated'),
             aes(x = value, y = ID, fill = interaction(group, Status)),
             stat = 'identity', position = 'identity', alpha = 1)+
    facet_grid(. ~ sub, scales = 'free_x', space = 'fixed')+
    scale_fill_manual(
      values = c("N-Preterm.Non-treated"  = "orange",
                 "S-Preterm.Non-treated"  = "red",
                 "N-Preterm.Prob-treated" = "#FFB90F",
                 "S-Preterm.Prob-treated" = "brown4"),
      name = "Group",
      labels = c("N-Preterm.Non-treated"  = "N-Preterm (non-treated)",
                 "S-Preterm.Non-treated"  = "S-Preterm (non-treated)",
                 "N-Preterm.Prob-treated" = "N-Preterm (prob-treated)",
                 "S-Preterm.Prob-treated" = "S-Preterm (prob-treated)"))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 10),
          strip.text.y = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_markdown(size = 7),
          axis.text.x = element_text(angle = 60, vjust = 0.6, size = 8),
          panel.spacing = unit(0.2, "lines"),
          panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6),
          plot.margin = unit(c(0, -2000, 0, -2000), "pt"))+
    labs(x = "Fermentation product, log(mM+1)", y = "")+
    scale_y_discrete(position = 'right')
}

# ---- Identify responses of the probiotic treatment on SCFA restoration -------------------------------------------
pro_treat_response = all_scfa_pro_treat %>%
  select(ID, group, Diet, sub, Status, value) %>%
  pivot_wider(names_from = Status, values_from = value) %>%
  mutate(delta = `Prob-treated` - `Non-treated`,
         delta = ifelse(delta < resp_cutoff, 0, delta),
         present = ifelse(delta > resp_cutoff, "Present", "Absent"))

pro_treat_response_summary = pro_treat_response %>%
  group_by(Diet, sub) %>%
  summarise(n = n(),
            n_responded = sum(present == "Present"),
            pct_responded = 100 * n_responded / n,
            median_delta = median(delta), .groups = 'drop')

pro_treat_resp_plot = function(resp_df, diet_type) {
  ggplot(resp_df %>% filter(Diet == diet_type),
         aes(x = sub, y = ID, fill = present))+
    geom_tile(color = "white", width = 1.1, height = 0.9)+
    scale_fill_manual(values = c("Present" = "forestgreen", "Absent" = "lightgrey"),
                      labels = c("Present" = "Responded", "Absent" = "Non-responded"))+
    coord_fixed()+
    theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0, size = 6),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0, -2000, 0, -2000), "pt"))+
    scale_x_discrete(position = 'top')
}

breast_pro_treat_scfa_resp = scfa_pro_overlay_plot(all_scfa_pro_treat, 'Breast') +
  pro_treat_resp_plot(pro_treat_response, 'Breast') +
  plot_layout(guides = 'collect')
formula_pro_treat_scfa_resp = scfa_pro_overlay_plot(all_scfa_pro_treat, 'Formula') +
  pro_treat_resp_plot(pro_treat_response, 'Formula') +
  plot_layout(guides = 'collect')

ggsave(plot = breast_pro_treat_scfa_resp,
       filename = 'Figures/Pro_treat_scfa_resp_breast.jpeg',
       dpi = 1200, width = 24, height = 12, units = 'cm')
ggsave(plot = formula_pro_treat_scfa_resp,
       filename = 'Figures/Pro_treat_scfa_resp_formula.jpeg',
       dpi = 1200, width = 24, height = 12, units = 'cm')

# ---------------------------------------------------------------------------------------------------------------- #
# Metabolic activity of the SCFAs produced by the introduced probiotics in each individual. Fluxes are shown on a
# signed logarithmic scale because they span several orders of magnitude and can be negative when a strain consumes
# rather than secretes an acid. All seven acids are kept as fixed colour levels so that the panels share a single
# legend even when an individual produces only a subset of them, and the x axis text is retained only on the bottom
# row of the wrapped layout since every panel covers the same 0 to 12 h window.
# ---------------------------------------------------------------------------------------------------------------- #

pro_mflux_plots = function(flux_df, subject_ids, ncol = 2) {
  n_total = length(subject_ids)
  bottom_row = seq(n_total - ((n_total - 1) %% ncol), n_total)
  
  plot_list = list()
  for (i in seq_along(subject_ids)) {
    Ind = subject_ids[i]
    # The flux table carries the oligosaccharide exchanges as well, so the fermentation products are selected
    # explicitly before the signed log transform
    pro_scfa_data = flux_df %>%
      filter(ID == Ind, spec %in% names(pro_labels), rea %in% scfas) %>%
      rename(sub = rea) %>%
      select(spec, sub, time, mflux) %>%
      mutate(mflux = sign(mflux) * log(abs(mflux) + 1))
    
    # Fill in the acids this individual does not produce so that every panel carries the
    # same colour levels and patchwork collects a single legend
    pro_scfa_data = pro_scfa_data %>%
      mutate(spec = factor(spec, levels = names(pro_labels)),
             sub = factor(sub, levels = names(scfa_flux_labels)[
               match(scfa_levels, scfa_flux_labels)])) %>%
      complete(spec, sub, time = unique(time))
    
    p = ggplot(pro_scfa_data, aes(x = time, y = mflux, color = sub))+
      geom_line(linewidth = 0.7, na.rm = TRUE)+
      facet_wrap(~ spec, nrow = 1,
                 labeller = as_labeller(pro_labels))+
      scale_color_manual(values = scfa_colors, labels = scfa_flux_labels, name = "Fermentation product")+
      labs(x = 'Time (h)', y = Ind)+
      theme_bw()+
      theme(strip.text = element_text(size = 7, face = 'bold'),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            axis.text.y = element_text(size = 7),
            axis.title.y = element_text(size = 7, face = 'bold'),
            axis.title.x = element_blank())
    
    p = if (i %in% bottom_row) {
      p + theme(axis.text.x = element_text(size = 5))
    } else {
      p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    plot_list[[Ind]] = p
  }
  plot_list
}

sp_ids = meta$ID[meta$group == 'S-Preterm']
np_ids = meta$ID[meta$group == 'N-Preterm']

breast_sp_mflux = wrap_plots(pro_mflux_plots(breast_pro_ext$flux, sp_ids), ncol = 2, guides = 'collect')
breast_np_mflux = wrap_plots(pro_mflux_plots(breast_pro_ext$flux, np_ids), ncol = 2, guides = 'collect')
formula_sp_mflux = wrap_plots(pro_mflux_plots(formula_pro_ext$flux, sp_ids), ncol = 2, guides = 'collect')
formula_np_mflux = wrap_plots(pro_mflux_plots(formula_pro_ext$flux, np_ids), ncol = 2, guides = 'collect')

ggsave(plot = breast_sp_mflux, filename = "Figures/Breast_sepsis_pro_mflux_log.png",
       dpi = 1400, width = 36, height = 34, units = 'cm')
ggsave(plot = breast_np_mflux, filename = "Figures/Breast_preterm_pro_mflux_log.png",
       dpi = 1400, width = 36, height = 34, units = 'cm')
ggsave(plot = formula_sp_mflux, filename = "Figures/Formula_sepsis_pro_mflux_log.png",
       dpi = 1400, width = 36, height = 34, units = 'cm')
ggsave(plot = formula_np_mflux, filename = "Figures/Formula_preterm_pro_mflux_log.png",
       dpi = 1400, width = 36, height = 34, units = 'cm')


pro_flux_all = rbind(breast_pro_ext$flux %>% mutate(Diet = 'Breast'),
                     formula_pro_ext$flux %>% mutate(Diet = 'Formula')) %>%
  filter(spec %in% names(pro_labels)) %>%
  group_by(Diet, ID, spec, rea) %>%
  summarise(net = sum(mflux, na.rm = TRUE), .groups = 'drop')

pro_mflux_summary = pro_flux_all %>%
  filter(rea %in% scfas) %>%
  mutate(sub = factor(scfa_flux_labels[rea], levels = scfa_levels)) %>%
  group_by(Diet, sub, spec) %>%
  summarise(n_secretors = sum(net > 0),
            n_consumers = sum(net < 0),
            median_net = median(net),
            total_net = sum(net), .groups = 'drop')
write.table(pro_mflux_summary, file = "Supplementary_Table_S15.txt",
            row.names = FALSE, sep = '\t', quote = FALSE)

pro_hmo_flux = pro_flux_all %>%
  filter(rea %in% unname(hmo_reactions)) %>%
  mutate(sub = factor(names(hmo_reactions)[match(rea, hmo_reactions)], levels = hmo_levels)) %>%
  group_by(Diet, sub, spec) %>%
  summarise(n_consumers = sum(net < 0),
            median_uptake = median(-net), .groups = 'drop') %>%
  # every strain is kept, including those with no uptake at all, so that the table shows which strains
  # carry the oligosaccharide reactions and which do not
  complete(Diet, sub, spec = names(pro_labels),
           fill = list(n_consumers = 0, median_uptake = 0)) %>%
  mutate(spec = factor(spec, levels = names(pro_labels), labels = pro_labels))

supp_table_S17 = pro_hmo_flux %>%
  mutate(median_uptake = signif(median_uptake, 3)) %>%
  select(Medium = Diet, HMO = sub, Strain = spec,
         `Subjects with net uptake` = n_consumers,
         `Median uptake` = median_uptake) %>%
  arrange(Medium, HMO, desc(`Median uptake`))

write.table(supp_table_S17, "Supplementary_Table_S17.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------------------------------------------- #
# Net exchange flux of each probiotic strain for the fermentation products, under both nutritional conditions.
# Positive values are secretion and negative values uptake, so the cross-feeding between the lactate-producing
# bifidobacteria and the butyrate-producing C. butyricum can be read directly from the signs.
# ---------------------------------------------------------------------------------------------------------------- #
pro_mflux_plot_summary = pro_mflux_summary %>%
  mutate(Class = factor(scfa_class[as.character(sub)], levels = scfa_class_levels),
         spec = factor(spec, levels = c('PRO_BINF','PRO_BBIF','PRO_BB12',
                                        'PRO_LGG','PRO_FER','PRO_REU','PRO_CBU'))) %>%
  ggplot(aes(x = spec, y = median_net, fill = Diet))+
  geom_col(position = position_dodge(width = 0.6), width = 0.5)+
  geom_hline(yintercept = 0, colour = 'grey40', linewidth = 0.3)+
  ggh4x::facet_nested(~ Class + sub, scales = 'free_y')+
  scale_x_discrete(labels = pro_labels)+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  # Fluxes span several orders of magnitude in both directions and include zeros, so the axis is linear near zero
  # and logarithmic beyond it; sigma is set at the level below which a flux carries no weight in the comparison
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e3, base = 10),
                     breaks = c(-1e5, -1e4, 0, 1e4, 1e5, 1e6),
                     labels = scales::scientific)+
  theme_bw()+
  labs(x = "", y = "Median net flux")+
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.position = 'bottom')

ggsave(plot = pro_mflux_plot_summary,filename = 'Figures/pro_mflux_plot_summary.jpeg',
       dpi = 1200,width = 42,height = 13,units = 'cm')



# ---------------------------------------------------------------------------------------------------------------- #
# Cross-feeding network of the probiotic consortium. Each strain is linked to a short-chain fatty acid when it is a
# net secretor or a net consumer of that acid in the majority of individuals; edge width is scaled to the median
# net flux and direction indicates whether the acid is released or taken up. Links present in fewer than half of
# the individuals are omitted so that the network reflects consistent behaviour rather than isolated cases.
# ---------------------------------------------------------------------------------------------------------------- #
consistency_cutoff = 0.5   # a link is drawn when it holds in at least this fraction of individuals
n_preterm = length(preterm_ids)

crossfeed_edges = pro_mflux_summary %>%
  mutate(role = case_when(n_secretors >= consistency_cutoff * n_preterm ~ "Secretion",
                          n_consumers >= consistency_cutoff * n_preterm ~ "Uptake",
                          TRUE ~ NA_character_)) %>%
  filter(!is.na(role)) %>%
  mutate(from = ifelse(role == "Secretion", spec, sub),
         to   = ifelse(role == "Secretion", sub, spec),
         weight = abs(median_net),
         n_consistent = ifelse(role == "Secretion", n_secretors, n_consumers)) %>%
  select(Diet, from, to, role, weight, n_consistent, spec, sub)
write.table(supplement_table, file = "supplement_table.txt",
            row.names = FALSE, sep = '\t', quote = FALSE)


# Node table derived from the edges themselves so that it cannot fall out of step with the
# strain labels used upstream
crossfeed_nodes = data.frame(
  name = unique(c(crossfeed_edges$spec, crossfeed_edges$sub)),
  stringsAsFactors = FALSE) %>%
  mutate(type = ifelse(name %in% unique(crossfeed_edges$sub), "Metabolite", "Strain"))

crossfeed_plot = function(edge_df, node_df, diet_type) {
  e = edge_df %>% filter(Diet == diet_type) %>% select(from, to, role, weight, n_consistent)
  n = node_df %>% filter(name %in% c(e$from, e$to))
  stopifnot(all(c(e$from, e$to) %in% n$name))
  g = graph_from_data_frame(d = e, vertices = n, directed = TRUE)
  
  ggraph(g, layout = 'bipartite', types = (n$type == "Metabolite"))+
    geom_edge_arc(aes(width = log10(weight + 1), colour = role),
                  arrow = arrow(length = unit(2.5, 'mm'), type = 'closed'),
                  end_cap = circle(4, 'mm'), start_cap = circle(4, 'mm'),
                  strength = 0.15, alpha = 0.75)+
    geom_node_point(aes(shape = type, fill = type), size = 6, colour = 'grey20', stroke = 0.4)+
    geom_node_text(aes(label = name), repel = TRUE, size = 3)+
    scale_edge_width(range = c(0.3, 2.2), name = expression(log[10](flux)))+
    scale_edge_colour_manual(values = c("Secretion" = "#00AFBB", "Uptake" = "#E7B800"),
                             name = "")+
    scale_shape_manual(values = c("Strain" = 21, "Metabolite" = 22), name = "")+
    scale_fill_manual(values = c("Strain" = "white", "Metabolite" = "grey85"), name = "")+
    coord_flip()+
    theme_graph(base_family = 'sans')+
    labs(title = paste("Cross-feeding among the probiotic strains (", diet_type, ")", sep = ''))
}

breast_crossfeed_plot = crossfeed_plot(crossfeed_edges, crossfeed_nodes, 'Breast')
formula_crossfeed_plot = crossfeed_plot(crossfeed_edges, crossfeed_nodes, 'Formula')

ggsave('Figures/Crossfeed_network_breast.pdf', breast_crossfeed_plot,
       width = 20, height = 16, units = 'cm')
ggsave('Figures/Crossfeed_network_formula.pdf', formula_crossfeed_plot,
       width = 20, height = 16, units = 'cm')

# Edge list for the supplementary table
crossfeed_table = crossfeed_edges %>%
  select(Diet, strain = spec, metabolite = sub, role, n_consistent, median_net = weight) %>%
  arrange(Diet, strain, metabolite)


# ---------------------------------------------------------------------------------------------------------------- #
# Where the lactate that C. butyricum consumes comes from. The panel above covers the seven added strains only, so
# the resident members of the community are included here to show that the substrate was already being produced
# before the intervention; the bottleneck lay in the conversion step rather than in the supply.
# ---------------------------------------------------------------------------------------------------------------- #
lactate_source = rbind(breast_pro_ext$flux %>% mutate(Diet = 'Breast'),
                       formula_pro_ext$flux %>% mutate(Diet = 'Formula')) %>%
  filter(rea == "EX_lac_L(e)") %>%
  group_by(Diet, ID, spec) %>%
  summarise(net = sum(mflux, na.rm = TRUE), .groups = 'drop') %>%
  filter(net > 0) %>%
  group_by(Diet, spec) %>%
  summarise(n_subjects = n_distinct(ID),
            median_net = median(net), .groups = 'drop') %>%
  # the added strains are labelled and marked, the resident members keep their strain names
  mutate(origin = ifelse(spec %in% names(pro_labels), "Added Probiotics", "Resident microbe"),
         taxon = ifelse(spec %in% names(pro_labels),
                        pro_labels[spec],
                        spec %>% gsub("_(sp|subsp)_", " sp. ", .) %>% gsub("_", " ", .)))

# Only the strains that secrete lactate in at least five infants are shown, so the figure is not dominated by
# single-subject contributors
lactate_source_plot = lactate_source %>%
  filter(n_subjects >= 5) %>%
  ggplot(aes(x = median_net, y = reorder(taxon, median_net),
             colour = origin, size = n_subjects))+
  geom_point(alpha = 0.85)+
  facet_wrap(~ Diet, nrow = 1)+
  scale_x_continuous(labels = scales::scientific)+
  scale_colour_manual(values = c("Added Probiotics" = "#e31a1c",
                                 "Resident microb" = "grey20"), name = "")+
  scale_size_continuous(range = c(2, 8), name = "Individuals")+
  theme_bw()+
  labs(x = "Median net flux of (S)-lactate", y = "")+
  theme(axis.text.y = element_text(size = 16, face = 'italic'),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.spacing.x = unit(1.5, "cm"),
        legend.position = 'bottom')

ggsave('Figures/lactate_source_plot.jpeg', lactate_source_plot,
       dpi = 1200, width = 35, height = 15, units = 'cm')

supp_table_lactate_source = lactate_source %>%
  mutate(median_net = signif(median_net, 3)) %>%
  select(Medium = Diet, Taxon = taxon, Origin = origin,
         `Subjects with net secretion` = n_subjects,
         `Median net flux` = median_net) %>%
  arrange(Medium, desc(`Median net flux`))
write.table(supp_table_lactate_source, file = "Supplementary_Table_S16.txt",
            row.names = FALSE, sep = '\t', quote = FALSE)

# ================================================================================================================== #
# Part 7  Determinacy of the flux attribution
#
#   Manuscript: Results 2, closing paragraph; Figure S10
#   Three sensitivity analyses, run together because they share the constrained model
#   lists built in Parts 1 and 6. library() calls removed; Part 0 already loads them.
# ================================================================================================================== #

# ##################################################################################################################
# 7.1  Alternate-optima envelope of the exchange fluxes
# ##################################################################################################################

# ================================================================================================== #
# Sensitivity analysis 1: alternate-optima envelope of the exchange fluxes
#
# Purpose
#   The flux attribution reported in Results sections 2, 4 and 6 is read off a single optimal solution
#   returned by the LP solver at each time step. Linear programs of this kind normally admit a whole
#   face of alternate optima, so a reviewer is entitled to ask whether a different optimal solution
#   would have given a different attribution. This script answers that question at the level where it
#   can be answered exactly, by bounding what each exchange flux COULD have been rather than reporting
#   what one particular solution happened to give.
#
#   Two bounds are computed per model per medium.
#
#     capability : biomass is left free and the exchange flux is minimised and maximised over the
#                  entire feasible space allowed by the medium. A maximum of zero here is a hard
#                  negative: the strain cannot secrete that product under any optimal or suboptimal
#                  solution, at any time step, in any simulation. This is the bound that certifies
#                  statements of the form "strain X was the only producer of Y".
#
#     obligate   : biomass is pinned at its optimum and the exchange flux is minimised. A minimum
#                  above zero is a hard positive: whenever the strain grows at its maximum rate it
#                  must secrete that product. This certifies statements of the form "X produced Y".
#
#   Fluxes whose capability maximum is positive but whose obligate minimum is zero are the degenerate
#   cases. For those the reported median net flux is solution dependent and the manuscript should
#   report the attribution qualitatively rather than by magnitude.
#
# Prerequisites
#   Run after Part 1 and Part 5 of major_script.R, with these objects in the workspace:
#     breast_con_GEMs_list, formula_con_GEMs_list   resident models constrained to each medium
#     PRO_cons_breast, PRO_cons_formula             the seven probiotic models, same constraints
#     scfas, hmo_reactions, pro_labels              defined in Part 0
#
# Runtime
#   Roughly 15 exchanges x 2 directions x 2 modes x 2 media x ~156 models = 37k small LPs.
#   Expect minutes, not hours, on glpk. No re-simulation is required.
#
# Output
#   sensitivity/flux_envelope_raw.rds        full table, every model x medium x reaction
#   sensitivity/Supplementary_Table_S22.txt  the probiotic subset, for the manuscript
#   console verdicts for the four structural claims listed at the bottom
# ================================================================================================== #


dir.create('sensitivity', showWarnings = FALSE)

target_reactions = c(scfas, unname(hmo_reactions))

# -------------------------------------------------------------------------------------------------- #
# Core routine. Returns one row per requested reaction. A reaction absent from the model is reported
# as present = FALSE with NA bounds, which is a stronger negative than a zero bound and is kept
# distinct from it on purpose.
# -------------------------------------------------------------------------------------------------- #
exchange_envelope = function(model, reactions, mode = c('capability', 'obligate'),
                             tol = 1e-9) {
  mode = match.arg(mode)
  rid  = sybil::react_id(model)
  bm   = which(sybil::obj_coef(model) != 0)
  
  out = tibble(reaction = reactions, present = reactions %in% rid,
               min_flux = NA_real_, max_flux = NA_real_,
               biomass = NA_real_, status = NA_character_)
  
  opt = try(sybil::optimizeProb(model), silent = TRUE)
  if (inherits(opt, 'try-error') || opt@lp_stat != 5) {
    out$status = 'infeasible_base'
    return(out)
  }
  out$biomass = opt@lp_obj
  
  work = model
  if (mode == 'obligate') {
    if (opt@lp_obj < 1e-6) {
      out$status = 'no_growth'
      return(out)
    }
    # Pinning the lower bound rather than fixing both bounds keeps the LP feasible when the solver
    # returns an objective value a few ulp below the true optimum.
    work = sybil::changeBounds(work, react = bm, lb = opt@lp_obj * (1 - 1e-6))
  }
  
  for (i in which(out$present)) {
    idx = match(reactions[i], rid)
    m2  = sybil::changeObjFunc(work, react = idx, obj_coef = 1)
    lo  = try(sybil::optimizeProb(m2, lpdir = 'min'), silent = TRUE)
    hi  = try(sybil::optimizeProb(m2, lpdir = 'max'), silent = TRUE)
    if (!inherits(lo, 'try-error') && lo@lp_stat == 5) out$min_flux[i] = lo@lp_obj
    if (!inherits(hi, 'try-error') && hi@lp_stat == 5) out$max_flux[i] = hi@lp_obj
    out$status[i] = 'ok'
  }
  # Values within solver tolerance of zero are reported as zero. Leaving them at 1e-13 would make
  # every "cannot secrete" verdict below read as a vanishingly small positive.
  out = out %>%
    mutate(min_flux = ifelse(abs(min_flux) < tol, 0, min_flux),
           max_flux = ifelse(abs(max_flux) < tol, 0, max_flux))
  out
}

run_envelope_set = function(model_list, medium_label, reactions = target_reactions) {
  bind_rows(lapply(names(model_list), function(nm) {
    bind_rows(
      exchange_envelope(model_list[[nm]], reactions, 'capability') %>% mutate(mode = 'capability'),
      exchange_envelope(model_list[[nm]], reactions, 'obligate')   %>% mutate(mode = 'obligate')
    ) %>% mutate(model = nm, medium = medium_label, .before = 1)
  }))
}

# -------------------------------------------------------------------------------------------------- #
# The seven introduced strains first. These carry the four claims that the manuscript states
# categorically, so they are the ones that have to be certified.
# -------------------------------------------------------------------------------------------------- #
message('Probiotic strains, breast milk medium')
pro_env_breast  = run_envelope_set(PRO_cons_breast,  'Breast')
message('Probiotic strains, formula medium')
pro_env_formula = run_envelope_set(PRO_cons_formula, 'Formula')

# The seven strains constrained to the HMO-supplemented formula medium, which is the background of the
# 2 x 2 factorial. new_formula_hmo is read in Part 1 and is built by simulation/prepare_hmo_formula_diet.R.
PRO_cons_formula_hmo = GEMs_constraint(GEMs_list = PRO_GEMs, diet = new_formula_hmo)

pro_env_formula_hmo = run_envelope_set(PRO_cons_formula_hmo, 'Formula+HMO') %>%
  mutate(strain = unname(pro_labels[model]))

pro_env_formula_hmo %>%
  filter(reaction %in% unname(hmo_reactions), mode == 'capability') %>%
  group_by(strain) %>%
  summarise(n_present = sum(present),
            n_can_take_up = sum(present & min_flux < 0, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(n_can_take_up))

pro_envelope = bind_rows(pro_env_breast, pro_env_formula, pro_env_formula_hmo) %>%
  mutate(strain = unname(pro_labels[model]),
         sub = case_when(
           reaction %in% names(scfa_flux_labels) ~ unname(scfa_flux_labels[reaction]),
           reaction %in% unname(hmo_reactions)   ~ names(hmo_reactions)[match(reaction, hmo_reactions)],
           TRUE ~ reaction))

saveRDS(pro_envelope, 'sensitivity/flux_envelope_probiotics.rds')

# -------------------------------------------------------------------------------------------------- #
# The resident community models. Slower because there are ~149 of them, and only needed to support
# the community-level attribution statements in section 2. Set run_residents to FALSE to skip.
# -------------------------------------------------------------------------------------------------- #
run_residents = TRUE

if (run_residents) {
  message('Resident models, breast milk medium (', length(breast_con_GEMs_list), ' models)')
  res_env_breast  = run_envelope_set(breast_con_GEMs_list,  'Breast')
  message('Resident models, formula medium')
  res_env_formula = run_envelope_set(formula_con_GEMs_list, 'Formula')
  resident_envelope = bind_rows(res_env_breast, res_env_formula)
  saveRDS(resident_envelope, 'sensitivity/flux_envelope_residents.rds')
}

# ================================================================================================== #
# Verdicts on the four structural claims
# ================================================================================================== #

# Claim 1. "C. butyricum was the only butyrate producer among the seven probiotic strains."
# Certified if the capability maximum of EX_but(e) is zero, or the reaction is absent, for the other six.
claim_butyrate = pro_envelope %>%
  filter(reaction == 'EX_but(e)', mode == 'capability') %>%
  select(strain, medium, present, max_flux) %>%
  mutate(can_secrete = present & max_flux > 0) %>%
  arrange(medium, desc(can_secrete), strain)

# Claim 2. "None of the seven probiotic strains showed flux through the propionate exchange reaction."
claim_propionate = pro_envelope %>%
  filter(reaction == 'EX_ppa(e)', mode == 'capability') %>%
  select(strain, medium, present, min_flux, max_flux) %>%
  mutate(any_flux = present & (max_flux > 0 | min_flux < 0)) %>%
  arrange(medium, strain)

# Claim 3. "Only B. infantis and B. bifidum showed uptake fluxes for the eight HMO exchange reactions."
# Uptake is a negative flux, so the capability minimum is the bound that matters.
claim_hmo = pro_envelope %>%
  filter(reaction %in% unname(hmo_reactions), mode == 'capability') %>%
  group_by(strain, medium) %>%
  summarise(n_hmo_present = sum(present),
            n_hmo_can_take_up = sum(present & min_flux < 0, na.rm = TRUE),
            .groups = 'drop') %>%
  arrange(medium, desc(n_hmo_can_take_up))

# Claim 4. Which of C. butyricum's reported secretions are obligate rather than solution dependent.
claim_cbu_obligate = pro_envelope %>%
  filter(strain == 'C.butyricum', reaction %in% scfas) %>%
  select(sub, medium, mode, min_flux, max_flux) %>%
  pivot_wider(names_from = mode, values_from = c(min_flux, max_flux)) %>%
  mutate(verdict = case_when(
    max_flux_capability == 0                       ~ 'cannot secrete',
    min_flux_obligate  >  0                        ~ 'obligate secretion',
    max_flux_obligate  >  0                        ~ 'possible, solution dependent',
    TRUE                                           ~ 'no secretion at optimum'))

print(claim_butyrate,     n = 20)
print(claim_propionate,   n = 20)
print(claim_hmo,          n = 20)
print(claim_cbu_obligate, n = 20)

# -------------------------------------------------------------------------------------------------- #
# Degeneracy summary across all probiotic strain x reaction pairs. This is the number the manuscript
# needs in order to say how much of the attribution is structural. A pair is "determined" when the
# obligate interval has collapsed to a point, and "degenerate" when it has not.
# -------------------------------------------------------------------------------------------------- #
degeneracy_summary = pro_envelope %>%
  filter(mode == 'obligate', present, status == 'ok') %>%
  mutate(width = max_flux - min_flux,
         class = case_when(width < 1e-6            ~ 'determined',
                           min_flux > 0            ~ 'sign determined, magnitude free',
                           max_flux < 0            ~ 'sign determined, magnitude free',
                           TRUE                    ~ 'sign free')) %>%
  count(medium, class) %>%
  group_by(medium) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

print(degeneracy_summary)

# -------------------------------------------------------------------------------------------------- #
# Supplementary table. One row per strain x medium x compound, giving both envelopes.
# -------------------------------------------------------------------------------------------------- #
supp_table_S22 = pro_envelope %>%
  filter(mode == 'capability', status == 'ok' | !present) %>%
  select(Strain = strain, Medium = medium, Compound = sub, Reaction = reaction,
         `In model` = present, `Minimum flux` = min_flux, `Maximum flux` = max_flux) %>%
  arrange(Medium, Strain, Compound)

write.table(supp_table_S22, file = 'Supplementary_Table_S22.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

# ##################################################################################################################
# 7.2  Uniqueness of the minimised-total-flux solution under permutation of the reaction order
# ##################################################################################################################

# ================================================================================================== #
# Sensitivity analysis 1b: is the mtf solution unique?
#
# Why this is separate from sensitivity_01
#   The simulation solves with sec_obj = "mtf", so at each time step the flux distribution is the one
#   that minimises total flux among all distributions achieving the optimal biomass. That is a much
#   smaller solution set than the FBA optimal face, which is what the obligate mode of
#   sensitivity_01 bounds. Reporting the obligate interval as the degeneracy of the simulation would
#   therefore overstate it, sometimes by a lot.
#
#   The capability mode of sensitivity_01 is unaffected and remains the certificate for the
#   categorical claims: a maximum of zero with biomass free is a hard negative under any secondary
#   objective. This script handles the remaining question, which is whether the mtf flux values
#   themselves are pinned down.
#
# Method
#   An LP with a non-unique optimum has an optimal face with more than one vertex, and which vertex
#   the simplex method returns depends on the order in which it encounters the columns. Permuting the
#   reactions of a model produces a mathematically identical LP presented in a different column
#   order. If the mtf optimum is unique, every permutation returns the same exchange fluxes. If it is
#   not, permutations disagree, and the spread across permutations is a direct measure of how much
#   the reported flux could have differed.
#
#   Disagreement is proof of degeneracy. Agreement is strong evidence of uniqueness rather than proof,
#   since a finite set of permutations cannot enumerate every vertex. That distinction is worth
#   stating plainly in the manuscript rather than glossing.
#
# Prerequisites
#   Same workspace as sensitivity_01. Additionally requires that optimizeProb(m, algorithm = 'mtf')
#   runs, which it must, since the simulation depends on it.
#
# Runtime
#   n_perm x 2 media x ~156 models mtf solves. With n_perm = 6 that is about 1900 solves. Minutes.
#
# Output
#   sensitivity/mtf_uniqueness.rds
#   Supplementary_Table_S22b.txt
#   console verdict: the proportion of exchange fluxes that are uniquely determined
# ================================================================================================== #


dir.create('sensitivity', showWarnings = FALSE)

n_perm = 6              # first is the identity permutation, so five genuine re-orderings
perm_seed = 20260929
rel_tol = 1e-6          # two fluxes agree if they differ by less than this, relative to the larger

target_reactions = c(scfas, unname(hmo_reactions))

# -------------------------------------------------------------------------------------------------- #
# Reordering the columns of a modelorg. The stoichiometric matrix, the bounds, the objective and the
# reaction identifiers all have to move together. Fields that do not enter the LP are permuted too
# where they exist, so that the returned object stays internally consistent if inspected.
# -------------------------------------------------------------------------------------------------- #
permute_model = function(model, perm) {
  m = model
  sybil::S(m)          = sybil::S(model)[, perm, drop = FALSE]
  sybil::lowbnd(m)     = sybil::lowbnd(model)[perm]
  sybil::uppbnd(m)     = sybil::uppbnd(model)[perm]
  sybil::obj_coef(m)   = sybil::obj_coef(model)[perm]
  sybil::react_id(m)   = sybil::react_id(model)[perm]
  rn = try(sybil::react_name(model), silent = TRUE)
  if (!inherits(rn, 'try-error') && length(rn) == length(perm)) sybil::react_name(m) = rn[perm]
  rr = try(sybil::react_rev(model), silent = TRUE)
  if (!inherits(rr, 'try-error') && length(rr) == length(perm)) sybil::react_rev(m) = rr[perm]
  m
}

# -------------------------------------------------------------------------------------------------- #
# mtf flux of the target exchanges under one permutation. Fluxes are returned against the original
# reaction identifiers, so the permutation is invisible to the caller.
# -------------------------------------------------------------------------------------------------- #
mtf_exchange_flux = function(model, reactions, perm = NULL) {
  m = if (is.null(perm)) model else permute_model(model, perm)
  sol = try(sybil::optimizeProb(m, algorithm = 'mtf'), silent = TRUE)
  if (inherits(sol, 'try-error') || sol@lp_stat != 5) {
    return(tibble(reaction = reactions, flux = NA_real_, total_flux = NA_real_, ok = FALSE))
  }
  fd  = sybil::getFluxDist(sol)
  rid = sybil::react_id(m)
  idx = match(reactions, rid)
  tibble(reaction = reactions,
         flux = ifelse(is.na(idx), NA_real_, fd[idx]),
         total_flux = sol@lp_obj,
         ok = TRUE)
}

run_permutation_set = function(model_list, medium_label, reactions = target_reactions) {
  bind_rows(lapply(names(model_list), function(nm) {
    model = model_list[[nm]]
    nr = length(sybil::react_id(model))
    set.seed(perm_seed)
    perms = c(list(NULL), lapply(seq_len(n_perm - 1), function(i) sample.int(nr)))
    bind_rows(lapply(seq_along(perms), function(k) {
      mtf_exchange_flux(model, reactions, perms[[k]]) %>%
        mutate(perm_id = k, model = nm, medium = medium_label, .before = 1)
    }))
  }))
}

message('Probiotic strains, breast milk medium')
mtf_pro_breast  = run_permutation_set(PRO_cons_breast,  'Breast')
message('Probiotic strains, formula medium')
mtf_pro_formula = run_permutation_set(PRO_cons_formula, 'Formula')

mtf_pro = bind_rows(mtf_pro_breast, mtf_pro_formula) %>%
  mutate(strain = unname(pro_labels[model]))

run_residents = TRUE
if (run_residents) {
  message('Resident models, breast milk medium (', length(breast_con_GEMs_list), ' models)')
  mtf_res_breast  = run_permutation_set(breast_con_GEMs_list,  'Breast')
  message('Resident models, formula medium')
  mtf_res_formula = run_permutation_set(formula_con_GEMs_list, 'Formula')
  mtf_res = bind_rows(mtf_res_breast, mtf_res_formula)
} else {
  mtf_res = NULL
}

mtf_all = bind_rows(mtf_pro %>% mutate(set = 'Probiotic'),
                    if (!is.null(mtf_res)) mtf_res %>% mutate(set = 'Resident', strain = model))
saveRDS(mtf_all, 'sensitivity/mtf_uniqueness.rds')

# ================================================================================================== #
# Verdicts
# ================================================================================================== #

# -------------------------------------------------------------------------------------------------- #
# Sanity check first. Every permutation is the same LP, so the minimised total flux must agree across
# permutations to solver tolerance. If it does not, something in permute_model is wrong and nothing
# below can be trusted.
# -------------------------------------------------------------------------------------------------- #
total_flux_check = mtf_all %>%
  filter(ok) %>%
  distinct(set, medium, model, perm_id, total_flux) %>%
  group_by(set, medium, model) %>%
  summarise(spread = max(total_flux) - min(total_flux),
            rel_spread = spread / pmax(abs(max(total_flux)), 1e-12), .groups = 'drop')

cat('\n--- Sanity check: spread of the mtf objective across permutations ---\n')
print(summary(total_flux_check$rel_spread))
if (max(total_flux_check$rel_spread, na.rm = TRUE) > 1e-6) {
  warning('The mtf objective is not identical across permutations. ',
          'permute_model is not producing an equivalent LP; stop and fix before interpreting.')
}

# -------------------------------------------------------------------------------------------------- #
# Uniqueness per model x medium x reaction. A flux is treated as determined when its spread across
# permutations is negligible relative to its own magnitude, and separately when it is zero in every
# permutation, which is the case that matters most for the categorical claims.
# -------------------------------------------------------------------------------------------------- #
perm_ok = mtf_all %>%
  filter(ok) %>%
  distinct(set, medium, model, perm_id, total_flux) %>%
  group_by(set, medium, model) %>%
  mutate(ref = median(total_flux),
         keep = total_flux >= 0 &
           abs(total_flux - ref) <= pmax(1e-6, 1e-4 * abs(ref))) %>%
  ungroup() %>%
  select(set, medium, model, perm_id, keep)

mtf_all_clean = mtf_all %>%
  left_join(perm_ok, by = c('set','medium','model','perm_id')) %>%
  filter(is.na(keep) | keep)

cat('dropped permutations:', sum(!perm_ok$keep), 'of', nrow(perm_ok), '\n')

abs_tol = 1e-6          # comfortably above the GLPK feasibility tolerance of 1e-7

# (a stray reload of the uncleaned table was removed here; uniqueness is built from mtf_all_clean)

uniqueness = mtf_all_clean %>%
  filter(ok, !is.na(flux)) %>%
  group_by(set, strain, medium, reaction) %>%
  summarise(n_perm_ok = n(), f_min = min(flux), f_max = max(flux),
            f_ref = flux[perm_id == 1][1], .groups = 'drop') %>%
  mutate(across(c(f_min, f_max, f_ref), ~ ifelse(abs(.x) < abs_tol, 0, .x)),
         spread = f_max - f_min,
         scale = pmax(abs(f_min), abs(f_max)),
         rel_spread = ifelse(scale > 0, spread / scale, 0),
         sign_stable = sign(f_min) == sign(f_max),
         class = case_when(scale == 0 ~ 'always zero',
                           spread < abs_tol ~ 'determined',
                           rel_spread < 1e-6 ~ 'determined',
                           sign_stable ~ 'sign determined, magnitude varies',
                           TRUE ~ 'sign varies'),
         compound = case_when(
           reaction %in% names(scfa_flux_labels) ~ unname(scfa_flux_labels[reaction]),
           reaction %in% unname(hmo_reactions)   ~ names(hmo_reactions)[match(reaction, hmo_reactions)],
           TRUE ~ reaction))

uniqueness_summary = uniqueness %>% count(set, medium, class) %>%
  group_by(set, medium) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
print(uniqueness_summary, n = 40)
print(uniqueness %>% filter(class %in% c('sign determined, magnitude varies','sign varies')) %>%
        arrange(desc(rel_spread)) %>%
        select(set, strain, medium, compound, f_ref, f_min, f_max, rel_spread, class), n = 50)

# -------------------------------------------------------------------------------------------------- #
# The reactions that are not determined. These are the ones whose reported median net flux is
# solution dependent, and the manuscript should attribute them qualitatively rather than by
# magnitude. An empty table here is the best possible outcome.
# -------------------------------------------------------------------------------------------------- #
degenerate_exchanges = uniqueness %>%
  filter(class %in% c('sign determined, magnitude varies', 'sign varies')) %>%
  arrange(desc(rel_spread)) %>%
  select(set, strain, medium, compound, f_ref, f_min, f_max, rel_spread, class)

cat('\n--- Exchange fluxes that move between permutations ---\n')
print(degenerate_exchanges, n = 50)

# -------------------------------------------------------------------------------------------------- #
# The four claims again, now under mtf rather than over the FBA face.
# -------------------------------------------------------------------------------------------------- #
claim_butyrate_mtf = uniqueness %>%
  filter(set == 'Probiotic', compound == 'Butyrate') %>%
  select(strain, medium, f_ref, f_min, f_max, class) %>%
  arrange(medium, desc(f_ref))

claim_propionate_mtf = uniqueness %>%
  filter(set == 'Probiotic', compound == 'Propionate') %>%
  select(strain, medium, f_ref, f_min, f_max, class) %>%
  arrange(medium, strain)

claim_hmo_mtf = uniqueness %>%
  filter(set == 'Probiotic', reaction %in% unname(hmo_reactions)) %>%
  group_by(strain, medium) %>%
  summarise(n_uptake = sum(f_max < 0),
            n_any_nonzero = sum(class != 'always zero'),
            .groups = 'drop') %>%
  arrange(medium, desc(n_uptake))

cat('\n--- Butyrate, mtf flux per strain ---\n');   print(claim_butyrate_mtf, n = 20)
cat('\n--- Propionate, mtf flux per strain ---\n'); print(claim_propionate_mtf, n = 20)
cat('\n--- HMO uptake, mtf flux per strain ---\n'); print(claim_hmo_mtf, n = 20)

# -------------------------------------------------------------------------------------------------- #
# Supplementary table
# -------------------------------------------------------------------------------------------------- #
supp_table_S22b = uniqueness %>%
  filter(set == 'Probiotic') %>%
  select(Strain = strain, Medium = medium, Compound = compound, Reaction = reaction,
         `mtf flux` = f_ref, `Min across permutations` = f_min,
         `Max across permutations` = f_max, `Relative spread` = rel_spread,
         Classification = class) %>%
  arrange(Medium, Strain, Compound)

write.table(supp_table_S22b, file = 'Supplementary_Table_S22b.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

# ##################################################################################################################
# 7.3  Within-genus strain choice
# ##################################################################################################################

# ================================================================================================== #
# Sensitivity analysis 2, stage A: does the choice of strain within a genus change the repertoire?
#
# Purpose
#   Each genus-level ASV was assigned to the first unused strain of that genus in descending order of
#   abundance. For Bifidobacterium this matters more than for other genera, because the HMO conclusions
#   in sections 3, 4 and 6 depend on which oligosaccharide exchange reactions the assigned model
#   carries, and because the manuscript rests an internal control on B. lactis BB-12 lacking them.
#
#   Re-simulating under alternative assignments is expensive. This stage is the cheap screen that
#   decides whether it is necessary at all. It reads every Bifidobacterium model available in both
#   libraries and tabulates which of the eight quantified HMO exchanges and seven fermentation
#   product exchanges each one carries.
#
#   Decision rule for stage B (sensitivity_03):
#     If every candidate strain for a given ASV carries the same HMO exchange set, the assignment
#     cannot have changed any HMO result and no re-simulation is needed for that ASV.
#     If the candidates differ, that ASV enters stage B.
#
# Prerequisites
#   Run after Part 1 of major_script.R. Uses readMATmod_fix, hmo_reactions, scfas, hmo_dir,
#   ASV_GEM_map and the AGORA2 directory. Set agora_dir below to the AGORA2 .mat location.
#
# Output
#   sensitivity/bifido_repertoire.txt      every Bifidobacterium model x reaction
#   sensitivity/assignment_alternatives.txt per resident ASV, the candidates and whether they differ
#   console verdict on the BB-12 internal control
# ================================================================================================== #


dir.create('sensitivity', showWarnings = FALSE)

# agora_dir and hmo_dir are set in the path block of Part 0.

target_reactions = c(scfas, unname(hmo_reactions))

# -------------------------------------------------------------------------------------------------- #
# Reading a model only to ask which reactions it carries does not need the full sybil object, but
# readMATmod_fix is already validated against this library so it is reused rather than reimplemented.
# -------------------------------------------------------------------------------------------------- #
reaction_repertoire = function(path, reactions = target_reactions) {
  m = try(readMATmod_fix(path), silent = TRUE)
  if (inherits(m, 'try-error')) {
    return(tibble(reaction = reactions, present = NA, readable = FALSE))
  }
  tibble(reaction = reactions, present = reactions %in% sybil::react_id(m), readable = TRUE)
}

scan_library = function(dir_path, pattern, library_label) {
  files = list.files(dir_path, pattern = '\\.mat$', full.names = TRUE)
  files = files[grepl(pattern, basename(files), ignore.case = TRUE)]
  message(library_label, ': ', length(files), ' models matching "', pattern, '"')
  bind_rows(lapply(files, function(f) {
    reaction_repertoire(f) %>%
      mutate(model_file = basename(f),
             strain = sub('\\.mat$', '', basename(f)),
             library = library_label, .before = 1)
  }))
}

# -------------------------------------------------------------------------------------------------- #
# Stage A1. Every Bifidobacterium model in both libraries.
# -------------------------------------------------------------------------------------------------- #
bifido_hmo   = scan_library(hmo_dir,   '^Bifidobacterium', 'Shaaban HMO-extended')
bifido_agora = scan_library(agora_dir, '^Bifidobacterium', 'AGORA2 original')

bifido_repertoire = bind_rows(bifido_hmo, bifido_agora) %>%
  mutate(compound = case_when(
    reaction %in% names(scfa_flux_labels) ~ unname(scfa_flux_labels[reaction]),
    reaction %in% unname(hmo_reactions)   ~ names(hmo_reactions)[match(reaction, hmo_reactions)],
    TRUE ~ reaction),
    class = ifelse(reaction %in% unname(hmo_reactions), 'HMO', 'Fermentation product'))

bifido_wide = bifido_repertoire %>%
  filter(class == 'HMO') %>%
  select(library, strain, compound, present) %>%
  pivot_wider(names_from = compound, values_from = present) %>%
  mutate(n_hmo = rowSums(across(all_of(hmo_levels)), na.rm = TRUE)) %>%
  arrange(library, desc(n_hmo), strain)

write.table(bifido_wide, file = 'sensitivity/bifido_repertoire.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

print(bifido_wide, n = 100)

# -------------------------------------------------------------------------------------------------- #
# Stage A2. The BB-12 internal control.
#
# The manuscript argues that B. lactis BB-12 gained no colonisation advantage from HMO because it
# lacks the degradation pathways, in contrast with B. infantis and B. bifidum. That argument is
# sound only if the absence is a property of the organism rather than of the library. If no
# B. animalis reconstruction in either library carries any HMO exchange, the honest statement is
# that the library provides no HMO-extended B. animalis model, which is a limitation and must be
# written as one.
# -------------------------------------------------------------------------------------------------- #
animalis_check = bifido_wide %>%
  filter(grepl('animalis|lactis', strain, ignore.case = TRUE))

cat('\n--- B. animalis / B. lactis models and their HMO exchange counts ---\n')
print(animalis_check, n = 50)

if (nrow(animalis_check) == 0) {
  cat('\nVERDICT: no B. animalis model found in either library under this pattern.\n',
      'Check the file naming before concluding.\n')
} else if (all(animalis_check$n_hmo == 0)) {
  cat('\nVERDICT: no B. animalis reconstruction in either library carries any of the eight HMO\n',
      'exchanges. The BB-12 contrast is therefore a property of the model library, not a\n',
      'demonstrated property of the strain. State this in the Discussion as a limitation.\n')
} else {
  cat('\nVERDICT: at least one B. animalis reconstruction carries HMO exchanges. The BB-12\n',
      'result is strain specific and the manuscript must say which reconstruction was used\n',
      'and why, since an alternative choice would have changed the internal control.\n')
}

# -------------------------------------------------------------------------------------------------- #
# Stage A3. For each resident Bifidobacterium ASV, which strains were candidates at assignment time,
# and would any of them have given a different HMO repertoire?
#
# The assignment rule was: within a genus, take strains in descending abundance order and give each
# ASV the first strain not already taken. The candidate set for an ASV is therefore every strain of
# that genus present in the library. This reconstructs that set and asks whether the repertoire is
# invariant across it.
# -------------------------------------------------------------------------------------------------- #
resident_bifido = ASV_GEM_map %>%
  filter(grepl('Bifidobacterium', Genus, ignore.case = TRUE)) %>%
  select(ASV_ID, Genus, strain, renamed_gem)

cat('\n--- Resident Bifidobacterium ASVs and their assigned strains ---\n')
print(head(resident_bifido, 50))

# Repertoire of the strains actually assigned, taken from whichever library supplied each model.
assigned_repertoire = bifido_wide %>%
  filter(strain %in% resident_bifido$strain)

candidate_spread = bifido_wide %>%
  group_by(library) %>%
  summarise(n_models = n(),
            n_hmo_min = min(n_hmo), n_hmo_max = max(n_hmo),
            n_distinct_repertoires = n_distinct(do.call(paste, across(all_of(hmo_levels)))),
            .groups = 'drop')

cat('\n--- Spread of HMO repertoires among candidate Bifidobacterium strains ---\n')
print(candidate_spread)

if (nrow(assigned_repertoire) > 0 &&
    all(candidate_spread$n_distinct_repertoires == 1)) {
  cat('\nVERDICT: every candidate Bifidobacterium strain carries the same HMO exchange set.\n',
      'Strain choice within the genus cannot have changed any HMO result and stage B is not\n',
      'required for Bifidobacterium.\n')
} else {
  cat('\nVERDICT: candidate strains differ in their HMO exchange sets. Stage B re-simulation is\n',
      'required. Use sensitivity_03 with the alternative assignments written below.\n')
}

write.table(resident_bifido %>% left_join(bifido_wide, by = 'strain'),
            file = 'sensitivity/assignment_alternatives.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

# -------------------------------------------------------------------------------------------------- #
# Stage A4. The same screen for Bacteroidetes, which carries the other half of the HMO division of
# labour reported in section 3. Run this only if A3 shows that strain choice matters, since the
# Bacteroides genus has many more models and the scan is correspondingly slower.
# -------------------------------------------------------------------------------------------------- #
run_bacteroides = FALSE

if (run_bacteroides) {
  bacteroides_hmo = scan_library(hmo_dir, '^Bacteroides', 'Shaaban HMO-extended')
  bacteroides_wide = bacteroides_hmo %>%
    filter(reaction %in% unname(hmo_reactions)) %>%
    mutate(compound = names(hmo_reactions)[match(reaction, hmo_reactions)]) %>%
    select(strain, compound, present) %>%
    pivot_wider(names_from = compound, values_from = present) %>%
    mutate(n_hmo = rowSums(across(all_of(hmo_levels)), na.rm = TRUE)) %>%
    arrange(desc(n_hmo), strain)
  write.table(bacteroides_wide, file = 'sensitivity/bacteroides_repertoire.txt',
              sep = '\t', quote = FALSE, row.names = FALSE)
  print(bacteroides_wide, n = 100)
}

# ##################################################################################################################
# 7.4  Supplementary Figure S10
# ##################################################################################################################

# ================================================================================================== #
# Supplementary Figure S15: how far the exchange fluxes are determined by the models
#
# What it supports
#   The paragraph at the end of the second Results section reports six percentages, a breakdown of the
#   62 fluxes whose sign was not determined, and the mechanism behind the largest group of exceptions.
#   None of that is checkable from the text alone. This figure carries all three in one place, which
#   suits the material better than a table: the claims are about how one quantity moves between
#   solutions, and movement is easier to see than to read off rows.
#
#   (A) Composition of the four classes, by model set and medium. Supports the six percentages and the
#       statement that no flux of an introduced strain changed sign.
#   (B) The fluxes whose sign was not determined, by compound. Supports the 48 and the further 10.
#   (C) The mechanism for the lactate exchanges in four of the affected strains: across permutations
#       the two enantiomers exchange the flux while their sum stays flat.
#
# Prerequisites
#   uniqueness and mtf_all_clean in the workspace, both after the absolute-tolerance correction and
#   the removal of the four failed permutations.
#
# Output
#   Figures/figure_S15_uniqueness.jpeg
# ================================================================================================== #


class_levels = c('Always zero', 'Determined',
                 'Sign determined,\nmagnitude varies', 'Sign not determined')

# The first two classes carry no problem for the attribution and are drawn in greys; the last two are
# the ones the text qualifies, so they take the two warm colours.
class_colours = c('Always zero'                        = 'grey85',
                  'Determined'                         = 'grey55',
                  'Sign determined,\nmagnitude varies' = '#fdae61',
                  'Sign not determined'                = '#d73027')

relabel_class = function(x) {
  factor(recode(x,
                'always zero'                       = 'Always zero',
                'determined'                        = 'Determined',
                'sign determined, magnitude varies' = 'Sign determined,\nmagnitude varies',
                'sign varies'                       = 'Sign not determined'),
         levels = class_levels)
}

strip_model_name = function(x) sub('^ASV[0-9]+_', '', sub('\\.mat$', '', x))

# -------------------------------------------------------------------------------------------------- #
# Panel A. Class composition
# -------------------------------------------------------------------------------------------------- #
panelA_dat = uniqueness %>%
  mutate(Set = factor(ifelse(set == 'Probiotic', 'Introduced strains', 'Resident models'),
                      levels = c('Resident models', 'Introduced strains')),
         Class = relabel_class(class)) %>%
  count(Set, medium, Class) %>%
  group_by(Set, medium) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

panelA = ggplot(panelA_dat, aes(x = medium, y = pct, fill = Class))+
  geom_col(width = 0.65)+
  facet_wrap(~ Set)+
  scale_fill_manual(values = class_colours, name = NULL)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)))+
  theme_bw()+
  labs(x = "", y = "Exchange fluxes (%)")+
  theme(axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.height = unit(1.1, 'lines'),
        panel.grid.major.x = element_blank())

# -------------------------------------------------------------------------------------------------- #
# Panel B. Where the undetermined signs sit
# -------------------------------------------------------------------------------------------------- #
panelB_dat = uniqueness %>%
  filter(class == 'sign varies') %>%
  count(medium, compound) %>%
  mutate(compound = reorder(compound, n, sum))

panelB = ggplot(panelB_dat, aes(x = n, y = compound, fill = medium))+
  geom_col(width = 0.7)+
  scale_fill_manual(values = diet_colors, name = "Medium")+
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_bw()+
  labs(x = "Fluxes whose sign was not determined", y = "")+
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = 'bottom',
        panel.grid.major.y = element_blank())

# -------------------------------------------------------------------------------------------------- #
# Panel C. The lactate exchange mechanism
#
# The four strains named in the Results text. Each panel shows the two enantiomer fluxes across the
# six permutations together with their sum: the sum is flat while the two components exchange it,
# which is what makes the enantiomer assignment undetermined and the total reliable.
# -------------------------------------------------------------------------------------------------- #
focus_strains = c('ASV006_Staphylococcus_sp_HMSC070D05.mat',
                  'ASV035_Escherichia_coli_SE11.mat',
                  'ASV066_Lactobacillus_rhamnosus_ATCC_21052.mat',
                  'ASV024_Lactobacillus_acidophilus_ATCC_4796.mat')

# CHECK: if a name does not match, list the candidates with
#   mtf_all_clean %>% distinct(strain) %>% filter(grepl('Lactobacillus', strain))
stopifnot(all(focus_strains %in% unique(mtf_all_clean$strain)))

lac_long = mtf_all_clean %>%
  filter(ok, strain %in% focus_strains,
         reaction %in% c('EX_lac_D(e)', 'EX_lac_L(e)')) %>%
  mutate(Series = ifelse(reaction == 'EX_lac_D(e)', '(R)-Lactate', '(S)-Lactate'),
         Strain = strip_model_name(strain))

lac_total = lac_long %>%
  group_by(Strain, medium, perm_id) %>%
  summarise(flux = sum(flux), .groups = 'drop') %>%
  mutate(Series = 'Sum')

panelC_dat = bind_rows(lac_long %>% select(Strain, medium, perm_id, Series, flux), lac_total) %>%
  mutate(Series = factor(Series, levels = c('(R)-Lactate', '(S)-Lactate', 'Sum')))

series_colours = c('(R)-Lactate' = '#7b3294', '(S)-Lactate' = '#e66101', 'Sum' = 'grey20')

panelC = ggplot(panelC_dat, aes(x = factor(perm_id), y = flux,
                                colour = Series, group = Series, linetype = Series))+
  geom_line(linewidth = 0.7)+
  geom_point(size = 1.6)+
  facet_grid(medium ~ Strain, scales = 'free_y')+
  scale_colour_manual(values = series_colours, name = NULL)+
  scale_linetype_manual(values = c('(R)-Lactate' = 'solid', '(S)-Lactate' = 'solid',
                                   'Sum' = 'dashed'), name = NULL)+
  theme_bw()+
  labs(x = "Permutation of the reaction order", y = "Exchange flux")+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 8, face = 'italic'),
        strip.text.y = element_text(size = 10),
        legend.position = 'bottom')

# -------------------------------------------------------------------------------------------------- #
# Assembly
# -------------------------------------------------------------------------------------------------- #
figure_S15 = (panelA | panelB) / panelC +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = 'A')

ggsave(plot = figure_S15, filename = 'Figures/figure_S15_uniqueness.jpeg',
       dpi = 1200, width = 26, height = 20, units = 'cm')

# ================================================================================================================== #
# Part 8  Personalised nutrient supplementation
#
#   Manuscript: Results 5; Figure 4A to 4C; Tables S17 to S19
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Prediction of personalised nutrient supplements. For each preterm infant, the strains that secrete a short-chain
# fatty acid are examined in turn: the strain is placed on that individual's medium and the reduced costs of its
# exchange reactions are read from the dual solution, so that the substrates limiting its growth are identified
# whether or not they are already present in the medium. Relieving the limiting exchanges and re-solving exposes the
# constraints that only become binding once the first ones are lifted; three rounds are recorded but only the first
# is used downstream, because later rounds return interchangeable alternatives rather than genuine constraints.
# Metabolites the community already secretes, and oxygen, are excluded since neither represents a dietary addition.
# The prediction itself runs as a SLURM array over the 40 preterm infants (predict_supplements.R); the code below
# reads and summarises those results.
# ---------------------------------------------------------------------------------------------------------------- #

read_supplements = function(dir_path, max_round = supplement_round) {
  f = list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  ids = sub("_needed_mets\\.rds$", "", basename(f))
  do.call(rbind, lapply(seq_along(f), function(i) {
    d = readRDS(f[i])
    d = d[d$round <= max_round & !d$react %in% exclude_mets, , drop = FALSE]
    if (nrow(d) > 0) d$ID = ids[i]
    d
  }))
}

breast_supplements = read_supplements("supplement_pred_breast") %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = c('S-Preterm','N-Preterm')), Diet = 'Breast')
formula_supplements = read_supplements("supplement_pred_formula") %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = c('S-Preterm','N-Preterm')), Diet = 'Formula')

all_supplements = rbind(breast_supplements, formula_supplements) %>%
  mutate(Diet = factor(Diet, levels = c('Breast','Formula')))

# Per-individual supplement: the union of the substrates limiting that individual's SCFA producers
breast_mets_used = lapply(split(breast_supplements$react, breast_supplements$ID), unique)
formula_mets_used = lapply(split(formula_supplements$react, formula_supplements$ID), unique)
#saveRDS(breast_mets_used, 'breast_supplement_mets.rds')
#saveRDS(formula_mets_used, 'formula_supplement_mets.rds')
breast_mets_used = readRDS(file.path(data_dir, "derived/breast_supplement_mets.rds"))
formula_mets_used = readRDS(file.path(data_dir, "derived/formula_supplement_mets.rds"))

# ---- Summaries ---------------------------------------------------------------------------------------------------
supplement_per_subject = all_supplements %>%
  group_by(Diet, ID, group) %>%
  summarise(n_mets = n_distinct(react),
            n_producers = n_distinct(model), .groups = 'drop')

supplement_per_subject_summary = supplement_per_subject %>%
  group_by(Diet, group) %>%
  summarise(n = n(), median_mets = median(n_mets),
            q25 = quantile(n_mets, .25), q75 = quantile(n_mets, .75),
            max_mets = max(n_mets), .groups = 'drop')

# How widely is each substrate limiting? Counted both as strain-by-subject events and as subjects
supplement_freq = all_supplements %>%
  group_by(Diet, react) %>%
  summarise(n_events = n(),
            n_subjects = n_distinct(ID),
            pct_subjects = 100 * n_subjects / n_distinct(all_supplements$ID),
            .groups = 'drop') %>%
  arrange(Diet, desc(n_subjects))

# Overlap between the two feeding backgrounds within the same subject
supplement_overlap = data.frame(
  ID = names(breast_mets_used),
  n_breast = lengths(breast_mets_used),
  n_formula = lengths(formula_mets_used[names(breast_mets_used)]),
  n_shared = mapply(function(a, b) length(intersect(a, b)),
                    breast_mets_used, formula_mets_used[names(breast_mets_used)]),
  row.names = NULL) %>%
  mutate(jaccard = n_shared / (n_breast + n_formula - n_shared)) %>%
  left_join(meta %>% select(ID, group), by = 'ID')

# ---- Map reaction identifiers to metabolite names -----------------------------------------------------------------
# Names are taken from the arena media vector, which carries the readable label for every exchange reaction
met_name_lookup = function(reac_ids, mediac) {
  vapply(reac_ids, function(rx) {
    hit = which(mediac == rx)
    if (length(hit) > 0) names(mediac)[hit[1]] else NA_character_
  }, character(1))
}

supplement_table = supplement_freq %>%
  mutate(Met_name = met_name_lookup(react, breast_ext$mediac)) %>%
  select(Diet, Reac_ID = react, Met_name, n_subjects, pct_subjects, n_events)

write.table(supplement_table, file = "supplement_table.txt",
            row.names = FALSE, sep = '\t', quote = FALSE)

# ---- Plots -------------------------------------------------------------------------------------------------------

supplement_n_plot = ggplot(supplement_per_subject, aes(x = Diet, y = n_mets, fill = Diet))+
  geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA)+
  geom_point(aes(colour = group), position = position_jitter(width = 0.12),
             size = 1.5, alpha = 0.7)+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  scale_colour_manual(values = group_colors, name = "Group")+
  theme_bw()+
  labs(title = "Number of limiting substrates per individual",
       x = "", y = "Substrates")+
  theme(axis.text = element_text(size = 12), title = element_text(size = 12))

top_mets = supplement_freq %>%
  group_by(react) %>% summarise(total = sum(n_subjects), .groups = 'drop') %>%
  arrange(desc(total)) %>% head(20) %>% pull(react)

supplement_freq_plot = supplement_freq %>%
  filter(react %in% top_mets) %>%
  mutate(react = factor(react, levels = rev(top_mets))) %>%
  ggplot(aes(x = react, y = pct_subjects, fill = Diet))+
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8)+
  coord_flip()+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  theme_bw()+
  labs(title = "Substrates limiting the fermentation-product producers in the preterm infants",
       x = "", y = "Individuals affected (%)")+
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 11), title = element_text(size = 12))


# ---------------------------------------------------------------------------------------------------------------- #
# Effect of the personalised nutrient supplementation on predicted SCFA production. Preterm infants only. For each
# feeding background the supplemented and unsupplemented simulations of the same subject are compared with a paired
# Wilcoxon test, and the effect is contrasted between the two backgrounds. Individuals in which a given acid rises
# above the response threshold are classified as responders for that acid, so that the breadth of the response can
# be read per subject rather than only as a group median. Because the supplements were predicted from the SCFA
# producers of the untreated community, the flux attribution below establishes whether any gain comes from those
# producers or from other members of the community.
# ---------------------------------------------------------------------------------------------------------------- #

breast_supp_ext  = load_extracted('supp_breast', expected = 40)
formula_supp_ext = load_extracted('supp_formula', expected = 40)

# ---- SCFA before and after supplementation -----------------------------------------------------------------------
breast_supp_scfa = extract_scfa(breast_supp_ext, meta, 'Breast') %>%
  mutate(Status = 'Supplemented')
formula_supp_scfa = extract_scfa(formula_supp_ext, meta, 'Formula') %>%
  mutate(Status = 'Supplemented')

all_scfa_supp = rbind(
  breast_scfa  %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  formula_scfa %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  breast_supp_scfa, formula_supp_scfa) %>%
  mutate(Status = factor(Status, levels = c('Non-treated','Supplemented')),
         Diet = factor(Diet, levels = c('Breast','Formula')),
         group = factor(group, levels = c('S-Preterm','N-Preterm')),
         sub = factor(sub, levels = scfa_levels),
         ID = factor(ID, levels = unique(ID[order(group)])))

supp_effect_stat = all_scfa_supp %>%
  arrange(Diet, sub, Status, ID) %>%
  group_by(Diet, sub) %>%
  rstatix::wilcox_test(value ~ Status, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  group_by(Diet) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  ungroup()

supp_paired = all_scfa_supp %>%
  select(ID, group, Diet, sub, Status, value) %>%
  pivot_wider(names_from = Status, values_from = value) %>%
  mutate(Diff = Supplemented - `Non-treated`)

supp_effect_summary = supp_paired %>%
  group_by(Diet, sub) %>%
  summarise(n = sum(!is.na(Diff)),
            median_untreated = median(`Non-treated`, na.rm = TRUE),
            median_treated = median(Supplemented, na.rm = TRUE),
            median_diff = median(Diff, na.rm = TRUE),
            n_increased = sum(Diff > 0, na.rm = TRUE), .groups = 'drop') %>%
  left_join(supp_effect_stat %>% select(Diet, sub, p, p.adj, p.adj.signif),
            by = c('Diet','sub'))
write.table(supp_effect_summary,"Supplementary_Table_18.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

supp_interaction_stat = supp_paired %>%
  arrange(sub, Diet, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(Diff ~ Diet, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

# ---- Responders ---------------------------------------------------------------------------------------------------
supp_response = supp_paired %>%
  mutate(delta = ifelse(Diff < resp_cutoff, 0, Diff),
         present = ifelse(delta > resp_cutoff, "Present", "Absent"))

supp_response_summary = supp_response %>%
  group_by(Diet, sub) %>%
  summarise(n = n(),
            n_responded = sum(present == "Present"),
            pct_responded = 100 * n_responded / n,
            median_delta = median(delta), .groups = 'drop')

# Does the number of predicted supplements predict the size of the response?
supp_dose_response = supp_paired %>%
  left_join(supplement_per_subject %>% select(Diet, ID, n_mets), by = c('Diet','ID')) %>%
  group_by(Diet, sub) %>%
  summarise(rho = cor(n_mets, Diff, method = 'spearman', use = 'complete.obs'),
            p = cor.test(n_mets, Diff, method = 'spearman')$p.value,
            .groups = 'drop') %>%
  group_by(Diet) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  ungroup()

# ---- Which taxa respond to the added substrates? ------------------------------------------------------------------
# Net flux per taxon before and after supplementation; a positive change marks a taxon whose SCFA output rose
net_flux_table = function(extracted, diet_type, status) {
  extracted$flux %>%
    group_by(ID, spec, rea) %>%
    summarise(net_flux = sum(mflux, na.rm = TRUE), .groups = 'drop') %>%
    mutate(Diet = diet_type, Status = status)
}

supp_producer_shift = rbind(
  net_flux_table(breast_ext,       'Breast',  'Non-treated') %>% filter(ID %in% preterm_ids),
  net_flux_table(formula_ext,      'Formula', 'Non-treated') %>% filter(ID %in% preterm_ids),
  net_flux_table(breast_supp_ext,  'Breast',  'Supplemented'),
  net_flux_table(formula_supp_ext, 'Formula', 'Supplemented')) %>%
  mutate(rea = ifelse(rea %in% names(scfa_flux_labels), scfa_flux_labels[rea], rea)) %>%
  pivot_wider(names_from = Status, values_from = net_flux, values_fill = 0) %>%
  mutate(delta = Supplemented - `Non-treated`)

supp_responder_taxa = supp_producer_shift %>%
  filter(delta > 0) %>%
  group_by(Diet, rea, spec) %>%
  summarise(n_subjects = n_distinct(ID),
            median_delta = median(delta),
            total_delta = sum(delta), .groups = 'drop') %>%
  arrange(Diet, rea, desc(total_delta))

# ---- Per-individual comparison with a companion response panel ----------------------------------------------------
scfa_supp_overlay_plot = function(scfa_df, diet_type) {
  d = scfa_df %>% filter(Diet == diet_type)
  ggplot()+
    geom_bar(data = d %>% filter(Status == 'Supplemented'),
             aes(x = value, y = ID, fill = interaction(group, Status)),
             stat = 'identity', position = 'identity', alpha = 0.7)+
    geom_bar(data = d %>% filter(Status == 'Non-treated'),
             aes(x = value, y = ID, fill = interaction(group, Status)),
             stat = 'identity', position = 'identity', alpha = 1)+
    facet_grid(. ~ sub, scales = 'free_x', space = 'fixed')+
    scale_fill_manual(
      values = c("N-Preterm.Non-treated"  = "orange",
                 "S-Preterm.Non-treated"  = "red",
                 "N-Preterm.Supplemented" = "#FFCC99",
                 "S-Preterm.Supplemented" = "#F08080"),
      name = "Group",
      labels = c("N-Preterm.Non-treated"  = "N-Preterm (non-treated)",
                 "S-Preterm.Non-treated"  = "S-Preterm (non-treated)",
                 "N-Preterm.Supplemented" = "N-Preterm (supplemented)",
                 "S-Preterm.Supplemented" = "S-Preterm (supplemented)"))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 10),
          strip.text.y = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_markdown(size = 7),
          axis.text.x = element_text(angle = 60, vjust = 0.6, size = 8),
          panel.spacing = unit(0.2, "lines"),
          panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6),
          plot.margin = unit(c(0, -2000, 0, -2000), "pt"))+
    labs(x = "Fermentation product, log(mM+1)", y = "")+
    scale_y_discrete(position = 'right')
}

supp_resp_tile_plot = function(resp_df, diet_type) {
  ggplot(resp_df %>% filter(Diet == diet_type),
         aes(x = sub, y = ID, fill = present))+
    geom_tile(color = "white", width = 1.1, height = 0.9)+
    scale_fill_manual(values = c("Present" = "forestgreen", "Absent" = "lightgrey"),
                      labels = c("Present" = "Responded", "Absent" = "Non-responded"))+
    coord_fixed()+
    theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0, size = 6),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0, -2000, 0, -2000), "pt"))+
    scale_x_discrete(position = 'top')
}

breast_supp_scfa_resp = scfa_supp_overlay_plot(all_scfa_supp, 'Breast') +
  supp_resp_tile_plot(supp_response, 'Breast') +
  plot_layout(guides = 'collect')
formula_supp_scfa_resp = scfa_supp_overlay_plot(all_scfa_supp, 'Formula') +
  supp_resp_tile_plot(supp_response, 'Formula') +
  plot_layout(guides = 'collect')

ggsave(plot = breast_supp_scfa_resp,
       filename = 'Figures/Supp_treat_scfa_resp_breast.jpeg',
       dpi = 1200, width = 24, height = 12, units = 'cm')
ggsave(plot = formula_supp_scfa_resp,
       filename = 'Figures/Supp_treat_scfa_resp_formula.jpeg',
       dpi = 1200, width = 24, height = 12, units = 'cm')

# ---- Group-level comparison ---------------------------------------------------------------------------------------
supp_paired_plot = ggplot(all_scfa_supp, aes(x = Status, y = value))+
  geom_line(aes(group = ID, color = group), alpha = 0.3, linewidth = 0.4)+
  geom_point(aes(color = group), size = 1, alpha = 0.6)+
  geom_boxplot(width = 0.35, alpha = 0.25, outlier.shape = NA, fill = 'grey70')+
  facet_grid(Diet ~ sub, scales = 'free_y')+
  scale_color_manual(values = group_colors, name = "Group")+
  theme_bw()+
  labs(title = "Predicted fermentation products before and after nutrient supplementation",
       x = "", y = "Fermentation product, log(mM+1)")+
  theme(strip.text = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 9), title = element_text(size = 12))

supp_diff_plot = ggplot(supp_paired, aes(x = sub, y = Diff, fill = Diet))+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40')+
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1,
               position = position_dodge(width = 0.75))+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  theme_bw()+
  labs(title = "Change in predicted fermentation products after nutrient supplementation",
       x = "", y = "Supplemented minus untreated, log(mM+1)")+
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13), title = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.position = "bottom")
ggsave(plot = supp_diff_plot,filename = 'Figures/supp_diff_plot.jpeg',
       dpi = 1200,height = 13.4,width = 13.8,units = 'cm')

# ================================================================================================================== #
# Part 9  Combined probiotic and nutrient intervention
#
#   Manuscript: Results 5; Figure 4D; Table S20
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# Combined intervention. The seven probiotic strains are introduced at an initial abundance of 10 per strain into a
# medium carrying that individual's predicted nutrient supplements, so that the combined arm is the sum of the two
# single-agent arms rather than a separate design. Simulations run on the cluster through simulate_arena_syn.R; the
# analysis below reads their output and compares all four conditions.
# ---------------------------------------------------------------------------------------------------------------- #
breast_syn_ext  = load_extracted('syn_breast', expected = 40)
formula_syn_ext = load_extracted('syn_formula', expected = 40)

breast_syn_scfa = extract_scfa(breast_syn_ext, meta, 'Breast') %>%
  mutate(Status = 'Comb-treated')
formula_syn_scfa = extract_scfa(formula_syn_ext, meta, 'Formula') %>%
  mutate(Status = 'Comb-treated')

all_scfa_syn = rbind(
  breast_scfa  %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  formula_scfa %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  breast_syn_scfa, formula_syn_scfa) %>%
  mutate(Status = factor(Status, levels = c('Non-treated','Comb-treated')),
         Diet = factor(Diet, levels = c('Breast','Formula')),
         group = factor(group, levels = preterm_levels),
         sub = factor(sub, levels = scfa_levels),
         ID = factor(ID, levels = unique(ID[order(group)])))

syn_effect_stat = all_scfa_syn %>%
  arrange(Diet, sub, Status, ID) %>%
  group_by(Diet, sub) %>%
  rstatix::wilcox_test(value ~ Status, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  group_by(Diet) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  ungroup()

syn_paired = all_scfa_syn %>%
  select(ID, group, Diet, sub, Status, value) %>%
  pivot_wider(names_from = Status, values_from = value) %>%
  mutate(Diff = `Comb-treated` - `Non-treated`)

syn_effect_summary = syn_paired %>%
  group_by(Diet, sub) %>%
  summarise(n = sum(!is.na(Diff)),
            median_untreated = median(`Non-treated`, na.rm = TRUE),
            median_treated = median(`Comb-treated`, na.rm = TRUE),
            median_diff = median(Diff, na.rm = TRUE),
            n_increased = sum(Diff > 0, na.rm = TRUE), .groups = 'drop') %>%
  left_join(syn_effect_stat %>% select(Diet, sub, p, p.adj, p.adj.signif),
            by = c('Diet','sub'))

# Engraftment of the probiotic strains in the combined arm, for comparison with the probiotic-only arm
syn_engraft = rbind(engraftment(breast_syn_ext, 'Breast', pro_names),
                    engraftment(formula_syn_ext, 'Formula', pro_names)) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = preterm_levels),
         Diet = factor(Diet, levels = c('Breast','Formula')),
         species = factor(species, levels = pro_names))

# ---------------------------------------------------------------------------------------------------------------- #
# Comparison of the four conditions. The three interventions are contrasted against the untreated baseline and
# against each other, and the combined arm is tested for synergy by asking whether its effect exceeds the sum of
# the two single-agent effects.
# ---------------------------------------------------------------------------------------------------------------- #
scfa_compare_df = rbind(
  breast_scfa  %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  formula_scfa %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Non-treated'),
  breast_supp_scfa, formula_supp_scfa,
  breast_pro_scfa, formula_pro_scfa,
  breast_syn_scfa, formula_syn_scfa) %>%
  mutate(Status = factor(Status, levels = c('Non-treated','Supplemented',
                                            'Prob-treated','Comb-treated')),
         Diet = factor(Diet, levels = c('Breast','Formula')),
         group = factor(group, levels = preterm_levels),
         sub = factor(sub, levels = scfa_levels))

treat_colors = c("Non-treated" = "grey70", "Supplemented" = "#E7B800",
                 "Prob-treated" = "#00AFBB", "Comb-treated" = "#b2182b")

scfa_compare_stat = scfa_compare_df %>%
  arrange(Diet, sub, Status, ID) %>%
  group_by(Diet, sub) %>%
  rstatix::wilcox_test(value ~ Status, paired = TRUE, p.adjust.method = 'none',
                       ref.group = 'Non-treated') %>%
  ungroup() %>%
  group_by(Diet) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  ungroup()

# Effect of each intervention relative to the untreated baseline
all_effects = scfa_compare_df %>%
  select(ID, group, Diet, sub, Status, value) %>%
  pivot_wider(names_from = Status, values_from = value) %>%
  mutate(Nutrient = Supplemented - `Non-treated`,
         Probiotic = `Prob-treated` - `Non-treated`,
         Combined = `Comb-treated` - `Non-treated`,
         Expected = Nutrient + Probiotic,
         Synergy = Combined - Expected)

synergy_stat = all_effects %>%
  arrange(Diet, sub, ID) %>%
  group_by(Diet, sub) %>%
  rstatix::wilcox_test(Synergy ~ 1, mu = 0) %>%
  ungroup() %>%
  group_by(Diet) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj') %>%
  ungroup()

synergy_summary = all_effects %>%
  group_by(Diet, sub) %>%
  summarise(median_nutrient = median(Nutrient, na.rm = TRUE),
            median_probiotic = median(Probiotic, na.rm = TRUE),
            median_combined = median(Combined, na.rm = TRUE),
            median_expected = median(Expected, na.rm = TRUE),
            median_synergy = median(Synergy, na.rm = TRUE),
            n_synergistic = sum(Synergy > 0, na.rm = TRUE), .groups = 'drop') %>%
  left_join(synergy_stat %>% select(Diet, sub, p.adj, p.adj.signif), by = c('Diet','sub'))
write.table(synergy_summary,"Supplementary_Table_S19.txt",
            sep = '\t',row.names = FALSE,quote = FALSE)

# ---- Plots -------------------------------------------------------------------------------------------------------
scfa_treat_compare_plot = ggplot(scfa_compare_df, aes(x = Status, y = value, fill = Status))+
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6)+
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.35, colour = 'grey25')+
  facet_grid(Diet ~ sub, scales = 'free_y')+
  scale_fill_manual(values = treat_colors, name = "Condition")+
  theme_bw()+
  labs(x = "", y = "Fermentation product, log(mM+1)")+
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8))

effect_long = all_effects %>%
  select(ID, group, Diet, sub, Nutrient, Probiotic, Combined) %>%
  pivot_longer(cols = c(Nutrient, Probiotic, Combined),
               names_to = 'Intervention', values_to = 'Diff') %>%
  mutate(Intervention = factor(Intervention,
                               levels = c('Nutrient','Probiotic','Combined')))

intervention_effect_plot = ggplot(effect_long, aes(x = Intervention, y = Diff, fill = Intervention))+
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40')+
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.size = 0.8)+
  facet_grid(Diet ~ sub, scales = 'free_y')+
  scale_fill_manual(values = c("Nutrient" = "#E7B800", "Probiotic" = "#00AFBB",
                               "Combined" = "#b2182b"), name = "Intervention")+
  theme_bw()+
  labs(x = "", y = "Change from untreated, log(mM+1)")+
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8))

# Direction and magnitude per individual, so that decreases are visible alongside increases
effect_heatmap = ggplot(effect_long, aes(x = sub, y = ID, fill = Diff))+
  geom_tile(colour = "white")+
  facet_grid(Diet ~ Intervention)+
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, name = "Change in\nlog(mM+1)")+
  theme_minimal()+
  labs(x = "", y = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 5),
        panel.grid = element_blank())

# ================================================================================================================== #
# Part 10  Prebiotic and synbiotic intervention on the formula background
#
#   Manuscript: Results 6; Figure 4E to 4H; Tables S21, S22
# ================================================================================================================== #


# ---------------------------------------------------------------------------------------------------------------- #
# The formula medium carries no oligosaccharide, so transferring the thirteen HMO exchanges of the breast medium
# onto it isolates the prebiotic. Crossing that with the presence of the seven probiotic strains gives a 2 x 2
# design whose interaction term asks whether a prebiotic and the strains able to use it act together for more than
# the sum of their separate effects. The medium is built by prepare_hmo_formula_diet.R and the simulations run on
# the cluster through simulate_arena_hmo.R.
# ---------------------------------------------------------------------------------------------------------------- #
hmo_ext     = load_extracted('hmo_formula', expected = 40)
hmo_pro_ext = load_extracted('hmo_pro_formula', expected = 40)

hmo_scfa = extract_scfa(hmo_ext, meta, 'Formula') %>% mutate(Status = 'HMO')
hmo_pro_scfa = extract_scfa(hmo_pro_ext, meta, 'Formula') %>% mutate(Status = 'HMO + probiotic')

# The four cells of the factorial, all on the formula background
syn2x2 = rbind(
  formula_scfa     %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Formula'),
  hmo_scfa,
  formula_pro_scfa %>% mutate(Status = 'Probiotic'),
  hmo_pro_scfa) %>%
  mutate(Status = factor(Status, levels = c('Formula','HMO','Probiotic','HMO + probiotic')),
         group = factor(group, levels = preterm_levels),
         sub = factor(sub, levels = scfa_levels),
         Class = factor(scfa_class[as.character(sub)], levels = scfa_class_levels))

syn2x2_effects = syn2x2 %>%
  select(ID, group, sub, Class, Status, value) %>%
  pivot_wider(names_from = Status, values_from = value) %>%
  mutate(HMO_effect = HMO - Formula,
         Pro_effect = Probiotic - Formula,
         Both_effect = `HMO + probiotic` - Formula,
         Expected = HMO_effect + Pro_effect,
         Synergy = Both_effect - Expected)

# Each arm against the untreated formula baseline
syn2x2_stat = syn2x2 %>%
  arrange(sub, Status, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(value ~ Status, paired = TRUE, p.adjust.method = 'none',
                       ref.group = 'Formula') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

# Whether the combined effect exceeds the sum of the two single-agent effects
syn2x2_synergy_stat = syn2x2_effects %>%
  arrange(sub, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(Synergy ~ 1, mu = 0) %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

syn2x2_summary = syn2x2_effects %>%
  group_by(sub, Class) %>%
  summarise(median_formula = median(Formula),
            median_hmo = median(HMO_effect),
            median_pro = median(Pro_effect),
            median_both = median(Both_effect),
            median_expected = median(Expected),
            median_synergy = median(Synergy),
            n_synergistic = sum(Synergy > 0), .groups = 'drop') %>%
  left_join(syn2x2_synergy_stat %>% select(sub, p.adj, p.adj.signif), by = 'sub')
write.table(syn2x2_summary,file = 'Supplementary_table_S20.txt',
            sep = '\t',quote = F,row.names = TRUE)

# ---- HMO utilisation across the four cells -----------------------------------------------------------------------
hmo_use_2x2 = rbind(
  hmo_consumption(hmo_ext) %>% mutate(Status = 'HMO'),
  hmo_consumption(hmo_pro_ext) %>% mutate(Status = 'HMO + probiotic')) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = preterm_levels),
         sub = factor(sub, levels = hmo_levels),
         Status = factor(Status, levels = c('HMO','HMO + probiotic')),
         user = consumed_pct > utiliser_cutoff)

hmo_use_2x2_summary = hmo_use_2x2 %>%
  group_by(Status, sub) %>%
  summarise(n_users = sum(user, na.rm = TRUE),
            pct_users = 100 * n_users / n(),
            median_consumed = median(consumed_pct, na.rm = TRUE), .groups = 'drop')

hmo_use_2x2_stat = hmo_use_2x2 %>%
  arrange(sub, Status, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(consumed_pct ~ Status, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

# ---- Engraftment with and without the prebiotic ------------------------------------------------------------------
hmo_pro_engraft = rbind(
  engraftment(formula_pro_ext, 'Formula', pro_names) %>% mutate(Status = 'Probiotic'),
  engraftment(hmo_pro_ext, 'Formula + HMO', pro_names) %>% mutate(Status = 'HMO + probiotic')) %>%
  left_join(meta %>% select(ID, group), by = 'ID') %>%
  mutate(group = factor(group, levels = preterm_levels),
         species = factor(species, levels = pro_names),
         Status = factor(Status, levels = c('Probiotic','HMO + probiotic')))

hmo_pro_engraft_stat = hmo_pro_engraft %>%
  arrange(species, Status, ID) %>%
  group_by(species) %>%
  rstatix::wilcox_test(final ~ Status, paired = TRUE, p.adjust.method = 'none') %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

hmo_pro_engraft_summary = hmo_pro_engraft %>%
  group_by(Status, species) %>%
  summarise(median_final = median(final), .groups = 'drop') %>%
  pivot_wider(names_from = Status, values_from = median_final) %>%
  mutate(ratio = `HMO + probiotic` / Probiotic) %>%
  left_join(hmo_pro_engraft_stat %>% select(species, p.adj, p.adj.signif), by = 'species') %>%
  arrange(desc(ratio))
write.table(hmo_pro_engraft_summary, "Supplementary_Table_S21.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Plots -------------------------------------------------------------------------------------------------------
syn2x2_colors = c("Formula" = "grey70", "HMO" = "#E7B800",
                  "Probiotic" = "#00AFBB", "HMO + probiotic" = "#b2182b")

syn2x2_plot = ggplot(syn2x2, aes(x = Status, y = value, fill = Status))+
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6)+
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.35, colour = 'grey25')+
  ggh4x::facet_nested_wrap(~ Class + sub, scales = 'free_y', nrow = 1)+
  scale_fill_manual(values = syn2x2_colors, name = "Condition")+
  theme_bw()+
  labs(x = "", y = "Fermentation product, log(mM+1)",
       title = "Prebiotic and synbiotic intervention on the formula background")+
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8), plot.title = element_text(size = 11))

syn2x2_effect_long = syn2x2_effects %>%
  select(ID, group, sub, Class, HMO_effect, Pro_effect, Both_effect) %>%
  pivot_longer(cols = c(HMO_effect, Pro_effect, Both_effect),
               names_to = 'Arm', values_to = 'Diff') %>%
  mutate(Arm = factor(Arm, levels = c('HMO_effect','Pro_effect','Both_effect'),
                      labels = c('HMO','Probiotic','HMO + probiotic')))

syn2x2_effect_plot = ggplot(syn2x2_effect_long, aes(x = Arm, y = Diff, fill = Arm))+
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40')+
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.size = 0.8)+
  ggh4x::facet_nested_wrap(~ Class + sub, scales = 'free_y', nrow = 1)+
  scale_fill_manual(values = syn2x2_colors[-1], name = "Arm")+
  theme_bw()+
  labs(x = "", y = "Change from formula baseline, log(mM+1)")+
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8))


# ---- (S)-lactate, the only product with a significant synergy term -----------------------------------------------
# The four medians are computed per individual and are not additive, so the comparison is drawn per infant instead:
# each line runs from the effect predicted by adding the two single-agent arms to the effect actually observed when
# both are given, and its slope is that infant's synergy term.
lactate_synergy_data = syn2x2_effects %>%
  filter(sub == '(S)-Lactate') %>%
  select(ID, group, HMO_effect, Pro_effect, Expected, Both_effect, Synergy) %>%
  mutate(direction = factor(ifelse(Synergy > 0, "Positive", "Negative"),
                            levels = c("Positive","Negative")))

lactate_slope_data = lactate_synergy_data %>%
  pivot_longer(c(Expected, Both_effect), names_to = 'Arm', values_to = 'effect') %>%
  mutate(Arm = factor(Arm, levels = c('Expected','Both_effect'),
                      labels = c('Additive expectation','HMO + probiotic')))

lactate_slope_median = lactate_slope_data %>%
  group_by(Arm) %>%
  summarise(effect = median(effect), .groups = 'drop')

lactate_synergy_label = lactate_synergy_data %>%
  summarise(n_pos = sum(Synergy > 0), n = n(),
            med = median(Synergy)) %>%
  mutate(label = sprintf("Synergy: median %.2e, %d/%d positive", med, n_pos, n))

lactate_synergy_plot = ggplot(lactate_slope_data, aes(x = Arm, y = effect))+
  geom_line(aes(group = ID, colour = direction), alpha = 0.5, linewidth = 0.7)+
  geom_point(aes(colour = direction), size = 2, alpha = 0.7)+
  geom_line(data = lactate_slope_median, aes(group = 1),
            colour = 'grey15', linewidth = 1.7)+
  geom_point(data = lactate_slope_median, colour = 'grey15', size = 3.5)+
  geom_text(data = lactate_slope_median,
            aes(label = format(effect, digits = 3, scientific = TRUE),
                hjust = ifelse(Arm == 'Additive expectation', 1.18, -0.18)),
            colour = 'grey15', size = 4.5)+
  annotate("text", x = 1.5, y = Inf, vjust = 1.6, size = 4.7, colour = 'grey25',
           label = lactate_synergy_label$label)+
  scale_colour_manual(values = c("Positive" = "#b2182b", "Negative" = "#2166ac"),
                      name = "Synergy term")+
  scale_x_discrete(expand = expansion(mult = c(0.40, 0.30)))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14)))+
  theme_bw()+
  labs(x = "", y = "Change in (S)-lactate from baseline, log(mM+1)")+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = 'bottom')
ggsave(plot = lactate_synergy_plot, filename = 'Figures/lactate_synergy_plot.jpeg',
       dpi = 1200, width = 15, height = 14.4, units = 'cm')


# ---- The synergy term for all seven products ---------------------------------------------------------------------
syn_dat = syn2x2_synergy %>% mutate(sub = factor(sub, levels = scfa_levels))

panel_stats = syn_dat %>%
  group_by(Class, sub) %>%
  summarise(lo = min(Synergy, na.rm = TRUE), hi = max(Synergy, na.rm = TRUE),
            n_pos = sum(Synergy > 0), n = n(), .groups = 'drop') %>%
  mutate(rng = pmax(hi - lo, 1e-6),
         y_star = hi + 0.08 * rng,
         lab = sprintf('%d/%d', n_pos, n))

y_lab_pos = min(panel_stats$lo) - 0.10 * diff(range(syn_dat$Synergy, na.rm = TRUE))

syn2x2_sig = syn2x2_synergy_stat %>%
  filter(p.adj < 0.05) %>%
  select(sub, p.adj.signif) %>%
  left_join(panel_stats %>% select(Class, sub, y_star), by = 'sub')

syn2x2_synergy_plot = ggplot(syn_dat, aes(x = sub, y = Synergy))+
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40')+
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.size = 0.8, fill = 'grey75')+
  geom_text(data = syn2x2_sig, aes(x = sub, y = y_star, label = p.adj.signif),
            inherit.aes = FALSE, size = 5)+
  geom_text(data = panel_stats, aes(x = sub, y = y_lab_pos, label = lab),
            inherit.aes = FALSE, size = 3, colour = 'grey30')+
  ggh4x::facet_nested(~ Class + sub, scales = 'free_x', space = 'free_x')+
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12)))+
  theme_bw()+
  labs(x = "", y = "Synergy term (\u0394log(mM+1))")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 8),
        panel.spacing.x = unit(1.5, 'pt'))

ggsave(plot = syn2x2_synergy_plot, filename = 'Figures/syn2x2_synergy_plot.jpeg',
       dpi = 1200, width = 20, height = 12, units = 'cm')



hmo_use_2x2_plot = ggplot(hmo_use_2x2_summary, aes(x = sub, y = pct_users, fill = Status))+
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85)+
  scale_fill_manual(values = c("HMO" = "#E7B800", "HMO + probiotic" = "#b2182b"), name = "")+
  theme_bw()+
  labs(title = "HMO utilisation with and without the probiotic strains",
       x = "", y = paste0("Individuals with >", utiliser_cutoff, "% consumption (%)"))+
  theme(axis.text = element_text(size = 13), plot.title = element_blank(),
        axis.title.y = element_text(size = 15),legend.text = element_text(size = 15),
        legend.position = 'bottom')
ggsave(plot = hmo_use_2x2_plot,filename = 'Figures/hmo_use_2x2_plot.jpeg',
       dpi = 1200,height = 12,width = 13,units = 'cm')

hmo_consume_2x2_plot = ggplot(hmo_use_2x2, aes(x = sub, y = consumed_pct, fill = Status))+
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.size = 0.8,
               position = position_dodge(width = 0.75))+
  scale_fill_manual(values = c("HMO" = "#E7B800", "HMO + probiotic" = "#b2182b"), name = "")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))+
  theme_bw()+
  labs(x = "", y = "HMO consumed (%)")+
  theme(axis.text = element_text(size = 13), legend.position = 'bottom',
        axis.title.y = element_text(size = 15),legend.text = element_text(size = 15))
ggsave(plot = hmo_consume_2x2_plot, filename = 'Figures/hmo_consume_2x2_plot.jpeg',
       dpi = 1200, height = 12, width = 13, units = 'cm')

# ================================================================================================================== #
# Part 11  Robustness of the lactate results
#
#   Manuscript: Results 6
#   Both analyses follow from Part 7.2, which found the enantiomer split undetermined
#   for several resident strains. library() calls removed; Part 0 already loads them.
# ================================================================================================================== #

# ##################################################################################################################
# 11.1  Every lactate-dependent test repeated with the two enantiomers pooled
# ##################################################################################################################

# ================================================================================================== #
# Sensitivity analysis 4: are the lactate results robust to the enantiomer assignment?
#
# Why this is needed
#   sensitivity_01b showed that for a substantial number of resident strains the mtf solution does not
#   determine which lactate enantiomer carries the flux. The total is conserved across permutations of
#   the reaction order, but the split between EX_lac_D(e) and EX_lac_L(e) flips. The organisms
#   affected include Staphylococcus sp. HMSC070D05, Escherichia coli SE11, Lactobacillus rhamnosus,
#   Lactobacillus acidophilus, Ruminococcus gnavus and Veillonella, which between them are the
#   non-Bifidobacterium contributors that the manuscript attributes to one enantiomer or the other.
#
#   Because glpk is deterministic for a given problem, every infant's simulation resolves the tie the
#   same way. The resulting pattern is internally consistent across all 51 subjects and therefore
#   looks like a biological signal, when it reflects the solver rather than the network. Bifidobacterium
#   is not affected: its models carry only one of the two dehydrogenase routes, so its (S)-lactate
#   contribution is determined.
#
#   This script repeats every lactate-dependent analysis with the two enantiomers pooled. If the
#   conclusions hold on total lactate, they can be reported on total lactate and the enantiomer-level
#   attribution dropped. If they do not, the affected claims rest on the solver's choice and have to
#   be rewritten before submission.
#
# Prerequisites
#   Run after Parts 3 to 8 of major_script.R, with the extracted objects in the workspace:
#     breast_ext, formula_ext, hmo_ext, hmo_pro_ext, and the probiotic-arm extraction used to build
#     formula_pro_scfa. CHECK the name of that last object at the marked line before running.
#
# Output
#   Console comparison of every affected test, original against pooled.
#   sensitivity/lactate_pooling_comparison.txt
# ================================================================================================== #


dir.create('sensitivity', showWarnings = FALSE)

# -------------------------------------------------------------------------------------------------- #
# Pooled extraction. Identical to extract_scfa except that the two lactate exchanges are summed while
# the values are still concentrations. The log transform is applied afterwards, because log(a + 1)
# plus log(b + 1) is not log(a + b + 1) and pooling after the transform would be wrong.
# -------------------------------------------------------------------------------------------------- #
scfa_levels_pooled = c('Acetate', 'Propionate', 'Butyrate', 'Formate', 'Succinate', 'Lactate')
scfa_class_pooled  = c(Acetate = 'SCFA', Propionate = 'SCFA', Butyrate = 'SCFA', Formate = 'SCFA',
                       Succinate = 'Organic acid', Lactate = 'Organic acid')

extract_scfa_pooled = function(extracted, metadata, diet_type) {
  extracted$scfa_end %>%
    mutate(sub = ifelse(sub %in% names(name_map), name_map[sub], sub),
           sub = ifelse(sub %in% c('(R)-Lactate', '(S)-Lactate'), 'Lactate', sub)) %>%
    group_by(ID, sub) %>%
    summarise(value = sum(value), .groups = 'drop') %>%
    complete(ID, sub, fill = list(value = 0)) %>%
    left_join(metadata %>% select(ID, group), by = 'ID') %>%
    mutate(group = factor(group, levels = group_levels),
           sub = factor(sub, levels = scfa_levels_pooled),
           Class = factor(scfa_class_pooled[as.character(sub)],
                          levels = levels(factor(scfa_class_pooled))),
           value = log(value + 1),
           Diet = diet_type)
}

# -------------------------------------------------------------------------------------------------- #
# A sanity check before anything else. Pooling must leave the other five products untouched and must
# give a lactate value at least as large as either enantiomer alone.
# -------------------------------------------------------------------------------------------------- #
breast_scfa_pooled  = extract_scfa_pooled(breast_ext,  meta, 'Breast')
formula_scfa_pooled = extract_scfa_pooled(formula_ext, meta, 'Formula')

pool_check = breast_scfa_pooled %>%
  filter(sub != 'Lactate') %>%
  left_join(breast_scfa %>% select(ID, sub, value_orig = value), by = c('ID', 'sub')) %>%
  summarise(max_abs_diff = max(abs(value - value_orig), na.rm = TRUE))

cat('\n--- Sanity check: non-lactate products must be unchanged by pooling ---\n')
print(pool_check)
if (pool_check$max_abs_diff > 1e-10) {
  warning('Pooling altered a non-lactate product. Check name_map and the sub labels before ',
          'interpreting anything below.')
}

# ================================================================================================== #
# Section 2: group comparisons
# ================================================================================================== #

group_test = function(df, label) {
  df %>%
    group_by(sub) %>%
    rstatix::wilcox_test(value ~ group, p.adjust.method = 'none') %>%
    ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    mutate(analysis = label) %>%
    select(analysis, sub, group1, group2, p, p.adj)
}

sec2_original = bind_rows(
  group_test(breast_scfa  %>% filter(sub %in% c('(R)-Lactate', '(S)-Lactate')), 'Breast, enantiomers'),
  group_test(formula_scfa %>% filter(sub %in% c('(R)-Lactate', '(S)-Lactate')), 'Formula, enantiomers'))

sec2_pooled = bind_rows(
  group_test(breast_scfa_pooled  %>% filter(sub == 'Lactate'), 'Breast, pooled'),
  group_test(formula_scfa_pooled %>% filter(sub == 'Lactate'), 'Formula, pooled'))

cat('\n--- Section 2, lactate group comparisons ---\n')
print(bind_rows(sec2_original, sec2_pooled), n = 40)

# Medians alongside, since a p value alone does not say whether the direction changed
sec2_medians = bind_rows(
  breast_scfa  %>% filter(sub %in% c('(R)-Lactate','(S)-Lactate')) %>% mutate(Diet = 'Breast'),
  formula_scfa %>% filter(sub %in% c('(R)-Lactate','(S)-Lactate')) %>% mutate(Diet = 'Formula'),
  breast_scfa_pooled  %>% filter(sub == 'Lactate') %>% mutate(Diet = 'Breast'),
  formula_scfa_pooled %>% filter(sub == 'Lactate') %>% mutate(Diet = 'Formula')) %>%
  group_by(Diet, sub, group) %>%
  summarise(median = median(value), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = median)

cat('\n--- Section 2, medians by group ---\n')
print(sec2_medians, n = 20)

# ================================================================================================== #
# Section 6: the 2 x 2 factorial, which is where this matters most
# ================================================================================================== #

# CHECK THIS LINE. formula_pro_scfa in the main script is built from the probiotic arm on the formula
# background; substitute the name of the extracted object it came from.
formula_pro_ext_name = 'formula_pro_ext'
formula_pro_ext_obj  = get(formula_pro_ext_name)

hmo_scfa_pooled         = extract_scfa_pooled(hmo_ext, meta, 'Formula')         %>% mutate(Status = 'HMO')
hmo_pro_scfa_pooled     = extract_scfa_pooled(hmo_pro_ext, meta, 'Formula')     %>% mutate(Status = 'HMO + probiotic')
formula_pro_scfa_pooled = extract_scfa_pooled(formula_pro_ext_obj, meta, 'Formula') %>% mutate(Status = 'Probiotic')

syn2x2_pooled = rbind(
  formula_scfa_pooled %>% filter(ID %in% preterm_ids) %>% mutate(Status = 'Formula'),
  hmo_scfa_pooled,
  formula_pro_scfa_pooled,
  hmo_pro_scfa_pooled) %>%
  mutate(Status = factor(Status, levels = c('Formula', 'HMO', 'Probiotic', 'HMO + probiotic')),
         group = factor(group, levels = preterm_levels))

syn2x2_effects_pooled = syn2x2_pooled %>%
  select(ID, group, sub, Class, Status, value) %>%
  pivot_wider(names_from = Status, values_from = value) %>%
  mutate(HMO_effect = HMO - Formula,
         Pro_effect = Probiotic - Formula,
         Both_effect = `HMO + probiotic` - Formula,
         Expected = HMO_effect + Pro_effect,
         Synergy = Both_effect - Expected)

syn2x2_synergy_stat_pooled = syn2x2_effects_pooled %>%
  arrange(sub, ID) %>%
  group_by(sub) %>%
  rstatix::wilcox_test(Synergy ~ 1, mu = 0) %>%
  ungroup() %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance('p.adj')

synergy_comparison = bind_rows(
  syn2x2_synergy_stat %>% mutate(analysis = 'Original, seven products') %>%
    select(analysis, sub, n, statistic, p, p.adj, p.adj.signif),
  syn2x2_synergy_stat_pooled %>% mutate(analysis = 'Pooled lactate, six products') %>%
    select(analysis, sub, n, statistic, p, p.adj, p.adj.signif))

cat('\n--- Section 6, synergy term, original against pooled ---\n')
print(synergy_comparison, n = 40)

# The three numbers the text reports, recomputed on pooled lactate
synergy_key_numbers = syn2x2_effects_pooled %>%
  filter(sub == 'Lactate') %>%
  summarise(n = n(),
            n_positive = sum(Synergy > 0),
            median_synergy = median(Synergy),
            median_expected = median(Expected),
            median_both = median(Both_effect),
            ratio = median(Both_effect) / median(Expected),
            min_both = min(Both_effect), max_both = max(Both_effect),
            min_exp = min(Expected), max_exp = max(Expected))

cat('\n--- Section 6, key numbers on pooled lactate ---\n')
print(as.data.frame(synergy_key_numbers))

# The same numbers for (S)-lactate as currently reported, for side by side reading
synergy_key_numbers_S = syn2x2_effects %>%
  filter(sub == '(S)-Lactate') %>%
  summarise(n = n(),
            n_positive = sum(Synergy > 0),
            median_synergy = median(Synergy),
            median_expected = median(Expected),
            median_both = median(Both_effect),
            ratio = median(Both_effect) / median(Expected),
            min_both = min(Both_effect), max_both = max(Both_effect),
            min_exp = min(Expected), max_exp = max(Expected))

cat('\n--- Section 6, the same numbers as currently reported for (S)-lactate ---\n')
print(as.data.frame(synergy_key_numbers_S))

# ================================================================================================== #
# Section 4: probiotic intervention, and section 3: the paired diet comparison
#
# Neither reported a significant lactate result, so these are checked only to confirm that pooling
# does not create one. A newly significant result would matter as much as a lost one.
# ================================================================================================== #

paired_test = function(df, by, label) {
  df %>%
    group_by(sub) %>%
    rstatix::wilcox_test(as.formula(paste('value ~', by)), paired = TRUE, p.adjust.method = 'none') %>%
    ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    mutate(analysis = label) %>%
    select(analysis, sub, group1, group2, p, p.adj)
}

diet_paired_pooled = bind_rows(breast_scfa_pooled, formula_scfa_pooled) %>%
  paired_test('Diet', 'Section 3, breast against formula, pooled')

cat('\n--- Section 3, paired diet comparison on pooled lactate ---\n')
print(diet_paired_pooled %>% filter(sub == 'Lactate'))

# -------------------------------------------------------------------------------------------------- #
# Everything written out together
# -------------------------------------------------------------------------------------------------- #
comparison_out = bind_rows(
  sec2_original, sec2_pooled,
  synergy_comparison %>% transmute(analysis, sub = as.character(sub), group1 = NA_character_,
                                   group2 = NA_character_, p, p.adj),
  diet_paired_pooled)

write.table(comparison_out, file = 'sensitivity/lactate_pooling_comparison.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

# ##################################################################################################################
# 11.2  Organisms driving the (S)-lactate synergy, checked against the degeneracy classification
# ##################################################################################################################

# ================================================================================================== #
# Sensitivity analysis 5: which organisms drive the (S)-lactate synergy?
#
# Why this decides the question
#   sensitivity_04 showed that the synergy significant on (S)-lactate (adjusted p = 0.0074) is not
#   significant on pooled lactate (adjusted p = 0.062). sensitivity_01b showed that the enantiomer
#   split is not determined by mtf for a number of strains. The (S)-lactate result is defensible if
#   and only if the increase under the combined intervention comes from organisms whose enantiomer
#   assignment is determined, so this script decomposes the synergy term by organism and checks each
#   contributor against the degeneracy classification.
#
# Naming, which the first version of this script got wrong
#   The uniqueness table labels the seven introduced strains by their short display names (B.infantis)
#   while the flux tables label them by their model identifiers (PRO_BINF). Matching on a cleaned
#   organism string therefore silently failed for all seven, and since the introduced strains are the
#   largest contributors the verdict came out inverted. The mapping is applied explicitly below and a
#   check reports any contributor that could not be classified, so a future mismatch is visible rather
#   than silent.
#
# Prerequisites
#   formula_ext, hmo_ext, hmo_pro_ext and the probiotic-arm extraction, plus `uniqueness` in the
#   workspace. That object must be the version produced AFTER the absolute-tolerance correction; the
#   original classification in sensitivity_01b treated fluxes of order 1e-13 as sign changes and
#   would overstate the degenerate share here.
#
# Output
#   Console ranking of the organisms contributing to the (S)-lactate synergy, annotated with whether
#   each one's enantiomer assignment is determined, and the share of the total each class carries.
#   sensitivity/lactate_synergy_sources.txt
# ================================================================================================== #


dir.create('sensitivity', showWarnings = FALSE)

formula_pro_ext_name = 'formula_pro_ext'          # CHECK
formula_pro_ext_obj  = get(formula_pro_ext_name)

lac_L = 'EX_lac_L(e)'     # (S)-lactate
lac_D = 'EX_lac_D(e)'     # (R)-lactate

stopifnot(exists('uniqueness'), exists('pro_labels'))

# -------------------------------------------------------------------------------------------------- #
# Net lactate flux of every organism in each of the four arms. complete() fills the organisms absent
# from an arm with zero, which is what the synergy decomposition needs: a probiotic contributes
# nothing to the two arms it was not added to, and that zero is a real zero rather than a missing
# value.
# -------------------------------------------------------------------------------------------------- #
arm_flux = bind_rows(
  formula_ext$flux         %>% mutate(Arm = 'Formula'),
  hmo_ext$flux             %>% mutate(Arm = 'HMO'),
  formula_pro_ext_obj$flux %>% mutate(Arm = 'Probiotic'),
  hmo_pro_ext$flux         %>% mutate(Arm = 'Both')) %>%
  filter(rea %in% c(lac_L, lac_D), ID %in% preterm_ids) %>%
  group_by(Arm, ID, spec, rea) %>%
  summarise(net = sum(mflux, na.rm = TRUE), .groups = 'drop') %>%
  complete(Arm, ID, spec, rea, fill = list(net = 0))

# -------------------------------------------------------------------------------------------------- #
# The synergy term decomposed by organism. It is additive over organisms because the community flux
# of a product is the sum of its members' exchange fluxes, so each organism's contribution to the
# synergy is its own combined-arm flux minus what the two single-agent arms would have predicted.
# -------------------------------------------------------------------------------------------------- #
synergy_by_org = arm_flux %>%
  pivot_wider(names_from = Arm, values_from = net, values_fill = 0) %>%
  mutate(HMO_effect = HMO - Formula,
         Pro_effect = Probiotic - Formula,
         Both_effect = Both - Formula,
         Synergy = Both_effect - (HMO_effect + Pro_effect))

org_ranking = synergy_by_org %>%
  group_by(rea, spec) %>%
  summarise(n_subjects = n_distinct(ID[abs(Synergy) > 0]),
            total_synergy = sum(Synergy),
            median_synergy = median(Synergy),
            .groups = 'drop') %>%
  group_by(rea) %>%
  mutate(share = 100 * total_synergy / sum(total_synergy[total_synergy > 0])) %>%
  ungroup() %>%
  mutate(compound = ifelse(rea == lac_L, '(S)-Lactate', '(R)-Lactate'))

# ================================================================================================== #
# Classification of each contributor
# ================================================================================================== #

# -------------------------------------------------------------------------------------------------- #
# Two naming systems have to be reconciled. Resident models appear in the uniqueness table as file
# names (ASV006_Staphylococcus_sp_HMSC070D05.mat) and in the flux tables as organism names
# (Staphylococcus_sp_HMSC070D05); stripping the ASV prefix and the extension aligns them. The seven
# introduced strains appear as display labels (B.infantis) and as model identifiers (PRO_BINF), which
# no string cleaning can reconcile, so pro_labels is used directly.
# -------------------------------------------------------------------------------------------------- #
clean_name = function(x) {
  x %>% sub('^ASV[0-9]+_', '', .) %>% sub('\\.mat$', '', .) %>% tolower()
}

# Key used for matching: the display label for an introduced strain, the cleaned organism name
# otherwise. Applied identically to both sides.
to_key = function(x) {
  ifelse(x %in% names(pro_labels), unname(pro_labels[x]), clean_name(x))
}

degenerate_keys = uniqueness %>%
  filter(compound %in% c('(R)-Lactate', '(S)-Lactate'),
         class %in% c('sign varies', 'sign determined, magnitude varies')) %>%
  distinct(set, strain) %>%
  mutate(key = ifelse(set == 'Probiotic', strain, clean_name(strain))) %>%
  pull(key) %>% unique()

# Every key that the uniqueness table knows about at all, so that a contributor absent from it is
# reported as unclassified rather than defaulting to determined.
known_keys = uniqueness %>%
  filter(compound %in% c('(R)-Lactate', '(S)-Lactate')) %>%
  distinct(set, strain) %>%
  mutate(key = ifelse(set == 'Probiotic', strain, clean_name(strain))) %>%
  pull(key) %>% unique()

swap_keys = uniqueness %>%
  filter(compound %in% c('(R)-Lactate', '(S)-Lactate'), class == 'sign varies') %>%
  distinct(set, strain) %>%
  mutate(key = ifelse(set == 'Probiotic', strain, clean_name(strain))) %>%
  pull(key) %>% unique()

org_ranking = org_ranking %>%
  mutate(match_key = to_key(spec),
         enantiomer_assignment = case_when(
           match_key %in% swap_keys  ~ 'label unreliable',
           match_key %in% known_keys ~ 'label reliable',
           TRUE                      ~ 'unclassified'))

# -------------------------------------------------------------------------------------------------- #
# Matching audit. Any contributor landing in 'unclassified' carries weight that the verdict cannot
# account for, so its share is reported rather than being absorbed into either class.
# -------------------------------------------------------------------------------------------------- #
match_audit = org_ranking %>%
  filter(compound == '(S)-Lactate', total_synergy > 0) %>%
  group_by(enantiomer_assignment) %>%
  summarise(n_organisms = n(), total = sum(total_synergy), .groups = 'drop') %>%
  mutate(pct = round(100 * total / sum(total), 1))

cat('\n--- Matching audit: every positive contributor must be classified ---\n')
print(match_audit)

unclassified = org_ranking %>%
  filter(compound == '(S)-Lactate', total_synergy > 0, enantiomer_assignment == 'unclassified') %>%
  arrange(desc(total_synergy)) %>%
  select(spec, match_key, n_subjects, total_synergy, share)

if (nrow(unclassified) > 0) {
  cat('\n--- Contributors not found in the uniqueness table ---\n')
  print(unclassified, n = 30)
}

# ================================================================================================== #
# The decision
# ================================================================================================== #
top_S = org_ranking %>%
  filter(compound == '(S)-Lactate', total_synergy > 0) %>%
  arrange(desc(total_synergy)) %>%
  select(spec, n_subjects, total_synergy, share, enantiomer_assignment)

cat('\n--- Organisms contributing positively to the (S)-lactate synergy ---\n')
print(top_S, n = 25)

verdict = match_audit

cat('\n--- Share of the positive (S)-lactate synergy by assignment status ---\n')
print(verdict)

# -------------------------------------------------------------------------------------------------- #
# The classification of the largest contributors in full, since a single organism carrying most of
# the signal decides the question on its own and its class should be read directly rather than
# through an aggregate.
# -------------------------------------------------------------------------------------------------- #
top_detail = org_ranking %>%
  filter(compound == '(S)-Lactate', total_synergy > 0) %>%
  arrange(desc(total_synergy)) %>%
  slice_head(n = 6) %>%
  left_join(uniqueness %>%
              filter(compound %in% c('(R)-Lactate', '(S)-Lactate')) %>%
              mutate(key = ifelse(set == 'Probiotic', strain, clean_name(strain))) %>%
              select(key, medium, uniq_compound = compound, f_ref, f_min, f_max, class),
            by = c('match_key' = 'key')) %>%
  select(spec, share, uniq_compound, medium, f_ref, f_min, f_max, class)

cat('\n--- Enantiomer classification of the six largest contributors ---\n')
print(top_detail, n = 40)

write.table(org_ranking, file = 'sensitivity/lactate_synergy_sources.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

# ================================================================================================================== #
# Part 12  Supplement composition and cross-arm comparison
#
#   Manuscript: Results 5, supporting; Table S17
# ================================================================================================================== #


# ---- Composition of the predicted supplements --------------------------------------------------------------------
# The limiting substrates differ between the two media because the formula medium already supplies the trace
# nutrients that limit growth under breast milk, so amino acids become the bottleneck instead.
supp_recipe = rbind(breast_supplements  %>% mutate(Diet = 'Breast'),
                    formula_supplements %>% mutate(Diet = 'Formula')) %>%
  filter(round == 1) %>%
  distinct(Diet, ID, react) %>%
  count(Diet, react, name = 'n_infants')

supp_recipe_size = rbind(breast_supplements  %>% mutate(Diet = 'Breast'),
                         formula_supplements %>% mutate(Diet = 'Formula')) %>%
  filter(round == 1) %>%
  distinct(Diet, ID, react) %>%
  count(Diet, ID, name = 'n_substrates') %>%
  group_by(Diet) %>%
  summarise(median = median(n_substrates),
            min = min(n_substrates),
            max = max(n_substrates), .groups = 'drop')

supp_table_S17 = rbind(breast_supplements  %>% mutate(Diet = 'Breast'),
                       formula_supplements %>% mutate(Diet = 'Formula')) %>%
  filter(round == 1) %>%
  distinct(Diet, ID, group, react) %>%
  mutate(label = ifelse(react %in% names(mediac_names), mediac_names[react], react)) %>%
  group_by(Diet, ID, group) %>%
  summarise(n_substrates = n(),
            substrates = paste(sort(label), collapse = "; "), .groups = 'drop') %>%
  select(Medium = Diet, Infant = ID, Group = group,
         `Number of substrates` = n_substrates, Substrates = substrates) %>%
  arrange(Medium, Group, Infant)

write.table(supp_table_S17, "Supplementary_Table_S17.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

supp_recipe_top = supp_recipe %>%
  group_by(Diet) %>%
  slice_max(n_infants, n = 40) %>%
  ungroup() %>%
  mutate(label = ifelse(react %in% names(mediac_names), mediac_names[react], react))

supp_recipe_plot = ggplot(supp_recipe_top,
                          aes(x = n_infants, y = reorder(label, n_infants), fill = Diet))+
  geom_col(width = 0.7)+
  facet_wrap(~ Diet, ncol = 1, scales = 'free_y')+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  theme_bw()+
  labs(x = "Infants with this substrate in the predicted supplement", y = "")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.position = 'none')
ggsave(plot = supp_recipe_plot,filename = 'Figures/supp_recipe_plot.jpeg',
       dpi = 1600,height = 40,width = 16,units = 'cm')

# ---- Where the supplemented nutrients end up ---------------------------------------------------------------------
# Supplementation does not change the total output but redistributes flux between taxa, and most of the increase
# accrues to the family that already dominates the preterm gut.
supp_flux_class = supp_producer_shift %>%
  mutate(taxon_group = case_when(
    grepl("Bacteroides|Parabacteroides|Prevotella", spec) ~ "Bacteroidetes",
    grepl("Bifidobacterium", spec)                        ~ "Bifidobacterium",
    grepl("Escherichia|Enterobacter|Citrobacter|Yersinia|Klebsiella|Serratia|Pantoea", spec)
    ~ "Enterobacteriaceae",
    TRUE ~ "Other")) %>%
  filter(rea %in% scfa_levels) %>%
  pivot_longer(c(`Non-treated`, Supplemented),
               names_to = 'Arm', values_to = 'flux') %>%
  filter(flux > 0) %>%
  group_by(Diet, Arm, taxon_group) %>%
  summarise(total = sum(flux), .groups = 'drop')

supp_shift_class = supp_flux_class %>%
  mutate(Arm = factor(Arm,
                      levels = c('Non-treated','Supplemented'),
                      labels = c('Non-treated','Supp-treated')),
         taxon_group = factor(taxon_group,
                              levels = c("Enterobacteriaceae","Bacteroidetes",
                                         "Bifidobacterium","Other"))) %>%
  group_by(Diet, Arm) %>%
  mutate(pct = 100 * total / sum(total)) %>%
  ungroup()

supp_shift_plot = ggplot(supp_shift_class, aes(x = Arm, y = pct, fill = taxon_group))+
  geom_col(width = 0.6)+
  facet_wrap(~ Diet, nrow = 1)+
  scale_fill_manual(values = class_colors, name = "")+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  labs(x = "", y = "Share of net secretion flux (%)")+
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 11, face = 'italic'),
        legend.position = 'bottom')
ggsave(plot = supp_shift_plot,filename = 'Figures/supp_shift_plot.jpeg',
       dpi = 1200,width = 15,height = 13,units = 'cm')



# ---- The three intervention arms side by side --------------------------------------------------------------------
arm_compare = rbind(
  supp_effect_summary %>% mutate(Arm = 'Supplement'),
  pro_effect_summary  %>% mutate(Arm = 'Probiotic'),
  syn_effect_summary  %>% mutate(Arm = 'Combined')) %>%
  mutate(Arm = factor(Arm, levels = c('Supplement','Probiotic','Combined')),
         Class = factor(scfa_class[as.character(sub)], levels = scfa_class_levels))

arm_compare_plot = ggplot(arm_compare, aes(x = Arm, y = median_diff, fill = Diet))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6)+
  geom_hline(yintercept = 0, colour = 'grey40', linewidth = 0.3)+
  ggh4x::facet_nested(~ Class + sub, scales = 'free_y')+
  scale_fill_manual(values = diet_colors, name = "Diet")+
  theme_bw()+
  labs(x = "", y = "Median change from baseline, log(mM+1)")+
  theme(axis.text.x = element_text(size = 13, angle = 40, hjust = 1),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.position = 'bottom')
ggsave(plot = arm_compare_plot,filename = 'Figures/arm_compare_plot.jpeg',
       dpi = 1200,height = 13,width = 28,units = 'cm')

# ---- Engraftment with and without the oligosaccharides -----------------------------------------------------------
# Only the two strains reconstructed with HMO degradation pathways gain from the addition of the substrate; the
# third bifidobacterium, which lacks those reactions, behaves like the strains of other genera.
engraft_hmo_plot = hmo_pro_engraft_summary %>%
  mutate(species = factor(species, levels = c('PRO_BINF','PRO_BBIF','PRO_BB12',
                                              'PRO_LGG','PRO_FER','PRO_REU','PRO_CBU'))) %>%
  ggplot(aes(x = species, y = ratio, fill = p.adj < 0.05))+
  geom_col(width = 0.6)+
  geom_hline(yintercept = 1, linetype = 'dashed', colour = 'grey40')+
  scale_x_discrete(labels = pro_labels)+
  scale_fill_manual(values = c('TRUE' = 'goldenrod2', 'FALSE' = 'grey70'),
                    labels = c('TRUE' = 'p.adj < 0.05', 'FALSE' = 'ns'), name = "")+
  theme_bw()+
  labs(x = "", y = "Fold change in final abundance with HMO")+
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1, face = 'italic'),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = 'none')
ggsave(plot = engraft_hmo_plot,filename = 'Figures/engraft_hmo_plot.jpeg',
       dpi = 1200,height = 12.5, width = 13.6,units = 'cm')

# ================================================================================================================== #
# End of analysis
# ================================================================================================================== #