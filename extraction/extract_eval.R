# ================================================================================================================== #
# Reduce one simulation result to the quantities the analysis needs.
#
# A single Eval object holds the full state of every arena grid cell at every time step and occupies well over a
# hundred gigabytes once loaded, so the objects cannot be held in memory together. This script reads one object,
# extracts four small tables from it, writes them out and exits. The analysis script then works from the extracted
# tables alone and never touches the original objects again.
#
# Only the exchange reactions actually present in an arena are requested: the formula medium carries no HMO, and in
# communities where no organism holds the corresponding exchange reaction the substance is absent from the media
# vector altogether, which would otherwise abort the extraction.
#
#   Usage: Rscript extract_eval.R <ID> <condition>
#   where condition is one of no_treat_breast, no_treat_formula, pro_breast, pro_formula,
#   supp_breast, supp_formula, syn_breast, syn_formula, hmo_formula, hmo_pro_formula
# ================================================================================================================== #
.libPaths("/home/sli/personal/sli/R_libs/4.3")
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript extract_eval.R <ID> <condition>")
}
Ind = args[1]
Cond = args[2]
suppressPackageStartupMessages({
  library(BacArena)
  library(dplyr)
  library(tidyr)
})

setwd("~/personal/sli/New_BacArena_project/")

scfas = c("EX_ac(e)","EX_ppa(e)","EX_succ(e)","EX_lac_D(e)","EX_lac_L(e)","EX_for(e)","EX_but(e)")
# Eight of the thirteen oligosaccharide exchanges carried by the breast medium are extracted: the six most
# frequently reported HMOs plus two further fucosylated species for which the model library provides comparable
# coverage (22 and 24 of the 191 modelled strains). The remaining five are represented by three or fewer models,
# or, in the case of free fucose, are degradation products rather than intact oligosaccharides.
hmo_reactions = c("2'-FL"  = "EX_2fuclac(e)",     "3'-FL"  = "EX_3fuclac(e)",
                  "DFL"    = "EX_dfuclac(e)",     "LDFT"   = "EX_lacdfucttr(e)",
                  "LNT"    = "EX_lacnttr(e)",     "LNnT"   = "EX_lacnnttr(e)",
                  "6'-SL"  = "EX_6slac(e)",       "3'-SL"  = "EX_3slac(e)")
flux_floor = 1e-9   # fluxes below this are numerical residue and are dropped to keep the tables small

# The file name suffix differs from the directory name for the untreated runs
suffix = sub("^no_treat_", "", Cond)
eval_file = file.path(paste0("Ind_eval_results_", Cond),
                      paste0(Ind, "_eval_result_", suffix, ".rds"))
if (!file.exists(eval_file)) stop("Missing: ", eval_file)

cat("Reading ", eval_file, "\n", sep = '')
e = readRDS(eval_file)

empty_sub = data.frame(sub = character(0), time = numeric(0), value = numeric(0),
                       ID = character(0), stringsAsFactors = FALSE)

# ---- 1. SCFA concentration at the final time step -----------------------------------------------------------------
scfa_present = intersect(scfas, e@mediac)
scfa_end = if (length(scfa_present) > 0) {
  plotSubCurve(e, mediac = scfa_present, ret_data = TRUE, unit = 'mM', useNames = TRUE) %>%
    subset(., time == max(time)) %>%
    select(sub, value) %>%
    mutate(ID = Ind)
} else {
  data.frame(sub = character(0), value = numeric(0), ID = character(0),
             stringsAsFactors = FALSE)
}

# ---- 2. HMO concentration over the whole simulation ---------------------------------------------------------------
hmo_present = intersect(unname(hmo_reactions), e@mediac)
hmo_curve = if (length(hmo_present) > 0) {
  plotSubCurve(e, mediac = hmo_present, ret_data = TRUE, unit = 'mM', useNames = TRUE) %>%
    select(sub, time, value) %>%
    mutate(ID = Ind)
} else {
  empty_sub
}

# ---- 3. Growth curve of every organism ----------------------------------------------------------------------------
growth = plotGrowthCurve(simlist = e, use_biomass = FALSE, ret_data = TRUE) %>%
  select(species, time, value) %>%
  mutate(ID = Ind)

# ---- 4. Exchange flux of every organism, for the fermentation products and the oligosaccharides --------------------
# Both sets are taken in one call so that the same table supports the producer attribution of the fermentation
# products and the question of which taxa actually consume the HMOs. Positive values are secretion and negative
# values uptake, so the oligosaccharide rows are read as uptake.
flux_reactions = intersect(c(scfas, unname(hmo_reactions)), e@mediac)
flux = if (length(flux_reactions) > 0) {
  plotReaActivity(simlist = e, reactions = flux_reactions, ret_data = TRUE) %>%
    filter(abs(mflux) > flux_floor) %>%
    select(spec, rea, time, mflux) %>%
    mutate(ID = Ind)
} else {
  data.frame(spec = character(0), rea = character(0), time = numeric(0),
             mflux = numeric(0), ID = character(0), stringsAsFactors = FALSE)
}

out = list(scfa_end = scfa_end, hmo_curve = hmo_curve, growth = growth, flux = flux)

# The media vector is identical across subjects of a condition, but a community missing some exchange reactions
# would give a short one, so the longest is kept
out_dir = paste0("extracted_", Cond)
dir.create(out_dir, showWarnings = FALSE)
mediac_file = file.path(out_dir, "_mediac.rds")

# Written once per condition. Concurrent array tasks can collide on this file, so a failed read is treated as
# "not yet written" rather than an error, and the write goes through a task-specific temporary file whose rename
# is atomic within the same filesystem.
need_write = TRUE
if (file.exists(mediac_file)) {
  existing = tryCatch(readRDS(mediac_file), error = function(err) NULL)
  if (!is.null(existing) && length(existing) >= length(e@mediac)) need_write = FALSE
}
if (need_write) {
  tmp = paste0(mediac_file, ".", Sys.getpid())
  saveRDS(e@mediac, tmp)
  file.rename(tmp, mediac_file)
}

rm(e); gc(verbose = FALSE)
saveRDS(out, file.path(out_dir, paste0(Ind, ".rds")))
cat("Finished ", Ind, " ", Cond, ": ",
    nrow(scfa_end), " SCFA rows, ", nrow(hmo_curve), " HMO rows, ",
    nrow(growth), " growth rows, ", nrow(flux), " flux rows\n", sep = '')