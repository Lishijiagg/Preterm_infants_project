.libPaths("/home/sli/personal/sli/R_libs/4.3")
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript simulate_arena.R <ID>")
}

Ind = args[1]
set.seed(sum(utf8ToInt(Ind))+1)

library(BacArena)
library(stringr)

meta = read.table("~/personal/sli/New_BacArena_project/metadata.tsv",
                  header = TRUE,row.names = 1,sep = '\t')

all_amount_data = read.table('~/personal/sli/New_BacArena_project/all_amount_data.txt',
                             sep = '\t',check.names = FALSE)

Diet = read.table("~/personal/sli/New_BacArena_project/new_breast_diet.txt",sep = '\t')

temp_meta <- subset(meta, ID == Ind)
temp_GEMs_amount <- all_amount_data[rownames(temp_meta)]
temp_GEMs_amount <- temp_GEMs_amount[temp_GEMs_amount != 0,,drop = FALSE]

ASV_GEM_map = read.table('~/personal/sli/New_BacArena_project/ASV_GEM_map_with_genus.txt',sep = '\t')

# define an empty arena object
ind_arena = Arena(n = 100,m = 100)

cons_GEMs_list = readRDS('~/personal/sli/New_BacArena_project/breast_all_GEMs_needed.rds')

add_bacs <- list()
for (i in names(cons_GEMs_list)) {
  add_bacs[[i]] <- Bac(model = cons_GEMs_list[[i]],limit_growth = TRUE,deathrate = 0.1,
                       minweight = runif(1,min = 0.2,max = 0.4),maxweight = runif(1,min = 0.8,max = 1),
                       growtype = "exponential")
}

# add selected personalized gut microbiota species to the defined empty arena object
for (asv in rownames(temp_GEMs_amount)) {
  temp_GEM_input <- ASV_GEM_map[ASV_GEM_map$ASV_ID == asv,]$renamed_gem
  ind_arena <- addOrg(object = ind_arena,
                      specI = add_bacs[[temp_GEM_input]],
                      amount = temp_GEMs_amount[asv,])
}

# add nutrients to the defined arena object
new_diet_ex <- as.character(Diet[,1])
new_diet_val <- as.numeric(Diet[,2])

for (met in 1:length(new_diet_ex)) {
  ind_arena <- addSubs(object = ind_arena,
                       smax = ifelse(new_diet_val[met]>1e-3,new_diet_val[met]/20,new_diet_val[met]),
                       mediac = new_diet_ex[met],
                       unit = 'mM',addAnyway = TRUE)
}

# run the arena simulation for 24 hours growth.
cat("Running simulation for", Ind, "...\n")
ind_eval = simEnv(ind_arena,time = 12,diffusion = TRUE)
dir.create("Ind_eval_results_no_treat_breast",showWarnings = FALSE)
saveRDS(ind_eval,file.path("Ind_eval_results_no_treat_breast",paste0(Ind,"_eval_result_breast.rds")))
cat("Finished simulation for",Ind,"\n")