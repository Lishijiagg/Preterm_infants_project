# ================================================================================================================== #
# Build the HMO-supplemented formula medium for the synbiotic arm.
#
# The formula medium carries no human milk oligosaccharide. The thirteen HMO exchange reactions of the breast medium
# are transferred onto it at the same concentrations, so that the two media differ in their oligosaccharide content
# and in nothing else. Together with the existing formula and formula-plus-probiotic runs this gives a 2 x 2 design
# in which the prebiotic and the strains able to use it can be tested for synergy.
#
# Run once in RStudio; writes new_formula_hmo_diet.txt.
# ================================================================================================================== #
setwd("~/personal/sli/New_BacArena_project/")

new_breast  = read.table("new_breast_diet.txt",  sep = '\t')
new_formula = read.table("new_formula_diet.txt", sep = '\t')

# Every oligosaccharide exchange the breast medium carries, not only the six quantified downstream
hmo_pattern = paste(c("2fuclac","3fuclac","lacnttr","lacnnttr","3slac","6slac",
                      "dfuclac","fuc_L","lacdfucttr","lacnfucpt_i","lacnhx",
                      "lntri_ii","neulacnttr_b"), collapse = "|")
hmo_rows = new_breast[grepl(hmo_pattern, new_breast[[1]]), ]
cat("HMO exchanges carried over from the breast medium: ", nrow(hmo_rows), "\n", sep = '')
print(hmo_rows[, 1:2])

already = intersect(hmo_rows[[1]], new_formula[[1]])
if (length(already) > 0) {
  cat("Already present in the formula medium and therefore not duplicated: ",
      paste(already, collapse = ", "), "\n", sep = '')
  hmo_rows = hmo_rows[!hmo_rows[[1]] %in% already, ]
}

new_formula_hmo = rbind(new_formula, hmo_rows)
cat("Formula medium: ", nrow(new_formula), " exchanges; with HMO: ",
    nrow(new_formula_hmo), "\n", sep = '')

write.table(new_formula_hmo, file = "new_formula_hmo_diet.txt",
            row.names = TRUE, col.names = TRUE, sep = '\t')
