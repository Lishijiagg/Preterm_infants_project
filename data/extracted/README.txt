Per-infant tables produced by extraction/extract_eval.R from the BacArena Eval objects.

One directory per simulation condition, one .rds per infant, plus _mediac.rds holding the name
mapping for the medium components. Each .rds contains four tables: end-point fermentation product
concentrations, the HMO time course, the growth curve of every organism, and the exchange flux of
every organism at every time step.

extracted_no_treat_breast    baseline, breast milk medium, 51 infants
extracted_no_treat_formula   baseline, formula medium, 51 infants
extracted_pro_breast         seven probiotic strains, breast milk medium, 40 preterm infants
extracted_pro_formula        seven probiotic strains, formula medium, 40 preterm infants
extracted_supp_breast        personalised nutrient supplementation, breast milk medium
extracted_supp_formula       personalised nutrient supplementation, formula medium
extracted_syn_breast         nutrients and probiotics combined, breast milk medium
extracted_syn_formula        nutrients and probiotics combined, formula medium
extracted_hmo_formula        HMO added to the formula medium
extracted_hmo_pro_formula    HMO and probiotics added to the formula medium

These are the direct input to GEMs-BacArena_analysis.R. The Eval objects they derive from are not
distributed because a single object exceeds 100 GB.
