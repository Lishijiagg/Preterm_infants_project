# Extras

Analyses that were explored but are not reported in the manuscript. They are kept for completeness.
They do not support any conclusion in the manuscript and have not been checked to the same standard
as the code in the analysis script.

sensitivity_03_reseed.R tests the stochastic stability of the flux attribution under repeated
simulation. It was never run. The within-genus screen in Part 7.3 of the analysis script found that
the candidate strains carry almost identical reaction repertoires, so re-simulation was not
warranted.

supplement_profile_permanova.R applies PERMANOVA to the nutrient requirement profiles of the two
preterm groups. The apparent group difference traced to a single model,
ASV004_Serratia_marcescens_Db11, which contributes a median of 81 substrates on its own in 22 of the
40 communities. The test measures a property of that reconstruction rather than a biological
difference, so it was removed from the Methods.
