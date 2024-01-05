# ProbUncertainty

R package to calculate uncertainty estimates of the difference in probability 

This R package provides functions to calculate
(1) the standard error of a difference between two probabilities using the delta method from a binomial GLM output, 
(2) profile likelihood-based confidence intervals of the difference in probability,
(3) score confidence intervals of the difference in probability

Please cite: "Classical test, linear models, and their extensions for the analysis of 2x2 contingency tables"

To install:
```
library(devtools)
install_github("https://github.com/rebebba/ProbUncertainty.git", build_vignettes = TRUE)

library(ProbUncertainity)
vignette("Intro_to_ProbUncertainity", package="ProbUncertainity")
```
