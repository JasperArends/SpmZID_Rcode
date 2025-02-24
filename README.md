# Spearman's rho for zero-inflated data
Implementations based on Arends, J.R.M. and Perrone, E. (2025) *Spearman's rho for bi-variate zero-inflated data* [in preparation].

## Required libraries
`dplyr` `Metrics` `boot` `copula` `ggplot2` `latex2exp`

## Files
**Copulas.R**
Set of copulas that can be used for simulations.

**Estimation_functions.R**
Estimators for Spearman's rho adjusted to zero-inflated distributions.

**gen_sample.R**
Generates a sample of zero-inflated Poisson random variables for an arbitrary copula and zero-inflated exponential random variables.

**bounds.R**
Functions for the computation and estimatation of bounds on zero-inflated continuous and discrete random variables.

**DEMO.R**
Example of how the estimators can be used given a data set and how a simulation study can be initialised. The folder "Figures" contains the codes that were used to create the figures and tables from the article. These can be used in addition to the demo file for a more detailed simulation study.
