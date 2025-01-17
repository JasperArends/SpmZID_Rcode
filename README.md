# Spearman's rho for zero-inflated data
Implementations based on Arends, J.R.M. and Perrone, E. (n.d.) *Spearman's rho for bi-variate zero-inflated data*.

## Required libraries
`dplyr` `Metrics` `boot` `copula` `ggplot2` `latex2exp`

## Files
**Copulas.R**
Set of useful copulas that can be used for simulations.

**Estimation_functions.R**
Estimators for Spearman's rho adjusted to zero-inflated distributions.

**GenSample.R**
Generates a sample of zero-inflated Poisson random variables for an arbitrary copula and zero-inflated exponential random variables.

**Bound.R**
Functions for the computation and estimatation of bounds on zero-inflated continuous and discrete random variables.
