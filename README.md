# Spearman's rho for zero-inflated data
Implementations based on Arends, J.R.M. and Perrone, E. (2025) *Spearman's rho for bi-variate zero-inflated data*. In preparation.

## Required libraries
`dplyr` `Metrics` `boot` `copula` `ggplot2` `latex2exp`

## Files
The folder "Figures" contains the codes that were used to create the figures and tables from the article. These can be used in addition to the demo file for a more detailed simulation study.

**DEMO.R**
Example of how the estimators can be used given a data set and how a simulation study can be initialised.

**asymptotics.R**
Analyse the asymptotic behaviour or our estimator using bootstrap and a Monte-Carlo simulation.

**bounds.R**
Functions for the computation and estimatation of bounds on zero-inflated continuous and discrete random variables.

**copulas.R**
Set of copulas that can be used for simulations.

**gen_sample.R**
Generates a sample of zero-inflated Poisson random variables for an arbitrary copula and zero-inflated (exponential) random variables.

**main_estimators.R**
Estimators for Spearman's rho adjusted to zero-inflated distributions.

**main_sim.R**
Simulation functions.

## References
1. De Greef, N., van den Heuvel, E. and Zhuozhao, Z. (2020). Correlation estimation for bivariate zero-inflated discrete data. Master's thesis, Eindhoven University of Technology.
2. Mesfioui, M. and Trufin, J. (2017). Bounds on multivariate Kendall's tau and Spearman's rho for zero-inflated continuous variables and their application to insurance. *Methodology and Computing in Applied Probability*, 24, 1051-1059.
3. Safari-Katesari, H., Samadi, S.Y. and Zaroudi, S. (2020). Modelling count data via copulas. *Statistics*, 54(6), 1329-1355. doi: doi.org/10.1080/02331888.2020.1867140.
