######################################################################
# ESTIMATIONG OF VARIANCE
######################################################################
library(boot)

source("main_estimators.R")
source("main_sim.R")

######################
# Monte-Carlo
# simulations
######################
# Estimates the variance of all estimators for Spearman's rho based on `sim` simulations

spm_var.MC <- function(lambda1, lambda2, pi1, pi2, cop_par, N, sim, copula) {
  df_sim <- corr_bzip_main(lambda1, lambda2, pi1, pi2, cop_par, N, sim, copula)
  var_est <- sapply(df_sim[,2:4], function(x) var(x))
  
  return (var_est)
}


######################
# Bootstrap
######################

spm_var.boot <- function(lambda1, lambda2, pi1, pi2, cop_par, N, sim, copula, R) {
  
  spm_Arends_boot <- function(X, d) {
    return ( spm_Arends(X[d,1], X[d,2]) )
  }
  
  # Draw sample
  X <- gen_pois(lambda1, lambda2, pi1, pi2, cop_par, copula, N)
  
  # Perform bootstrap
  rhoS_boot <- boot(X, statistic = spm_Arends_boot, R = R)$t
  
  return ( var(rhoS_boot) )
}
