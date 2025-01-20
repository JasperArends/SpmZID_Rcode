######################################################################
# DEMO
######################################################################

######################
# LOAD LIBRARIES
# AND DEPENDENCIES
######################
library(dplyr)
library(Metrics) # MSE

source("Copulas.R")
source("GenSample.R")
source("Estimator_functions.R")
source("bounds.R"


######################################################################
# ESTIMATOR COMPARISON - ZERO-INFLATED POISSON DISTRIBUTION
# We compare the performance of estimators using the Fréchet-copula
# for several parameter choices.
######################################################################

######################
# DEFINE PARAMETERS
######################
# Parameters for X-distribution
lambda_F <- 2 # Poisson parameter
pi_F <- 0.5 # Inflation at zero

# Parameters for Y-distribution
lambda_G <- 2 # Poisson parameter
pi_F <- 0.5 # Inflation at zero

# Copula
copula <- copula_frechet # See 'Copulas.R' file for available copulas
cop_par <- 0.5

# Simulations
N <- 150 # Sample size
sim <- 300 # Number of simulations

######################
# PERFORM SIMULATIONS
######################
df_sim <- corr_bzip_main(lambda1, lambda2, pi1, pi2, cop_par, N, sim, copula) # Simulations
df_stat <- corr_bzip_stats(lambda1, lambda2, pi1, pi2, cop_par, N, sim, copula, df_sim) # Statistics

# Print results
head(df_stat)


######################################################################
# APPLICATION
# For the data stored in `A` and `B`, the estimators are computed and
# compared. The proposed estimator is then adjusted using the
# estimates of the bounds. 
######################################################################

######################
# LOAD DATA
######################
# Example of how data should be read into the variables `A` and `B`.
# The data is loaded as the .csv file `data.csv`.
# df_data <- read.csv('data.csv')
# A <- df_data[,1]
# B <- df_data[,2]

######################
# COMPUTE ESTIMATES
######################
spm_estimates <- list(
  spm_Arends(A, B), # Proposed estimator
  spm_Mesfioui(A, B), # Estimator by Mesfioui and Trufin (2022)
  cor(A, B, method="spearman")
)
names(spm_estimates) <- c("Arends", "Mesfioui", "cor")

######################
# DATA STATISTICS
######################
N <- length(A) # Sample size
p1 <- sum(A == 0) / N # Probability of 0 in A
p2 <- sum(B == 0) / N # Probability of 0 in B

######################
# COMPUTE BOUNDS
# (upper and lower
# bounds respectively)
######################
bnd_cont <- spm_bounds_cont_Arends_est(A, B) # Zero-inflated continuous data
bnd_disc <- spm_bounds_disc_Arends_est(A, B) # Zero-inflated discrete data

######################
# COMPUTE ADJUSTED
# ASSOCIATION MEASURE
######################
# Zero-inflated continuous data
if (spm_estimates$Arends > 0)
  spm_estimates$Arends / bnd_cont[1]
else
  spm_estimates$Arends / bnd_cont[2]

# Zero-inflated discrete data
if (spm_estimates$Arends > 0)
  spm_estimates$Arends / bnd_disc[1]
else
  spm_estimates$Arends / bnd_disc[2]
