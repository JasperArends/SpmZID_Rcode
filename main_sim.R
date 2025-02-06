######################################################################
# MAIN SIMULATION FUNCTIONS
######################################################################

# LOAD LIBRARIES
library(Metrics) # mse function

# Simulations
corr_bzip_main <- function(lambda1, lambda2, pi1, pi2, cop_par, N, sim, copula) {
  
  # Initialise data frame
  df_sim <- data.frame(matrix(nrow=sim, ncol=4))
  colnames(df_sim) <- c("ID", "cor", "Mesfioui", "Arends")
  
  # Perform Monte-Carlo simulations
  for (nx in 1:sim) {
    # Generate sample
    X <- gen_pois(lambda1, lambda2, pi1, pi2, cop_par, copula, N)
    
    # Make estimates
    df_sim[nx,] <- c(nx,
                     cor(X[,1], X[,2], method="spearman"),
                     spm_Mesfioui(X[,1], X[,2]),
                     spm_Arends(X[,1], X[,2]))
  }
  
  return (df_sim)
}

# Statistics on simulations
corr_bzip_stats <- function(lambda1, lambda2, pi1, pi2, cop_par, N, sim, copula, df_sim) {
  
  # Initialise data frame
  df_stat <- data.frame(matrix(nrow=3, ncol=6))
  colnames(df_stat) <- c("Method", "Mean", "Variance", "MSE", "Left_95CI", "Right_95CI")
  
  # Calculate true value
  spm_true <- calculate_Safari_rho_Spearman(lambda1, lambda2, pi1, pi2, cop_par, copula)
  
  # Calculate statistics
  df_stat[1,] <- c("cor",
                   mean(df_sim$cor),
                   var(df_sim$cor),
                   mse(spm_true, df_sim$cor),
                   quantile(df_sim$cor, 0.025)[[1]],
                   quantile(df_sim$cor, 0.975)[[1]])
  
  df_stat[2,] <- c("Mesfioui",
                   mean(df_sim$Mesfioui),
                   var(df_sim$Mesfioui),
                   mse(spm_true, df_sim$Mesfioui),
                   quantile(df_sim$Mesfioui, 0.025)[[1]],
                   quantile(df_sim$Mesfioui, 0.975)[[1]])
  
  df_stat[3,] <- c("Arends",
                   mean(df_sim$Arends),
                   var(df_sim$Arends),
                   mse(spm_true, df_sim$Arends),
                   quantile(df_sim$Arends, 0.025)[[1]],
                   quantile(df_sim$Arends, 0.975)[[1]])
  
  return (df_stat)
}
