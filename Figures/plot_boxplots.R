# Here, we compare the estimators from Nasri and Rémillard (2024) and
# the proposed estimator for zero-inflated continuous and discrete
# data, plotting the results as boxplots.

######################################################################
# SIMULATION BOXPLOTS
######################################################################
library(ggplot2)
library(latex2exp)

######################
# DEFINE PARAMETERS
######################

# Parameters for X-variable
lambda_F <- 2 # Poisson parameter
pi_F <- 0.5 # 1 - inflation at zero

# Parameters for Y-variable
lambda_G <- 2 # Poisson parameter
pi_G <- 0.5 # 1 - inflation at zero

# Copula
copula <- copula_frechet
cop_par <- 0.5

# Simulations
sim <- 1000 # Number of simulations

######################
# FIXED SAMPLE SIZE
######################

N <- 150 # Sample size

# Simulations
df_sim <- corr_bzip_main(lambda_F, lambda_G, pi_F, pi_G, cop_par, N, sim, copula)

# True value
rhoS_true <- calculate_Safari_rho_Spearman(lambda_F, lambda_G, pi_F, pi_G,
                                           cop_par, copula)

# Create boxplot
df_plot <- reshape2::melt(df_sim[,2:4], id.var=c())
ggplot(df_plot, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  geom_hline(yintercept=rhoS_true, linetype=2) +
  scale_fill_manual(breaks=c("cor", "Nasri", "Arends"),
                    values=c("#fbb4ae","#fbb4ae", "#b3cde3")) +
  labs(x="", y=TeX("$\\rho_S$ estimate")) +
  theme_classic() +
  theme(legend.position="none")

######################
# CHANGING SAMPLE SIZE
######################

# Initialise data frame with results
df_sim <- data.frame(matrix(nrow=0, ncol=5))
colnames(df_sim) <- c("N", "ID", "cor", "Nasri", "Arends")

for (N in c(50, 100, 500, 1000, 5000, 10000)) {
  
  # Simulations
  df_simN <- corr_bzip_main(lambda_F, lambda_G, pi_F, pi_G, cop_par, N, sim, copula)
  
  # Store results
  df_sim <- rbind(df_sim,
                  cbind(N, df_simN))
  
}

rhoS_true <- calculate_Safari_rho_Spearman(lambda_F, lambda_G, pi_F, pi_G,
                                           cop_par, copula)

# Create boxplot
df_sim$N <- as.factor(df_sim$N)
ggplot(df_sim, aes(x=N, y=Arends)) +
  geom_boxplot(fill="#b3cde3") +
  geom_hline(yintercept=rhoS_true, linetype=2) +
  labs(x="N", y=TeX("$\\hat{\\rho}_A$")) +
  theme_classic()
