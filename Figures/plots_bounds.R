######################################################################
# PLOTS REGARDING THE BOUNDS ON SPEARMAN'S RHO
######################################################################
source("Bounds.R")

######################################################################
# COMPARISON BETWEEN BOUNDS FOR ZERO-INFLATED CONTINUOUS DISTRIBUTIONS
# Fix one level of inflation and vary the other, compare the bounds
# from Mesfioui and Trufin (2017) and the proposed ones in Proposition
# 4.1.
######################################################################

# Fix one inflation parameter
p2 <- 0.5 #

# Compute upper bounds
df_upper <- data.frame(p1=seq(0, 1, 0.02), Arends=p2, Mesfioui=0)
for (px in 1:nrow(df_upper)) {
  df_upper$Arends[px] <- spm_bounds_cont_Arends(1 - df_upper$p1[px], 1 - p2)[1]
  df_upper$Mesfioui[px] <- spm_bounds_Mesfioui(1 - df_upper$p1[px], 1 - p2)[1]
}

# Compute lower bounds
df_lower <- data.frame(p1=seq(0, 1, 0.02), Arends=p2, Mesfioui=0)
for (px in 1:nrow(df_upper)) {
  df_lower$Arends[px] <- spm_bounds_cont_Arends(1 - df_lower$p1[px], 1 - p2)[2]
  df_lower$Mesfioui[px] <- spm_bounds_Mesfioui(1 - df_lower$p1[px], 1 - p2)[2]
}

# Plot upper bounds
df_plot <- reshape2::melt(df_upper, id="p1")
ggplot(df_plot, aes(x=p1, y=value,
                    col=variable, linetype=variable)) +
  geom_line() +
  scale_color_manual(breaks=c("Arends", "Mesfioui"), values=c("#C81919","black")) +
  scale_linetype_manual(breaks=c("Arends", "Mesfioui"), values=c(2, 1)) +
  labs(x=TeX("$1 - \\pi_F$"), y=TeX("$\\rho_{max}$")) +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.position="none")

# Plot lower bounds
df_plot <- reshape2::melt(df_lower, id="p1")
ggplot(df_plot, aes(x=p1, y=value,
                    col=variable, linetype=variable)) +
  geom_line() +
  scale_color_manual(breaks=c("Arends", "Mesfioui"), values=c("#C81919","black")) +
  scale_linetype_manual(breaks=c("Arends", "Mesfioui"), values=c(2, 1)) +
  labs(x=TeX("$1 - \\pi_F$"), y=TeX("$\\rho_{min}$")) +
  ylim(-1, 0) +
  theme_classic() +
  theme(legend.position="none")

######################################################################
# BOUNDS ESTIMATION FOR ZERO-INFLATED CONTINUOUS DATA
# Estimates of the bounds for zero-inflated continuous distributions
# with Fréchet-copula and fixed copula-parameter. Estimates are
# plotted as means with a 95% confidence interval.
######################################################################

# Simulation parameters
N <- 150 # Sample size
sim <- 1000 # Number of simulations
p2 <- 0.5 # Inflation at zero
alpha <- 0.5 # Copula parameter

# Initialise data frame for results
df_bounds_ci <- data.frame(matrix(nrow=0, ncol=9))
colnames(df_bounds_ci) <- 
  c("p1", "max_true", "min_true",
    "max_avg", "max_lw", "max_up",
    "min_avg", "min_lw", "min_up")

# Compute / estimate upper bounds
for (p1 in seq(0, 1, 0.01)) {
  df_bounds <- data.frame(matrix(nrow=sim, ncol=2))
  colnames(df_bounds) <- c("upper", "lower")
  for (nx in 1:sim) {
    # Draw sample
    U <- gen_Frechet(N, alpha)
    
    # Estimate probabilities of inflation at zero
    p1_est <- sum(U[,1] <= p1)/N
    if (is.na(p1_est))
      p1_est <- 0
    p2_est <- sum(U[,2] <= p2)/N
    
    # Compute corresponding bounds
    df_bounds[nx,] <- spm_bounds_Arends(1 - p1_est, 1 - p2_est)
  }
  
  df_bounds_ci[nrow(df_bounds_ci)+1,] <- 
    c(pi1,
      spm_bounds_Arends(1 - p1, 1 - p2),
      mean(df_bounds$upper),
      quantile(df_bounds$upper, c(0.025, 0.975)),
      mean(df_bounds$lower),
      quantile(df_bounds$lower, c(0.025, 0.975)))
}

ggplot(df_bounds_ci) +
  geom_ribbon(aes(x=p1, ymin=max_lw, ymax=max_up), alpha=.15) +
  geom_line(aes(x=p1, y=max_avg), col="black") +
  geom_line(aes(x=p1, y=max_true), col="#C81919", linetype=2) +
  theme_classic() +
  labs(x=TeX("$1 - \\pi_F"), y=TeX("$\\rho_{max}")) +
  ylim(0, 1) +
  theme(axis.line = element_line(size=.3),
        axis.ticks = element_line(size=.3))

ggplot(df_bounds_ci) +
  geom_ribbon(aes(x=p1, ymin=min_lw, ymax=min_up), alpha=.15) +
  geom_line(aes(x=p1, y=min_avg), col="black") +
  geom_line(aes(x=p1, y=min_true), col="#C81919", linetype=2) +
  theme_classic() +
  labs(x=TeX("$1 - \\pi_F"), y=TeX("$\\rho_{min}")) +
  ylim(-1, 0) +
  theme(axis.line = element_line(size=.3),
        axis.ticks = element_line(size=.3))


######################################################################
# BOXPLOTS OF BOUNDS FOR SPECIFIC VALUES OF INFLATION PROBABILITIES
# Simulations to estimate the bounds for zero-inflated continuous
# distributions for fixed values of p1 and p2, and varying sample
# sizes.
######################################################################

# Simulation parameters
p1 <- 0.5 # Inflation at zero for X-variable
p2 <- 0.5 # Inflation at zero for Y-variable
sim <- 1000 # Number of simulations

# Initialise data frames for results
df_bounds_est <- data.frame(matrix(nrow=0, ncol=3))
colnames(df_bounds_est) <- c("N", "upper", "lower")

for (N in c(50, 100, 500, 1000, 5000, 10000)) {
  for (nx in 1:sim) {
    # Draw sample
    U <- gen_Frechet(N, alpha)
    
    # Estimate probabilities of inflation
    p1_est <- sum(U[,2] <= p1)/N
    p2_est <- sum(U[,2] <= p2)/N
    
    # Estimate bounds
    df_bounds_est[nrow(df_bounds_est)+1,] <- 
      c(N,
        spm_bounds_Arends(1-p1_est, 1-p2_est))
  }
}

# Boxplot for upper bound
upper_true <- spm_bounds_Arends(1 - p1, 1 - p2)[1] # True value
df_bounds_est$N <- as.factor(df_bounds_est$N)

ggplot(df_bounds_est) +
  geom_boxplot(aes(x=N, y=upper), lwd=.3, outlier.size=.2, fill='grey') +
  geom_hline(yintercept=upper_true, linetype=2, lwd=.3) +
  labs(x="N", y=TeX("$\\rho_{max}$")) +
  theme_classic() +
  theme(axis.line = element_line(size=.3),
        axis.ticks = element_line(size=.3),
        plot.title = element_text(hjust=.5, size=10)) +
  ylim(0, 1)

# Boxplot for lower bound
lower_true <- spm_bounds_Arends(1 - p1, 1 - p2)[2] # True value
df_bounds_est$N <- as.factor(df_bounds_est$N)

ggplot(df_bounds_est) +
  geom_boxplot(aes(x=N, y=lower), lwd=.3, outlier.size=.2, fill='grey') +
  geom_hline(yintercept=lower_true, linetype=2, lwd=.3) +
  theme_classic() +
  theme(axis.line = element_line(size=.3),
        axis.ticks = element_line(size=.3),
        plot.title = element_text(hjust=.5, size=10)) +
  labs(x="N", y=TeX("$\\rho_{min}$")) +
  ylim(-1, 0)
