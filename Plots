######################################################################
# SPEARMAN'S RHO FOR BI-VARIATE ZERO-INFLATED DATA
# PLOTS
######################################################################
library(reshape2)
library(ggplot)
library(latex2exp)


######################################################################
# BOUNDS FOR ZERO-INFLATED CONTINUOUS DISTRIBUTIONS
######################################################################
source("Bounds.R")

p2 <- 0.5 # Fix p2, vary p1

# Compute upper bounds
df_upper <- data.frame(p1=seq(0, 1, 0.02), Arends=p2, Mesfioui=0)
for (px in 1:nrow(df_upper)) {
  df_upper$Arends[px] <- spm_bounds_Arends(1 - df_upper$p1[px], 1 - p2)[1]
  df_upper$Mesfioui[px] <- spm_bounds_Mesfioui(1 - df_upper$p1[px], 1 - p2)[1]
}

# Compute lower bounds
df_lower <- data.frame(p1=seq(0, 1, 0.02), Arends=p2, Mesfioui=0)
for (px in 1:nrow(df_upper)) {
  df_lower$Arends[px] <- spm_bounds_Arends(1 - df_lower$p1[px], 1 - p2)[2]
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
