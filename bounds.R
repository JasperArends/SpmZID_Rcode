######################################################################
# BOUNDS ON SPEARMAN'S RHO
# for zero-inflated data
######################################################################
source("main_estimators.R")
source("copulas.R")

######################################################################
# ZERO-INFLATED CONTINUOUS DISTRIBUTIONS
# Based on Arends, J.R.M. and Perrone, E. (n.d.). Spearman's rho for
# zero inflated data.
######################################################################

# True value
spm_bounds_cont_Arends <- function(pi1, pi2) {
  p1 <- 1 - pi1
  p2 <- 1 - pi2
  
  pm <- max(p1, p2)
  upper <- (1 - pm)^3 + 3 * pm * (1 - pm)
  
  if (1 - p1 - p2 < 0) {
    lower <- -3 * (1 - p1) * (1 - p2)
  } else
    lower <- -(1 - p1 - p2)^3 - 3 * (1 - p1 - p2 + p1 * p2) * (p1 + p2)
  return ( c(upper, lower) )
}

# Estimator version
spm_bounds_cont_Arends_est <- function(A, B) {
  pi1 <- sum(A > 0)/length(A)
  pi2 <- sum(B > 0)/length(B)
  
  return ( spm_bounds_Arends(pi1, pi2) )
}

######################################################################
# APPROXIMATE BOUNDS FOR FULLY DISCRETE DATA
# Based on Arends, J.R.M. and Perrone, E. (n.d.). Spearman's rho for
# zero inflated data.
######################################################################
spm_bound_disc_Arends <- function(lambda1, lambda2, pi1, pi2) {
  return ( c(spm_up_bound_disc_Arends(lambda1, lambda2, pi1, pi2),
             spm_lw_bound_disc_Arends(lambda1, lambda2, pi1, pi2)) )
}

spm_bound_disc_Arends_est <- function(A, B) {
  return ( spm_up_bound_disc_Arends_est(A, B),
           spm_lw_bound_disc_Arends_est(A, B) )
}

spm_up_bound_disc_Arends <- function(lambda1, lambda2, pi1, pi2) {
  p1 <- (1 - pi1) + pi1 * dpois(0, lambda1) # P(X = 0)
  p2 <- (1 - pi2) + pi2 * dpois(0, lambda2) # P(Y = 0)
  
  # Assume p1 <= p2
  if (p1 > p2) { # Swap distributions
    lambdac <- lambda1
    lambda1 <- lambda2
    lambda2 <- lambdac
    
    pic <- pi1
    pi1 <- pi2
    pi2 <- pic
    
    p1c <- p1
    p1 <- p2
    p2 <- p1c
  }
  
  if (max(p1, p2) == 1) # Only zeros occur for one variable
    return ( 0 )
  else {
    
    # Find s_tilde
    s <- 0
    while (TRUE) {
      Fs <- (1 - pi1) + pi1 * ppois(s, lambda1)
      if (Fs > p2)
        break
      s <- s + 1
    }
    Fsm <- (1 - pi1) + pi1 * ppois(s - 1, lambda1)
    
    # Find u_tilde
    u <- 0
    while (TRUE) {
      Gu <- (1 - pi2) + pi2 * ppois(u, lambda2)
      if (Gu > Fs)
        break
      u <- u + 1
    }
    Gum <- (1 - pi2) + pi2 * ppois(u - 1, lambda2)
    
    rho_S11 <- 1
    rhomax <- (1 - p2)^3 * rho_S11 + 3 * (1 - p2) * p2 +
      3 * (p2 - Fsm) * (Fs * (p2 - Gu - Gum) + Gu * Gum)
    
    return (rhomax)
  }
}

spm_lw_bound_disc_Arends <- function(lambda1, lambda2, pi1, pi2) {
  p1 <- (1 - pi1) + pi1 * dpois(0, lambda1) # P(X = 0)
  p2 <- (1 - pi2) + pi2 * dpois(0, lambda2) # P(Y = 0)
  
  if (1 - p1 - p2 <= 0) {
    rhomin <- -3 * (1 - p1) * (1 - p2)
  } else {
    p11 <- 1 - p1 - p2
    
    # Find s'
    s <- 0
    while (TRUE) {
      Fs <- (1 - pi1) + pi1 * ppois(s, lambda1)
      if (Fs + p2 - 1 > 0)
        break
      s <- s + 1
    }
    Fsm <- (1 - pi1) + pi1 * ppois(s - 1, lambda1)
    
    # Find u'
    u <- 0
    Gum <- 0
    while (TRUE) {
      Gu <- (1 -  pi2) + pi2 * ppois(u, lambda2)
      if (Fsm + Gu - 1 > 0 | Gum == 1)
        break
      Gum <- Gu
      u <- u + 1
    }
    
    # Find t'
    t <- 0
    while (TRUE) {
      Gt <- (1 - pi2) + pi2 * ppois(t, lambda2)
      if (Gt + p1 - 1 > 0)
        break
      t <- t + 1
    }
    Gtm <- (1 - pi2) + pi2 * ppois(t - 1, lambda2)
    
    # Find v'
    v <- 0
    Fvm <- 0
    while (TRUE) {
      Fv <- (1 - pi1) + pi1 * ppois(v, lambda1)
      if (Fv + Gtm - 1 > 0 | Fvm == 1)
        break
      Fvm <- Fv
      v <- v + 1
    }
    W <- c(1 - p1 - p2,
           1 - Fsm - Gtm,
           1 - Fs - p2,
           1 - Fsm - p2,
           1 - p1 - Gt,
           1 - p1 - Gtm,
           1 - Fsm - Gum,
           1 - Fvm - Gtm)
    
    if (s == v) {
      I1 <- 1 - p2
      I2 <- 1 - p1
    } else {
      I1 <- Fv
      I2 <- Gu
    }
    
    rhomin <- W[1]^3 * (-1) + 3 * W[1] * (p1 * p2 - p1 - p2) - 3 * p1 * p2 +
      -3 * W[3] * ( p2 * W[4] + (Gum - p2)^2 + W[7] * (I2 - 2 * p2 + Gum) ) +
      -3 * W[5] * ( p1 * W[6] + (Fvm - p1)^2 + W[8] * (I1 - 2 * p1 + Fvm) )
  }
  
  return (rhomin)
}

spm_up_bound_disc_Arends_est <- function(A, B) {
  N <- length(A)
  
  p1 <- sum(A == 0)/N
  p2 <- sum(B == 0)/N
  
  # Assume p1 <= p2
  if (p1 > p2) { # Swap distributions
    C <- A
    A <- B
    B <- C
    
    p1c <- p1
    p1 <- p2
    p2 <- p1c
  }
  
  if (max(p1, p2) == 1) # Only zeros occur for one variable
    return ( 0 )
  else {
    
    Fd <- cumsum( as.vector(table(A)) ) /N
    Gd <- cumsum( as.vector(table(B)) ) /N
    
    # Find s_tilde
    s <- which(Fd > p2)[1]
    Fs <- Fd[s]
    Fsm <- Fd[s - 1]
    
    # Find u_tilde
    u <- which(Gd > Fs)[1]
    Gu <- Gd[u]
    Gum <- Gd[u - 1]
    
    rho_S11 <- 1
    rhomax <- (1 - p2)^3 * rho_S11 + 3 * (1 - p2) * p2 +
      3 * (p2 - Fsm) * (Fs * (p2 - Gu - Gum) + Gu * Gum)
    
    return (rhomax)
  }
}

spm_lw_bound_disc_Arends_est <- function(A, B) {
  N <- length(A)
  
  p1 <- sum(A == 0)/N
  p2 <- sum(B == 0)/N
  
  if (1 - p1 - p2 <= 0) {
    rhomin <- -3 * (1 - p1) * (1 - p2)
  } else {
    p11 <- 1 - p1 - p2
    
    Fd <- cumsum( as.vector(table(A)) ) /N
    Gd <- cumsum( as.vector(table(B)) ) /N
    
    # Find s_tilde'
    s <- which(Fd + p2 - 1 > 0)[1]
    Fs <- Fd[s]
    Fsm <- Fd[s - 1]
    
    # Find u_tilde'
    u <- which(Fsm + Gd - 1 > 0)[1]
    Gu <- Gd[u]
    Gum <- Gd[u - 1]
    
    # Find t_tilde'
    t <- which(Gd + p1 - 1 > 0)[1]
    Gt <- Gd[t]
    Gtm <- Gd[t - 1]
    
    # Find v_tilde'
    v <- which(Fd + Gtm - 1 > 0)[1]
    Fv <- Fd[v]
    Fvm <- Fd[v - 1]
    
    W <- c(1 - p1 - p2,
           1 - Fsm - Gtm,
           1 - Fs - p2,
           1 - Fsm - p2,
           1 - p1 - Gt,
           1 - p1 - Gtm,
           1 - Fsm - Gum,
           1 - Fvm - Gtm)
    
    if (s == v) {
      I1 <- 1 - p2
      I2 <- 1 - p1
    } else {
      I1 <- Fv
      I2 <- Gu
    }
    
    rhomin <- W[1]^3 * (-1) + 3 * W[1] * (p1 * p2 - p1 - p2) - 3 * p1 * p2 +
      -3 * W[3] * ( p2 * W[4] + (Gum - p2)^2 + W[7] * (I2 - 2 * p2 + Gum) ) +
      -3 * W[5] * ( p1 * W[6] + (Fvm - p1)^2 + W[8] * (I1 - 2 * p1 + Fvm) )
  }
  
  return (rhomin)
}

######################################################################
# ZERO-INFLATED CONTINUOUS DISTRIBUTIONS
# Based on Mesfioui, M., Trufin, J., and Zuyderhoff, P. (2022). Bounds
# on Spearman's rho when at least one random variable is discrete.
# European Actuarial Journal, 12, 321-348.
######################################################################
spm_bounds_Mesfioui <- function(pi1, pi2) {
  p1 <- 1 - pi1
  p2 <- 1 - pi2
  v <- sqrt((1 - p1^3) * (1 - p2^3))
  return ( c(v, -v))
}

spm_bounds_Mesfioui_est <- function(A, B) {
  pi1 <- sum(A > 0)/length(A)
  pi2 <- sum(B > 0)/length(B)
  
  return ( spm_bounds_Mesfioui(pi1, pi2))
}

######################################################################
# TRUE BOUNDS FOR DISCRETE DATA
# Based on Safari-Katesari, H., Samadi, S.Y. and Zaroudi, S. (2020).
# Modelling count data via copulas. Statistics, 54(6), 1329-1355.
# doi.org/10.1080/02331888.2020.1867140.
######################################################################
spm_up_bound_Safari <- function(lambda1, lambda2, pi1, pi2)
  calculate_Safari_rho_Spearman(lambda1, lambda2, pi1, pi2, 0, copula_upper)

spm_lw_bound_Safari <- function(lambda1, lambda2, pi1, pi2)
  calculate_Safari_rho_Spearman(lambda1, lambda2, pi1, pi2, 0, copula_lower)
