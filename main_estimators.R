######################################################################
# SPEARMAN'S RHO ESTIMATORS
# Adjusted for zero-inflated data
######################################################################
library(dplyr)

######################################################################
# ESTIMATOR rho_M
# Based on Mesfiou, M. and Trufin, J. (2017). Bounds on multivariate
# Kendall's tau and Spearman's rho for zero-inflated continuous
# variables and their application to insurance. Methodology and
# Computing in Applied Probability, 24, 1051-1059.
######################################################################
spm_Mesfioui <- function(A, B) {
  N <- length(A)
  
  p1 <- sum(A == 0) / N
  p2 <- sum(B == 0) / N
  
  Q <- 0
  
  for (nx in 1:N) {
    X <- A[nx]
    Y <- B[nx]
    
    F1 <- sum(A[A > 0] <= X) / sum(A > 0)
    F2 <- sum(B[B > 0] <= Y) / sum(B > 0)
    
    xi1 <- p1 * (X == 0) + 2 * (p1 + (1 - p1) * F1) * (X > 0) - 1
    xi2 <- p2 * (Y == 0) + 2 * (p2 + (1 - p2) * F2) * (Y > 0) - 1
    
    Q <- Q + (1 + xi1) * (1 + xi2) + (1 - xi1) * (1 - xi2)
  }
  
  Q <- 6 * Q /( 4 * N ) - 3
  
  return(Q)
}


######################################################################
# ESTIMATOR rho_A
# Based on Arends, J.R.M. and Perrone, E. (n.d.). Spearman's rho for
# zero inflated data.
######################################################################
spm_Arends <- function(A, B) {
  
  N <- length(A)
  
  p00 <- sum(A == 0 & B == 0) / N # P(X = 0, Y = 0)
  p01 <- sum(A == 0 & B > 0) / N  # P(X = 0, Y > 0)
  p10 <- sum(A > 0 & B == 0) / N  # P(X > 0, Y = 0)
  p11 <- sum(A > 0 & B > 0) / N   # P(x > 0, Y > 0)
  
  ######################################
  # Margins
  p0m <- p00 + p01 # P(X = 0)
  p1m <- p10 + p11 # P(X > 0)
  pm0 <- p00 + p10 # P(Y = 0)
  pm1 <- p01 + p11 # P(Y > 0)
  
  
  ######################################
  # Split data
  X10 <- which(A > 0 & B == 0)
  X01 <- which(A == 0 & B > 0)
  X11 <- which(A > 0 & B > 0)
  
  n11 <- length(X11)
  n10 <- length(X10)
  n01 <- length(X01)
  
  # Initialise probabilities
  p1s <- 0 # P(X10 > X11)
  p1c <- 0 # P(X10 = X11)
  p2s <- 0 # P(Y01 > Y11)
  p2c <- 0 # P(Y01 = Y11)
  
  for (ix in X10) {
    p1s <- p1s + sum(A[ix] > A[X11])
    p1c <- p1c + sum(A[ix] == A[X11])
  }
  p1s <- p1s / max(n10 * n11, 1)
  p1c <- p1c / max(n10 * n11, 1)
  
  for (ix in X01) {
    p2s <- p2s + sum(B[ix] > B[X11])
    p2c <- p2c + sum(B[ix] == B[X11])
  }
  p2s <- p2s / max(n01 * n11, 1)
  p2c <- p2c / max(n01 * n11, 1)
  
  D <- p11 * (p10 * (1 - 2*p1s - p1c) + p01 * (1 - 2*p2s - p2c))
  
  ######################################
  # Estimate P(C | C1111) - P(D | C1111)
  # Probabilities are stored in pC1, pC2, pC3 and pC4 such that pC1[1]
  # and pC1[2] correspond to concordance, discordance and ties for the
  # first case, where X3 > 0 and Y2 > 0 (similarly for pC2, pC3 and pC4).
  pC2 <- c(0, 0) # X3 = 0, Y2 > 0
  pC3 <- c(0, 0) # X3 > 0, Y2 = 0
  pC4 <- c(0, 0) # X3 = 0, Y2 = 0
  
  nX <- c(0, 0, 0, 0)
  nY <- c(0, 0, 0, 0)
  nXdb <- c(0, 0)
  for (ix in X11) {
    X1 <- A[ix]
    Y1 <- B[ix]
    
    nX[1] <- sum(X1 > A[X11])
    nX[2] <- sum(X1 > A[X10])
    nX[3] <- sum(X1 < A[X11])
    nX[4] <- sum(X1 < A[X10])
    nY[1] <- sum(Y1 > B[X11])
    nY[2] <- sum(Y1 > B[X01])
    nY[3] <- sum(Y1 < B[X11])
    nY[4] <- sum(Y1 < B[X01])
    
    pC2[1] <- pC2[1] + nX[1] * nY[2] + nX[3] * nY[4]
    pC2[2] <- pC2[2] + nX[1] * nY[4] + nX[3] * nY[2]
    
    pC3[1] <- pC3[1] + nX[2] * nY[1] + nX[4] * nY[3]
    pC3[2] <- pC3[2] + nX[2] * nY[3] + nX[4] * nY[1]
    
    pC4[1] <- pC4[1] + nX[2] * nY[2] + nX[4] * nY[4]
    pC4[2] <- pC4[2] + nX[2] * nY[4] + nX[4] * nY[2]
  }
  
  # Concordance - discordance
  pC <- c((pC2[1] - pC2[2]),
          (pC3[1] - pC3[2]),
          (pC4[1] - pC4[2]))
  # Convert to probabilities
  cnt <- c(as.double(n11) * (n11 - 1) * n01,
           as.double(n11) * n10 * (n11 - 1),
           as.double(n11) * n10 * n01)
  pC <- pC / pmax(cnt, 1) # Avoid division by 0
  
  rho11 <- cor(A[X11], B[X11], method="spearman")/3
  if (is.na(rho11))
    rho11 <- 0
  
  pPos <- p11 * (p11^2 * rho11 + p01 * p11 * pC[1] +
                   p11 * p10 * pC[2] + p01 * p10 * pC[3])
  
  # Spearman's rho estimate
  rho_est <- 3 * pPos + 3 * (p00 * p11 - p10 * p01) + 3*D
  
  return (rho_est)
}


######################################################################
# TRUE VALUE OF SPEARMAN'S RHO FOR DISCRETE DATA
# Based on Safari-Katesari, H., Samadi, S.Y. and Zaroudi, S. (2020).
# Modelling count data via copulas. Statistics, 54(6), 1329-1355.
# doi.org/10.1080/02331888.2020.1867140.
# Implementation by de Greef, N., van den Heuvel, E. and Zhuozhao, Z.
# (2020). Correlation estimation for bivariate zero-inflated discrete
# data. Master's thesis, Eindhoven University of Technology
######################################################################
calculate_Safari_rho_Spearman <- function(lambda1, lambda2, pi1, pi2, cop_par, copula){
  y1 <- 0
  y2 <- 0
  part1 <- 0
  part2 <- 0
  nruns <- 50
  for(y1 in 0:nruns){
    f1y1 <- (y1 == 0) * ((1 - pi1) + pi1 * exp(-lambda1)) + (y1 > 0) * (pi1 * exp(-lambda1 + log(lambda1)*y1 - lgamma(y1 + 1) ))
    f2y1 <- (y1 == 0) * ((1 - pi2) + pi2 * exp(-lambda2)) + (y1 > 0) * (pi2 * exp(-lambda2 + log(lambda2)*y1 - lgamma(y1 + 1)))
    for(y2 in 0:nruns){
      f2y2 <- (y2 == 0) * ((1 - pi2) + pi2 * exp(-lambda2)) + (y2 > 0) * (pi2 * (exp(-lambda2 + log(lambda2)*y2 - lgamma(y2 + 1))))
      F1 <- 1 - pi1 + pi1 * ppois(y1, lambda1)
      F2 <- 1 - pi2 + pi2 * ppois(y2, lambda2)
      F1m <- (y1 > 0) * (1 - pi1 + pi1 * ppois(y1 - 1, lambda1))
      F2m <- (y2 > 0) * (1 - pi2 + pi2 * ppois(y2 - 1, lambda2))
      h <- copula(F1, F2, cop_par) - copula(F1m, F2, cop_par) - copula(F1, F2m, cop_par) + copula(F1m, F2m, cop_par)
      part1 <- part1 + h * ( (1-F1) * (1-F2) + F1m * F2m - (1/2) * f1y1 * f2y2 )
    }
    part2 <- part2 + f1y1^2 + f2y1^2
  }
  safari_rho_Spearman <- 6 * part1 + 3 * part2 - 3
  return(safari_rho_Spearman)
}
