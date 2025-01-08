######################################################################
# GENERATE SAMPLE FROM ZERO-INFLATED DISTRIBUTIONS
######################################################################
library(copula)

# Generate sample with Poisson margins with any copula
gen_pois <- function(lambda1, lambda2, pi1, pi2, cop_par, copula, N) {
  A <- c()
  B <- c()
  for (ID in 1:N) {
    x1 <- rbinom(1, 1, pi1) * rpois(1, lambda1)
    A <- c(A,x1)
    F1 <- (1 - pi1) + pi1 * ppois(x1, lambda1) # P(X <= x)
    MPROB1 <- (1 - pi1) * I(x1 == 0) + pi1 * dpois(x1, lambda1) # P(X == x)
    F1M <- ifelse(x1 == 0, 0, (1 - pi1) + pi1 * ppois(x1 - 1, lambda1)) # P(X < x)
    stop <- 0
    tmp <- 0
    uu <- runif(1)
    while (stop != 1) {
      F2 <- (1 - pi2) + pi2 * ppois(tmp, lambda2) # P(Y <= tmp)
      cond_F2 <- (copula(F1, F2, cop_par) - copula(F1M, F2, cop_par)) / MPROB1 # P(V <= v | U = u)
      if(cond_F2 > uu) {
        stop <- 1
        x2 <- tmp
        B <- c(B,x2)
      }
      tmp <- tmp + 1
    }
  }
  
  return(cbind(A, B))
}

# Zero-inflated from the Frechet copula
gen_Frechet <- function(N, alpha) {
  # Copula C(u, v) = (1 - alpha) * Π(u, v) + alpha * M(u, v)
  U <- matrix(nrow=N, ncol=2) # Sample
  U[,1] <- runif(N)
  
  V <- runif(N) # Second entry for independence copula
  Z <- runif(N) # Decision variable
  U[,2] <- (Z <= alpha) * U[,1] + (Z > alpha) * V
  return (U)
}

# Zero-inflated from the Frechet copula with exponential margins
gen_FrecExp <- function(lambda1, lambda2, pi1, pi2, alpha, N) {
  U <- gen_Frechet(N, alpha)
  p1 <- 1 - pi1
  p2 <- 1 - pi2
  X <- matrix(nrow=N, ncol=2)
  X[,1] <- -log(1- pmax(0, U[,1] - p1)/(1 - p1))/lambda1
  X[,2] <- -log(1- pmax(0, U[,2] - p2)/(1 - p2))/lambda2
  return (X)
}
