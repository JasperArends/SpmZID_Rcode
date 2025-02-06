######################################################################
# COPULAS
######################################################################

# Upper Frechet-Hoeffding copula
copula_upper <- function(u, v, cop_par) {
  min(u, v)
}

# Lower Frechet-Hoeffding copula
copula_lower <- function(u, v, cop_par) {
  max(u + v - 1, 0)
}

# Mixed Frechet copula
copula_frechet <- function(u, v, cop_par) {
  ifelse(rep(cop_par > 0, length(u)),
         (1 - cop_par) * u * v + cop_par * pmin(u, v),
         (1 + cop_par) * u * v - cop_par * pmax(u + v - 1, 0)
  )
}

# Frank copula
copula_frank <- function(u, v, theta) {
  -(1/theta) * log( 1 + (exp(-theta * u) - 1) * (exp(-theta * v) - 1)/(exp(-theta) - 1) )
}
