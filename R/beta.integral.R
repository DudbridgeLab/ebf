# Integral of beta(1, b) density raised to power of npoints

beta.integral <- function(x, npoints) {
  pgamma(1, npoints+1, 1/x*npoints, lower=F)
  #ntegrate(function(b) b^npoints * (1-x)^(npoints*(b-1)), 1, Inf)$value
  #integrate(function(b) dbeta(x, npoints, npoints*b), 1, Inf)$value
}

