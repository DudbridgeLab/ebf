# Integral of product of two beta(1, b) densities

beta.product.integral <- function(x, y, npoints) {
  integrate(function(b) exp(dbeta(x, 1, b, log=T) + dgamma(b, npoints+1, 1/y*npoints, log=T)),
            1, Inf)$value
  #integrate(function(b) b^(npoints+1) * (1-x)^(b-1) * (1-y)^(npoints*(b-1)), 1, Inf)$value
  #integrate(function(b) dbeta(x, 1, b) * dbeta(y, npoints, npoints*b), 1, Inf)$value
}

