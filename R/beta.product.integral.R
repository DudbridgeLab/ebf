# Integral of product of beta likelihood and gamma posterior
# Just up to the 99.9% quantile of the posterior

beta.product.integral <- function(x, y, npoints, f) {
  if (f) integrate(function(b) exp(dbeta(x, 1, b, log=T) +
                                     dgamma(b, npoints+1, npoints/y, log=T)),
                   #1, qgamma(0.999, npoints+1, npoints/y))$value
                   1, Inf)$value
  else exp((npoints+1)*log(npoints/y) - log(1-x)- (npoints + 2)*log(npoints/y-log(1-x)) +
             log(npoints+1))
}

