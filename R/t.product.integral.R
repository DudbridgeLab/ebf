# Integral of product of t densities

t.product.integral <- function(x1, s1, df1, x2, s2, df2, xmin=-Inf, xmax=Inf) {
  # what this should be doing
  # integrate(function(x)
  #  dt((x1-x)/s1, df1) / s1 *
  #    dt((x2-x)/s2, df2) / s2,
  #  xmin, xmax, stop.on.error=FALSE) $value

  # faster coding
  if (xmin == xmax) return(0)
  if (df1 == Inf) df1 = 1e8
  if (df2 == Inf) df2 = 1e8
  integrate(function(x)
    (1 + ((x1-x)/s1)^2 / df1) ^ ((-df1-1) / 2) *
      (1 + ((x2-x)/s2)^2 / df2) ^ ((-df2-1) / 2),
    xmin, xmax, stop.on.error=FALSE) $value / s1 / s2 *
   exp(lgamma((df1+1) / 2) + lgamma((df2+1) / 2) -
    log(df1*df2) / 2 - log(pi) - lgamma(df1 / 2) - lgamma(df2 / 2))

}
