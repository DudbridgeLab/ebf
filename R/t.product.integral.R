# Integral of product of t densities

t.product.integral <- function(x1, s1, x2, s2, df, xmin, xmax) {
  integrate(function(x)
    dt((x1-x)/s1, df) / s1 *
      dt((x2-x)/s2, df) / s2,
    xmin, xmax, stop.on.error=FALSE) $value
}
