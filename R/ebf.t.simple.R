# Calculate simple t EBFs

ebf.t.simple <- function(x, se, xmin, xmax, df, complement=FALSE) {
  # posterior mean likelihood
  area1 = rep(0, length(x))
  bias = rep(0, length(x))
  for(i in 1:length(x)) {
    if (complement == FALSE) {
      area1[i] = t.product.integral(x[i], se[i], df[i], x[i], se[i], df[i], xmin, xmax)
    } else {
      area1[i] = t.product.integral(x[i], se[i], df[i], x[i], se[i], df[i], -Inf, xmin) +
        t.product.integral(x[i], se[i], df[i], x[i], se[i], df[i], xmax, Inf)
    }
    # bias
    if (df[i] <= length(t.bias))
      bias[i] = ebf::t.bias[df[i]]
    else
      bias[i] = 0.5

  }

  # normalising term
  if (complement == FALSE ) {
    area2 = pt((xmax-x)/se, df) - pt((xmin-x)/se, df)
  } else {
    area2 = pt((xmin-x)/se, df) + pt((xmax-x)/se, df, lower=F)
  }

    # EBF
  area1 / area2 / exp(bias * area2)
}

