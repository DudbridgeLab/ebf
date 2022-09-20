# Calculate simple t EBFs

ebf.t.simple <- function(x, se, xmin, xmax, df, complement=FALSE) {
  # expand into a vector
  se = rep(0, length(x)) + se

  # posterior mean likelihood
  if (complement == FALSE) {
    area1 = sapply(1:length(x), function(i) {
      t.product.integral(x[i], se[i], x[i], se[i], df, xmin, xmax)
    })
  } else {
    area1 = sapply(1:length(x), function(i) {
      t.product.integral(x[i], se[i], x[i], se[i], df, -Inf, xmin) +
        t.product.integral(x[i], se[i], x[i], se[i], df, xmax, Inf)
    })
  }

  # normalising term
  if (complement == FALSE ) {
    area2 = pt((xmax-x)/se, df) - pt((xmin-x)/se, df)
  } else {
    area2 = pt((xmin-x)/se, df) + pt((xmax-x)/se, df, lower=F)
  }

  # bias
  if (df<=nrow(ebf::t.bias)) bias = ebf::t.bias$bias[df] else bias = 0.5

  # EBF
  area1 / area2 / exp(bias * area2)
}

