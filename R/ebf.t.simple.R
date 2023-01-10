# Calculate simple t EBFs

ebf.t.simple <- function(x, s, xmin, xmax, df, complement=FALSE) {
  df[is.infinite(df)] = 1e8

  # posterior mean likelihood
  area1 = rep(0, length(x))
  if (complement == FALSE) {
    area1 = exp(lgamma((df+1)/2)*2 + lgamma(df+0.5) -
                  lgamma(df/2)*2 - lgamma(df+1) - (log(df) + log(pi))/2) *
      (pt((xmax-x)/s * sqrt(1/df + 2), 2*df+1) -
         pt((xmin-x)/s * sqrt(1/df + 2), 2*df+1)) /s
    tails = is.infinite(xmin) + is.infinite(xmax)
  } else {
    area1 = exp(lgamma((df+1)/2)*2 + lgamma(df+0.5) -
                  lgamma(df/2)*2 - lgamma(df+1) - (log(df) + log(pi))/2) *
      (pt((xmin-x)/s * sqrt(1/df + 2), 2*df+1) +
         pt((xmax-x)/s * sqrt(1/df + 2), 2*df+1, lower=F)) / s
    tails = is.finite(xmin) + is.finite((xmax))
  }


  # bias
  bias = rep(0, length(x))
  for(i in 1:length(x)) {
    if (df[i] <= length(t.bias))
      bias[i] = t.bias[df[i]] * tails/2
    else
      bias[i] = 0.5 * tails/2
  }

  # normalising term
  if (complement == FALSE ) {
    area2 = pt((xmax-x)/s, df) - pt((xmin-x)/s, df)
  } else {
    area2 = pt((xmin-x)/s, df) + pt((xmax-x)/s, df, lower=F)
  }

  # EBF
  area1 / area2 / exp(bias)
}

