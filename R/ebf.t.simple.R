# Calculate simple t EBFs

ebf.t.simple <- function(x, se, xmin, xmax, df, complement=FALSE) {
  df[is.infinite(df)] = 1e8

  # posterior mean likelihood
  area1 = rep(0, length(x))
  if (complement == FALSE) {
    area1 = exp(lgamma((df+1)/2)*2 + lgamma(df+0.5) -
                  lgamma(df/2)*2 - lgamma(df+1) - (log(df) + log(pi))/2) *
      (pt((xmax-x)/se * sqrt(1/df + 2), 2*df+1) -
         pt((xmin-x)/se * sqrt(1/df + 2), 2*df+1)) /se
  } else {
    area1 = exp(lgamma((df+1)/2)*2 + lgamma(df+0.5) -
                  lgamma(df/2)*2 - lgamma(df+1) - (log(df) + log(pi))/2) *
      (pt((xmin-x)/se * sqrt(1/df + 2), 2*df+1) +
         pt((xmax-x)/se * sqrt(1/df + 2), 2*df+1, lower=F)) / se
  }


  # bias
  bias = rep(0, length(x))
  for(i in 1:length(x)) {
    if (df[i] <= length(t.bias))
      bias[i] = t.bias[df[i]]
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

