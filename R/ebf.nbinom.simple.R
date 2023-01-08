# Calculate simple negative binomial EBFs

ebf.nbinom.simple <- function(x, size, xmin, xmax, shape, complement=FALSE) {

  # posterior mean likelihood
  area1 = rep(0, length(x))
  bias = rep(0, length(x))
  for(i in 1:length(x)) {
    if (complement == FALSE) {
      area1[i] = pbeta(xmax, 2*size[i]-shape, 2*x[i]+shape) -
        pbeta(xmin, 2*size[i]-shape, 2*x[i]+shape)
    } else {
      area1[i] = pbeta(xmin, 2*size[i]-shape, 2*x[i]+shape) +
        pbeta(xmax, 2*size[i]-shape, 2*x[i]+shape, lower=F)
    }

    area1[i] = exp(log(area1[i]) + lchoose(x[i]+size[i]-1, size[i]-1) +
                     lbeta(2*size[i]-shape, 2*x[i]+shape))# -
                     #lbeta(size[i]-shape, x[i]+shape))

    # bias
    if (shape == 1 & size[i] <= length(nbinom.bias) &
        (complement == FALSE & xmin == 0 & xmax == 1 |
         complement == TRUE & xmin == xmax)) {
      nbinom.bias[size[i]]
    }
    else {
      bias[i] = compute.nbinom.bias(size[i], xmin, xmax, shape, complement)
    }
    bias[i] = compute.nbinom.bias(size[i], xmin, xmax, shape, complement)
  }

  # normalising term
  if (complement == FALSE ) {
    area2 = pbeta(xmax, size[i]-shape, x+shape) -
      pbeta(xmin, size[i]-shape, x+shape)
  } else {
    area2 = pbeta(xmin, size[i]-shape, x+shape) +
      pbeta(xmax, size[i]-shape, x+shape, lower=F)
  }
area2=1

  # EBF
  area1 / area2 / exp(bias)
}

