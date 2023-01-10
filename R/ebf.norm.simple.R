# Calculate simple normal EBFs

ebf.norm.simple <- function(x, s, xmin, xmax, complement=FALSE) {

  if (complement == FALSE) {
    area1 = pnorm(xmax, x, s/sqrt(2)) - pnorm(xmin, x, s/sqrt(2))
    area2 = pnorm(xmax, x, s) - pnorm(xmin, x, s)
    tails = is.infinite(xmin) + is.infinite(xmax)
  } else {
    area1 = pnorm(xmin, x, s/sqrt(2)) + pnorm(xmax, x, s/sqrt(2), lower=F)
    area2 = pnorm(xmin, x, s) + pnorm(xmax, x, s, lower=F)
    tails = is.finite(xmin) + is.finite(xmax)
  }

  # EBF
  area1 / area2 / (2*s*sqrt(pi)) / exp(0.5*tails/2)
}
