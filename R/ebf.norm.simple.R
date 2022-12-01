# Calculate simple normal EBFs

ebf.norm.simple <- function(x, se, xmin, xmax, complement=FALSE) {

  if (complement == FALSE) {
    area1 = pnorm(xmax, x, se/sqrt(2)) - pnorm(xmin, x, se/sqrt(2))
    area2 = pnorm(xmax, x, se) - pnorm(xmin, x, se)
    tails = is.infinite(xmin) + is.infinite(xmax)
  } else {
    area1 = pnorm(xmin, x, se/sqrt(2)) + pnorm(xmax, x, se/sqrt(2), lower=F)
    area2 = pnorm(xmin, x, se) + pnorm(xmax, x, se, lower=F)
    tails = is.finite(xmin) + is.finite(xmax)
  }

  # EBF
  area1 / area2 / (2*se*sqrt(pi)) / exp(0.5*tails/2)
}
