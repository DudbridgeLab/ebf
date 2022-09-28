# Calculate simple normal EBFs

ebf.norm.simple <- function(x, se, xmin, xmax, complement=FALSE) {

  if (complement == FALSE) {
    # posterior mean likelihood
    area1 = pnorm(xmax, x, se/sqrt(2)) - pnorm(xmin, x, se/sqrt(2))

    # normalising term
    area2 = pnorm(xmax, x, se) - pnorm(xmin, x, se)
  } else {
    area1 = pnorm(xmin, x, se/sqrt(2)) + pnorm(xmax, x, se/sqrt(2), lower=F)
    area2 = pnorm(xmin, x, se) + pnorm(xmax, x, se, lower=F)
  }

  # EBF
  area1 / area2 / (2*se*sqrt(pi)) / exp(0.5*area2)
}
