# Calculate simple normal EBFs

ebf.norm.simple <- function(x, se, xmin, xmax, complement=FALSE) {
  # posterior mean likelihood
  area1 = pnorm(xmax, x, se/sqrt(2)) - pnorm(xmin, x, se/sqrt(2))

  # normalising term
  area2 = pnorm(xmax, x, se) - pnorm(xmin, x, se)
  if (complement == TRUE) {
    area1 = 1-area1
    area2 = 1-area2
  }

  # EBF
  area1 / area2 / (2*se*sqrt(pi)) / exp(0.5*area2)
}
