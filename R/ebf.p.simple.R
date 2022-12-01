# Calculate simple EBFs for p-values

ebf.p.simple <- function(x) {

  q = log(1-x)
  numer = (-1/(2*q^3) + 1/q^2 - 1/q) # factor of 1/2 omitted
  denom = 1/q^2 - 1/q

 # bias is about log(5/2)
  # for small p, we get approximately 1/(10p) !
  numer/denom / 5

}
