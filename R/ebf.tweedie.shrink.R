# Calculate shrinkage EBFs with Tweedie formula
# Does not adjust for interval hypotheses,
# which is done in ebf.t.shrink and ebf.t.shrink
#
# Requires the splines library

ebf.tweedie.shrink <- function(x, se, nbin=100) {
  # transform to z scores (or t statistics)
  z = x/se

  # histogram of z scores
  h = hist(z, breaks=nbin)
  bincounts = h$counts
  midpoints = h$mids

  # fit splines by Poisson regression
  fit = glm(bincounts ~ splines::ns(midpoints, df=7), family=poisson)

  # estimate of posterior mean likelihoods
  odp.tweedie = exp(predict.glm(fit, newdata=data.frame(midpoints=z))) /
    diff(midpoints)[1] / length(x)

  # transform back
  odp.tweedie / se
}
