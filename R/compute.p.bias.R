#' Compute EBF bias for p-values
#'
#' Computes the bias in the log posterior marginal likelihood for p-values.
#'
#' The p-values are assumed to follow a Beta(1, \code{shape}) distribution,
#' where \code{shape}>=1.
#' The bias turns out to be about log(5/2) = 0.92 for most values of \code{shape}.
#'
#' @param shape Second parameter for Beta distribution of p-values.
#'
#' @param nsample Number of random samples taken to estimate bias.
#' If 0 (default), analytic calculation is made.
#'
#' @import stats
#'
#' @return The bias for the given shape parameter.
#' This is about log(5/2) = 0.92 for the bulk of the distribution.
#'
#' @template references
#'
#' @export

compute.p.bias <- function(shape=1, nsample=0) {

  bias = rep(0, length(shape))
  for(i in 1:length(shape)) {

    # get it exactly if numerical precision allows
    if (nsample == 0) {
      bias[i] = cubature::cubintegrate(function(x) {
        q.obs = log(1-x[1])
        q.rep = log(1-x[1]) + log(1-x[2])
        pbf.obs = -1/4/q.obs^3 + 1/2/q.obs^2 - 1/2/q.obs
        pbf.rep = -2/q.rep^3 + 2/q.rep^2 - 1/q.rep
        dbeta(x[1], 1, shape[i]) * dbeta(x[2], 1, shape[i]) *
          (log(pbf.obs) - log(pbf.rep))
      }, c(0,0), c(1,1))$integral
    }

    # for large values, get it by Monte Carlo
    else {
      bias.i = rep(0, nsample)
      p.obs = rbeta(nsample, 1, shape[i])
      p.rep = rbeta(nsample, 1, shape[i])
      q.obs = log(1-p.obs)
      q.rep = log(1-p.obs) + log(1-p.rep)
      pbf.obs = -1/4/q.obs^3 + 1/2/q.obs^2 - 1/2/q.obs
      pbf.rep = -2/q.rep^3 + 2/q.rep^2 - 1/q.rep
      w = which(pbf.obs>0 & pbf.rep>0)
      bias.i[w] =
        (log(pbf.obs[w]) - log(pbf.rep[w]))
      bias[i] = mean(bias.i[w])
    }
  }
  bias
}

