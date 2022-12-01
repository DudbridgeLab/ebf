# Calculate binomial EBFs with shrinkage to non-parametric posterior

ebf.nbinom.npml <- function(x, size, index, xmin, xmax, shape, complement=FALSE,
                            points, nsupport, tol, nboot) {
  npoints = length(points)

  # total size for mixture component distributions
  if (nboot == 0) np.size = sum(size[points])
  else np.size = 1

  # posterior marginal likelihood
  pml = rep(0, length(x))
  pml.numer = rep(0, length(x))
  pml.denom = rep(0, length(x))

  # bootstrap loop
  for(boot in 0:nboot) {
    if (boot == 0) x.boot = x[points]
    else {
      if (npoints == 1) x.boot = rnbinom(npoints, size[points], np.support.mle)
      else x.boot = rnbinom(npoints, size[points],
                            sample(np.support.mle, npoints, replace=T, prob=np.weights.mle))
    }

    # initialise
    np.support = seq(min(size[points] / (x.boot + size[points]), na.rm=T),
                     max(size[points] / (x.boot + size[points]), na.rm=T),
                     length=nsupport)
    np.weights = rep(1/nsupport, nsupport)
    llhd = 0
    newllhd = 0
    density.matrix = matrix(0, nrow=npoints, ncol=nsupport)

    ### EM algorithm to maximise likelihood
    while(llhd==0 | abs((newllhd-llhd)/llhd) > tol) {
      llhd = newllhd

      # prob of each data point for each support point
      # data in rows, supports in columns
      for(j in 1:nsupport) {
        # analytic
        if (nboot == 0) {
          density.matrix[,j] = exp(lchoose(x.boot + size[points] - 1, size[points] - 1) +
                                     lbeta(size[points] + np.size + shape,
                                           x.boot + (1 - np.support[j]) / np.support[j] * np.size + shape) -
                                     lbeta(np.size + shape,
                                           (1 - np.support[j]) / np.support[j] * np.size + shape))
        }
        #  bootstrap
        else {
          density.matrix[,j] = dnbinom(x.boot, size[points], np.support[j])
        }
      }

      # prob of each data point marginal over support points
      likelihood = density.matrix %*% np.weights

      # posterior supports from MLE
      np.support = (t(size[points] / likelihood) %*% density.matrix) /
        (t((size[points]+x.boot) / likelihood) %*% density.matrix)

      # posterior weights from Bayes formula
      new.weights = (t(1/likelihood) %*% density.matrix) * np.weights
      # normalised
      np.weights = new.weights / sum(new.weights)

      np.support =  as.vector(np.support)
      np.weights =  as.vector(np.weights)
      newllhd = sum(log(likelihood))
    }
    if (boot == 0) {
      np.support.mle = np.support
      np.weights.mle = np.weights
    }

    ### shrinkage EBFs
    # analytic
    if (nboot == 0) {
      for(i in index) {
        area1 = rep(0, nsupport)
        for(j in 1:nsupport) {
          if (complement == FALSE)
            area = pbeta(xmax, size[i] + np.size + shape,
                         x[i] + (1 - np.support[j]) / np.support[j] * np.size + shape) -
              pbeta(xmin, size[i] + np.size + shape,
                    x[i] + (1 - np.support[j]) / np.support[j] * np.size + shape)
          else
            area = pbeta(xmin, size[i] + np.size + shape,
                         x[i] + (1 - np.support[j]) / np.support[j] * np.size + shape) +
              pbeta(xmax, size[i] + np.size + shape,
                    x[i] + (1 - np.support[j]) / np.support[j] * np.size + shape, lower=F)

          area1[j] = exp(log(area) + lchoose(x[i] + size[i] - 1, size[i] - 1) +
                           lbeta(size[i] + np.size + shape,
                                 x[i] + (1 - np.support[j]) / np.support[j] * np.size + shape) -
                           lbeta(np.size + shape,
                                 (1 - np.support[j]) / np.support[j] * np.size + shape))
        }

        # normalising terms
        if (complement == FALSE)
          area2 = pbeta(xmax, np.size + shape,
                        (1 - np.support) / np.support * np.size + shape) -
            pbeta(xmin, np.size + shape,
                  (1 - np.support) / np.support * np.size + shape)
        else
          area2 = pbeta(xmin, np.size + shape,
                        (1 - np.support) / np.support * np.size + shape) +
            pbeta(xmax, np.size +shape,
                  (1 - np.support) / np.support * np.size + shape, lower=F)

        pml.numer[i] = pml.numer[i] + sum(np.weights * area1)
        pml.denom[i] = pml.denom[i] + sum(np.weights * area2)
      }
    }
    # bootstrap
    else {
      if (complement == FALSE)
        area2 = which(np.support >= xmin & np.support <= xmax)
      else
        area2 = which(np.support <= xmin | np.support >= xmax)
      if (length(area2) > 0) {
        for(i in index) {
          pml.numer[i] = pml.numer[i] +
            sum(np.weights[area2] * dnbinom(x[i], size[i], np.support[area2]))
          pml.denom[i] = pml.denom[i] + sum(np.weights[area2])
        }
      }
    }
  } # bootstrap loop

  pml.numer = pml.numer / (nboot+1)
  pml.denom =  pml.denom / (nboot+1)

  for(i in index) {
    # bias
    if (shape == 1 & size[i] <= length(nbinom.bias) &
        (complement == FALSE & xmin == 0 & xmax == 1 |
         complement == TRUE & xmin == xmax)) {
      bias = nbinom.bias[size[i]]
    }
    else {
      bias = compute.nbinom.bias(size[i], xmin, xmax, shape, complement)
    }

    pml[i] = pml.numer[i] / pml.denom[i] /
      exp(bias * 2 / (npoints+1))
  }

  pml
}
