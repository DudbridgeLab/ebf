# Calculate Poisson EBFs with shrinkage to non-parametric posterior

ebf.poisson.npml <- function(x, interval, index, xmin, xmax, complement=FALSE,
                             points, nsupport, tol, nboot) {
  npoints = length(points)

  # total intervals for mixture component distributions
  if (nboot == 0) np.interval = sum(interval[points])
  else np.interval = interval[points]

  # posterior marginal likelihood
  pml = rep(0, length(x))
  pml.numer = rep(0, length(x))
  pml.denom = rep(0, length(x))

  # bootstrap loop
  for(boot in 0:nboot) {
    if (boot == 0) x.boot = x[points]
    else {
      if (npoints == 1) x.boot = rpois(npoints, np.support.mle) * interval[points]
      else x.boot = rpois(npoints,
                          sample(np.support.mle, npoints, replace=T, prob=np.weights.mle) *
                            interval[points])
    }

    # initialise
    np.support = seq(min(x.boot/interval[points], na.rm=T),
                     max(x.boot/interval[points], na.rm=T),
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
          density.matrix[,j] = dnbinom(x.boot, np.support[j]*np.interval + 1,
                                       np.interval/(np.interval+interval[points]))
        }
        # bootstrap
        else {
          density.matrix[,j] = dpois(x.boot, np.support[j]*interval[points])
        }
      }
      # prob of each data point marginal over support points
      likelihood = density.matrix %*% np.weights

      # posterior supports from MLE
      np.support = (t(x.boot / likelihood) %*% density.matrix) /
        (t(interval[points] / likelihood) %*% density.matrix)

      # posterior weights from Bayes formula
      new.weights = (t(1/likelihood) %*% density.matrix) * np.weights
      # normalised
      np.weights = new.weights / sum(new.weights)

      np.support = as.vector(np.support)
      np.weights = as.vector(np.weights)
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
        if (complement == FALSE)
          area1 = dnbinom(x[i], np.support*np.interval + 1,
                          np.interval / (np.interval+interval[i])) *
            (pgamma(xmax, x[i] + np.support*np.interval + 1,
                    interval[i] + np.interval) -
               pgamma(xmin, x[i] + np.support*np.interval + 1,
                      interval[i]+np.interval))
        else
          area1 = dnbinom(x[i], np.support*np.interval + 1,
                          np.interval / (np.interval+interval[i])) *
            pgamma(xmin, x[i] + np.support*np.interval + 1,
                   interval[i] + np.interval) +
            dnbinom(x[i], np.support*np.interval + 1,
                    np.interval / (np.interval+interval[i])) *
            pgamma(xmax, x[i] + np.support*np.interval + 1,
                   interval[i] + np.interval, lower=F)

        # normalising terms
        if (complement == FALSE)
          area2 = pgamma(xmax, np.support[j]*np.interval + 1, np.interval) -
            pgamma(xmin, np.support[j]*np.interval + 1, np.interval)
        else
          area2 = mean(pgamma(xmin, np.support[j]*np.interval + 1, np.interval) +
                         pgamma(xmax, np.support[j]*np.interval + 1, np.interval, lower=F))
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
            sum(np.weights[area2] * dpois(x[i], np.support[area2]*interval[i]))
          pml.denom[i] = pml.denom[i] + sum(np.weights[area2])
        }
      }
    }

  } # bootstrap loop

  pml.numer = pml.numer / (nboot+1)
  pml.denom =  pml.denom / (nboot+1)

  if (xmax==Inf) bias = 0.5
  else bias = 0

  for(i in index) {
    pml[i] = pml.numer[i] / pml.denom[i] /
      exp(bias * 2 / (npoints+1))
  }

  pml
}
