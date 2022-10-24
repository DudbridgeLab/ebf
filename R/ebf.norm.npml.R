# Calculate normal EBFs with shrinkage to non-parametric posterior

ebf.norm.npml <- function(x, se, index, xmin, xmax, complement=FALSE,
                          points, nsupport, tol, nboot) {
  npoints = length(points)

  # variance of mixture components
  if (nboot == 0) np.var = 1 / (sum(1 / se[points]^2))
  else np.var = 0

  # posterior marginal likelihood
  pml = rep(0, length(x))
  pml.numer = rep(0, length(x))
  pml.denom = rep(0, length(x))

  # bootstrap loop
  for(boot in 0:nboot) {
    if (boot == 0) x.boot = x[points]
    else {
      if (npoints == 1) x.boot = x[points]
      else x.boot = sample(np.support.mle, npoints, replace=T, prob=np.weights.mle)
      x.boot = x.boot + rnorm(npoints, sd=se[points])
    }

    # initialise
    np.support = seq(min(x.boot, na.rm=T), max(x.boot, na.rm=T),
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
        density.matrix[,j] = dnorm(x.boot, np.support[j], sqrt(se[points]^2 + np.var))
      }

      # prob of each data point marginal over support points
      likelihood = density.matrix %*% np.weights

      # posterior supports from MLE
      np.support = (t(x.boot / (se[points]^2 + np.var) / likelihood) %*% density.matrix) /
        (t(1 / (se[points]^2 + np.var) / likelihood) %*% density.matrix)

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
    if (nboot == 0)  {
      for(i in index) {
        # Laplace approximation
        ivw.var = 1/se[i]^2 + 1/np.var
        ivw.mean = (x[i]/se[i]^2 + np.support/np.var) / ivw.var

        laplace = dnorm(x[i], ivw.mean, se[i]) *
          dnorm(np.support, ivw.mean, sqrt(np.var)) *
          sqrt(2 * pi / ivw.var)

        if (complement == FALSE) {
          area1 = pnorm(xmax, ivw.mean, 1/sqrt(ivw.var)) -
            pnorm(xmin, ivw.mean, 1/sqrt(ivw.var))
          area2 = pnorm(xmax, np.support, sqrt(np.var)) -
            pnorm(xmin, np.support, sqrt(np.var))
        }
        else {
          area1 = pnorm(xmin, ivw.mean, 1/sqrt(ivw.var)) +
            pnorm(xmax, ivw.mean, 1/sqrt(ivw.var), lower=F)
          area2 = pnorm(xmin, np.support, sqrt(np.var)) +
            pnorm(xmax, np.support, sqrt(np.var), lower=F)
        }

        pml.numer[i] = pml.numer[i] + sum(np.weights * laplace * area1)
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
            sum(np.weights[area2] * dnorm(x[i], np.support[area2], se[i]))
          pml.denom[i] = pml.denom[i] + sum(np.weights[area2])
        }
      }
    }


  } # bootstrap loop

  pml.numer = pml.numer / (nboot+1)
  pml.denom = pml.denom / (nboot+1)
  for(i in index)
    pml[i] =  pml.numer[i] / pml.denom[i] /
    exp(pml.denom[i] / (npoints + 1))

  pml
}
