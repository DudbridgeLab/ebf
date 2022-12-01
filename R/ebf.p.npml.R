# Calculate EBFs for p-values with shrinkage to non-parametric posterior

ebf.p.npml <- function(x, index, points, nsupport, tol, nboot) {
  npoints = length(points)

  # shape of mixture components
  if (nboot == 0) np.shape = npoints
  else np.shape = 1

  # posterior marginal likelihood
  pml = rep(0, length(x))

  # bootstrap loop
  for(boot in 0:nboot) {
    if (boot == 0) x.boot = x[points]
    else {
      if (npoints==1) x.boot = rbeta(1, 1, np.support.mle)
      else x.boot = rbeta(npoints,
                          1, sample(np.support.mle, npoints, replace=T, prob=np.weights.mle))
    }

    # initialise
    np.support = seq(min(x.boot, na.rm=T), max(x.boot, na.rm=T),
                     length=nsupport)
    np.weights = rep(1/nsupport, nsupport)
    llhd = 0
    newllhd = 0
    density.matrix = matrix(0, nrow=npoints, ncol=nsupport)

    ### EM algorithm to maximise likelihood
    start = TRUE
    while(start | abs((newllhd-llhd)) > tol) {
      start = FALSE
      llhd = newllhd

      # prob of each data point for each support point
      # data in rows, supports in columns
      for(j in 1:nsupport) {
        # analytic
        if (FALSE & nboot == 0) {
          for(i in 1:npoints)
            density.matrix[i,j] = beta.product.integral(x.boot[i], np.support[j], npoints)# /
#              beta.integral(np.support[j], npoints)
        }
        # bootstrap
        else {
          density.matrix[,j] = dbeta(x.boot, 1, np.support[j])
        }
      }

      # prob of each data point marginal over support points
      likelihood = density.matrix %*% np.weights

      # posterior supports from MLE
      #analytic
      if (FALSE & nboot == 0) {
        density.vector = rep(0,npoints)
        for(j in 1:nsupport) {
          # expected data count
          prob = density.matrix[,j]*np.weights[j]/likelihood
          # maximising expected data likelihood
          np.support[j] = optimise(function(y) {
            denom = beta.integral(y,npoints)
            for(i in 1:npoints)
              density.vector[i] = beta.product.integral(x.boot[i], y, npoints) / denom
            w = which(!is.na(density.vector))
            sum(log(density.vector[w]) * prob[w])
          }, c(0,1), maximum=T)$maximum
        }
      }
      # bootstrap
      else {
        np.support = (t(-1/sum(log(1 - x.boot)) / likelihood) %*% density.matrix)
        np.support[np.support<1] = 1
      }
      # posterior weights from Bayes formula
      new.weights = (t(1/likelihood) %*% density.matrix) * np.weights
      # normalised
      np.weights = new.weights / sum(new.weights)
      #print(np.support)

      np.support =  as.vector(np.support)
      np.weights =  as.vector(np.weights)
      newllhd = sum(log(likelihood))
    }
    if (boot == 0) {
      np.support.mle = np.support
      np.weights.mle = np.weights
      print(np.support)
    }

    ### shrinkage EBFs
    # analytic
    if (nboot == 0)  {
      for(j in 1:nsupport) {
        denom = beta.integral(np.support[j],npoints)
        for(i in index) {
          ratio = beta.product.integral(x[i], np.support[j], npoints) / denom
          if (is.na(ratio)) ratio = 1
          pml[i] = pml[i] + np.weights[j] * ratio
        }
      }
    }
    # bootstrap
    else {
      for(i in index) {
        pml[i] = pml[i] +
          sum(np.weights * dbeta(x[i], 1, np.support))
      }
    }

  } # bootstrap loop

  pml = pml / (nboot+1)
  for(i in index)
    pml[i] =  pml[i] / exp(2 / (npoints + 1) * log(5/2))

  pml
}
