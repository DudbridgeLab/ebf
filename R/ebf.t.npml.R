# Calculate t EBFs with shrinkage to non-parametric posterior

ebf.t.npml <- function(x, se, df, index, xmin, xmax, complement=FALSE,
                       points, nsupport, tol, nboot) {
  npoints = length(points)

  # variance of mixture components
  np.sd = 1 / sqrt(sum(1 / se[points]^2))
  np.df = sum(df[points]) + npoints - 1


  # posterior marginal likelihoods
  pml = rep(0, length(x))
  pml.numer = rep(0, length(x))
  pml.denom = rep(0, length(x))

  # bootstrap loop
  for(boot in 0:nboot) {
    if (boot == 0) x.boot = x[points]
    else {
      if (npoints == 1) x.boot = x[points]
      else x.boot = sample(np.support.mle, npoints, replace=T, prob=np.weights.mle)
      x.boot = x.boot + rt(npoints, df[points]) * se[points]
    }

    # initialise
    np.support = seq(min(x.boot, na.rm=T), max(x.boot, na.rm=T),
                     length=nsupport)
    np.weights = rep(1/nsupport, nsupport)
    llhd = 0
    newllhd = 0
    density.matrix = matrix(0, nrow=npoints, ncol=nsupport)

    # EM algorithm to maximise likelihood
    while(llhd==0 | abs((newllhd-llhd)/llhd) > tol) {
      llhd = newllhd
      # prob of each data point for each support point
      # data in rows, supports in columns
      for(j in 1:nsupport) {
        # analytic
        if (nboot == 0) {
          for(i in 1:npoints)
            density.matrix[i,j] =
              t.product.integral(x.boot[i], se[points[i]], df[points[i]],
                                 np.support[j], np.sd, np.df, -Inf, Inf)
        }
        # bootstrap
        else {
          density.matrix[,j] = dt((x.boot - np.support[j]) / se[points], df[points]) /
            se[points]
        }
      }

      # prob of each data point marginal over support points
      likelihood = density.matrix %*% np.weights

      # posterior supports
      density.vector = rep(0,npoints)
      for(j in 1:nsupport) {
        # expected data count
        prob = density.matrix[,j]*np.weights[j]/likelihood
        # maximising expected data likelihood
        # analytic
        if (nboot == 0){
          np.support[j] = optim(np.support[j], function(y) {
            for(i in 1:npoints)
              density.vector[i] = t.product.integral(x.boot[i], se[points[i]], df[points[i]],
                                                     y, np.sd, np.df, -Inf, Inf)
            -sum(log(density.vector) * prob)
          }, method="BFGS")$par
        }
        # bootstrap
        else {
          np.support[j] = optim(np.support[j], function(y) {
            density.vector = dt((x.boot - y) / se[points], df[points], log=T) - log(se[points])
            -sum(density.vector*prob)
          }, method="BFGS")$par
        }
      }
      np.support =  as.vector(np.support)

      # posterior weights
      new.weights = (t(1/likelihood) %*% density.matrix) * np.weights
      # normalised
      np.weights = new.weights / sum(new.weights)
      np.weights =  as.vector(np.weights)

      newllhd = sum(log(likelihood))
    }
    if (boot == 0) {
      np.support.mle = np.support
      np.weights.mle = np.weights
    }

    # shrinkage EBFs
    # analytic
    if (nboot == 0) {
      for(i in index) {
        area1 = rep(0, nsupport)
        for(j in 1:nsupport) {
          if (complement == FALSE) {
            area1[j] = t.product.integral(x[i], se[i], df[i],
                                          np.support[j], np.sd, np.df,
                                          xmin, xmax)
          }
          else {
            area1[j] = t.product.integral(x[i], se[i], df[i],
                                          np.support[j], np.sd, np.df,
                                          -Inf, xmin) +
              t.product.integral(x[i], se[i], df[i],
                                 np.support[j], np.sd, np.df,
                                 xmax, Inf)
          }
        }

        # normalising terms
        if (complement == FALSE)
          area2 = pt((xmax - np.support) / np.sd, np.df) -
            pt((xmin - np.support) / np.sd, np.df)
        else
          area2 = pt((xmin - np.support) / np.sd, np.df) +
            pt((xmax - np.support) / np.sd, np.df, lower=F)

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
            sum(np.weights[area2] * dt((x[i] - np.support[area2]) / se[i], df[i]) / se[i])
          pml.denom[i] = pml.denom[i] + sum(np.weights[area2])
        }
      }
    }
  } # bootstrap loop

  pml.numer = pml.numer / (nboot+1)
  pml.denom = pml.denom / (nboot+1)
  for(i in index) {
    # bias correction
    if (df[i] <= length(t.bias))
      bias = t.bias[df[i]]
    else
      bias = 0.5
    pml[i] = pml.numer[i] / pml.denom[i] /
      exp(pml.denom[i] * bias * 2 / (npoints+1))
  }

  pml
}
