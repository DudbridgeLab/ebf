% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebf.binom.R
\name{ebf.binom}
\alias{ebf.binom}
\title{Empirical Bayes factors for binomial tests}
\usage{
ebf.binom(
  x,
  size,
  index = NULL,
  h0 = 0.5,
  h1 = NULL,
  shrink = FALSE,
  shape = 1,
  points = NULL,
  pi0 = 0
)
}
\arguments{
\item{x}{Vector of statistics.}

\item{size}{Vector containing the numbers of trials in each test.
The numbers of successes are contained in \code{x}.}

\item{index}{Vector of indices in \code{x}.  EBFs are calculated for
those elements only. If \code{NULL}, EBFs are calculated for all elements.}

\item{h0}{If a scalar, the value of a point null hypothesis.
If a vector with two elements, the lower and upper
bounds of the null hypothesis.}

\item{h1}{If a scalar, the value of a point alternative hypothesis.
If a vector with two elements, the lower and upper bounds of
the alternative hypothesis.  If \code{NULL} (default), the alternative
hypothesis is the complement of \code{h0}.}

\item{shrink}{If \code{TRUE}, shrinks EBFs to a posterior
distribution estimated from all the elements of \code{x}.}

\item{shape}{Parameter of the symmetric prior Beta distribution.}

\item{points}{Vector of indices in \code{x} to use in estimating shrinkage EBFs.}

\item{pi0}{Proportion of true null hypotheses in shrinkage EBFs.}
}
\value{
A data frame containing the following components:
  \itemize{
  \item{\code{index} Indices of elements in \code{x} for which EBFs are calculated.}
    \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
    \item{\code{ebf.units} Expressed as units of evidence, that is
     on the log scale with base 3.73.  One unit updates
     weaker belief to stronger belief, or weaker disbelief to weaker belief,
     according to a rational definition of weaker and stronger belief.
     See references for further details.}
    \item{\code{p} P-values if \code{h0} is a point null hypothesis.}
    \item{\code{p.log10} Expressed in -log10 units.}
    \item{\code{ebf.shrink} Shrunken EBFs if specified by \code{shrink==TRUE}.}
    \item{\code{ebf.shrink.units} Expressed as units of evidence.}
  }
}
\description{
Calculates empirical Bayes factors (EBFs) for one-sample binomial tests of proportions.
}
\details{
EBFs are calculated in favour of \code{h1}.

The default test has \code{h0=0.5}, with two-sided alternative.
For one-sided alternatives, use \code{h1=c(0.5,1)} or \code{h1=c(0,0.5)}.
To test higher values against lower values, use \code{h0=c(0,0.5)},
in which case \code{h1} defaults to \code{c(0.5,1)}.
To test lower values against higher values, use \code{h0=c(0.5,1)}.

The EBF includes bias adjustments to the log posterior marginal likelihoods.
When the hypothesis is \code{[0,1]}, pre-computed adjustments are used for \code{size} from 1 to 100.
For higher values of \code{size}, the asymptotic adjustment of 0.5 is used.

If a normal approximation is acceptable, use \code{\link{ebf.norm}}.

For shrinkage EBFs, each test is evaluated using the posterior distribution
obtained from every other test.
The calculation may be time-consuming when the number of tests is large.
To reduce the computation a subset of tests can be used instead,
with indices specified in \code{points}.  Each individual test is appended
to this list as it is evaluated.

Further reduction can be achieved by calculating EBFs only for the elements
of \code{x} given by \code{index}.  The posterior distribution is still calculated using all
the points in \code{x}, or those specified in \code{points}.

For test \code{i}, the contributions from tests \code{j!=i} can be weighted
by \code{pi0} to model the probability that the hypothesis is true for those
other tests.  The default \code{pi0=0} assumes that the hypothesis is true
for all tests.  If \code{pi0=1}, the single-test EBF is recovered.
}
\references{
Dudbridge F (2023) Empirical Bayes factors for common hypothesis tests.
arXiv:2301.11057
}
\author{
Frank Dudbridge
}
