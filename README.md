# ebf
### Empirical Bayes factors from summary statistics

An empirical Bayes factor (EBF) is one in which the distribution of parameters is estimated from the data in some way.
Here we use posterior distributions derived from vague priors, with adjustment for using the data twice:
first to estimate the posterior and then again to calculate the Bayes factor.  The concept is similar to Watanabe's
widely applicable information criterion (WAIC), and in some cases the EBF and WAIC are identical.

Informally an EBF is the Bayes factor we would expect in an independent replication of the experiment, using the posterior
from the present data as the new prior.  Inferences from EBFs are therefore expected to replicate.
But in classical testing scenarios the EBF can also have a fixed sampling distribution under
the null hypothesis.  Indeed an approximate general conversion from p-values is EBF=10p.
The EBF is therefore interpretable under frequentist, Bayesian and information perspectives,
while standing itself as a measure of evidence.

The following argument gives a calibration of EBF to p-values that reflects intuitions (Dudbridge 2023).  Define weaker belief to be
a probability that is easily changed by a Bayes factor, and stronger belief as probability that is less easily changed.
The sharpest distinction between weaker and stronger belief is where the third derivative of the logistic function is zero.
This corresponds to probability 0.789, and with that boundary a Bayes factor of 3.73 would update any weaker belief to a stronger belief.
For several common tests, an EBF of 3.73 corresponds to a p-value of around 0.05.
EBFs of roughly 4, 15, 50 and 200 correspond to increasingly strong levels of evidence.

In multiple testing situations, the distribution can be estimated from the ensemble of tests, giving a version of
Storey's optimal discovery procedure (ODP) for simultaneous inference.

The EBF is a promising framework for the objective measurement and interpretation of statistical evidence.
This package implements EBF calculations where the "data" are statistics from classical hypothesis tests.
At the point at which we would calculate a p-value from a test statistic, this package provides a corresponding EBF.

To install within R, `devtools::install_github("DudbridgeLab/ebf").`

For main documentation, `help(ebf::ebf)`.

## Quick start
Test whether a parameter `x` is zero, assuming its standard error `s` is known

`ebf.norm(x, s)`

One-sided test

`ebf.norm(x, s, h1=c(0, Inf))`

Compare x>0 to x<0

`ebf.norm(x, s, h0=c(-Inf,0), h1=c(0, Inf))`

T-test with estimated standard error

`ebf.t(x, s, df)`

Binomial test of `k` successes in `n` trials, null hypothesis of probability 0.5

`ebf.binom(k, n)`

Null hypothesis of probability 0.25

`ebf.binom(k, n, h0=0.25)`

Chi-square test

`ebf.chisq(chisq, df)`

## Citations

Dudbridge F (2023) Empirical Bayes factors for common hypothesis tests. arXiv:2301.11057

