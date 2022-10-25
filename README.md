# ebf
### Empirical Bayes factors from summary statistics

An empirical Bayes factor (EBF) is one in which the distribution of parameters is estimated from the data in some way.
Here we use posterior distributions derived from vague priors, with adjustment for using the data
first to estimate the posterior and then again to calculate the Bayes factor.  The concept is similar to that of Watanabe's
widely applicable information criterion (WAIC), and in some cases the EBF and WAIC are identical.

Informally an EBF is the Bayes factor we would expect in an independent replication of the experiment, using the posterior
from the present data as the new prior.  But in classical testing scenarios the EBF can also have a fixed sampling distribution under
the null hypothesis.  It is therefore interpretable under frequentist, Bayesian and information perspectives,
while standing itself as a measure of evidence.

Furthermore the following argument gives a calibration of EBF to p-values that reflects intuitions.  Define weaker belief to be
a probability that is easily changed by a Bayes factor, and stronger belief as probability that is less easily changed.
The sharpest distinction between weaker and stronger belief is at probability 0.789 (Dudbridge 2022), and with that boundary 
a Bayes factor of 3.73 would update any weaker belief to a stronger belief.  For several common tests, an EBF of 3.73 corresponds to a p-value
of around 0.05.  EBFs of roughly 4, 15, 50 and 200 correspond to increasingly strong levels of evidence.

In simultaneous testing situations, the distribution can be estimated from the ensemble of tests, giving an automatic "correction"
for multiple testing.  This is now similar to empirical Bayes methods, notably by Efron, used for obtaining
posterior probabilities of hypotheses.  Here however we focus on the evidence in the data and formally separate it from the inferential step.

The EBF is a promising framework for the objective measurement and interpretation of statistical evidence.
This package implements EBF calculations where the "data" are statistics from classical hypothesis tests.
At the point at which we would calculate a p-value from a test statistic, this package provides a corresponding EBF.

To install within R, `devtools::install_github("DudbridgeLab/ebf").`

For main documentation, `help(ebf::ebf)`.

## Quick start
Test whether a parameter `x` is zero, assuming its standard error is known

`ebf.norm(x, se)`

One-sided test

`ebf.norm(x, se, h1=c(0, Inf))`

Compare x>0 to x<0

`ebf.norm(x, se, h0=c(-Inf,0), h1=c(0, Inf))`

T-test with estimated standard error

`ebf.t(x, se, df)`

Binomial test of `k` successes in `n` trials, null hypothesis of probability 0.5

`ebf.binom(k, n)`

Null hypothesis of probability 0.25

`ebf.binom(k, n, h0=0.25)`

Chi-square test

`ebf.chisq(chisq, df)`

## Citation

Dudbridge F (2022) Units of evidence and expected Bayes factors for objective reporting of statistical evidence (submitted).
