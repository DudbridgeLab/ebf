# ebf
Empirical Bayes factors from summary statistics.

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
