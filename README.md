
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NAVAECI: Non-Asymptotically Valid and Asymptotically Exact (NAVAE) Confidence Intervals

## Installation

You can install the development version of NAVAECI from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("AlexisDerumigny/NAVAECI")
```

## Confidence interval for the mean

``` r
library(NAVAECI)

n = 4000
x = rexp(n, rate = 1/3)  # so the mean is 3
Navae_ci_mean(x, bound_K = 9, alpha = 0.1)
#> Call: Navae_ci_mean(data = x, alpha = 0.1, bound_K = 9)
#> 
#> CLT-based confidence interval:
#>      5 %     95 % 
#> 2.946852 3.106489 
#> 
#> NAVAE confidence interval:
#>      5 %     95 % 
#> 2.904051 3.149290
```

## Confidence interval for the linear regression

``` r
n = 4000
X1 = rnorm(n, sd = 1)
true_eps = rnorm(n)
Y = 2 + 8 * X1 + true_eps

Navae_ci_ols(Y, X1, K_xi = 3, a = 1.1)
#> Call: Navae_ci_ols(Y = Y, X = X1, a = 1.1, K_xi = 3)
#> 
#> CLT-based confidence intervals:
#>              2.5 %   97.5 % estimate     length
#> intercept 1.962180 2.024766 1.993473 0.06258577
#> X1        7.987254 8.050967 8.019111 0.06371274
#> 
#> NAVAE confidence intervals:
#>              2.5 %   97.5 % regime estimate    length
#> intercept 1.717172 2.269774    Edg 1.993473 0.5526014
#> X1        7.742453 8.295769    Edg 8.019111 0.5533160
```
