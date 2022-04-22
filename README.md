
bsreg: Bayesian Spatial Regression Models
=======

[![CRAN](https://www.r-pkg.org/badges/version/bsreg)](https://cran.r-project.org/package=bsreg)
[![month](https://cranlogs.r-pkg.org/badges/bsreg)](https://www.r-pkg.org/pkg/bsreg)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/bsreg)](https://www.r-pkg.org/pkg/bsreg)

Estimation of Bayesian spatial models.

Installation
-------

**bsreg** is available on [CRAN](https://CRAN.R-project.org/package=bsreg). The development version can be installed from GitHub.
``` r
install.packages("bsreg")
devtools::install_github("nk027/bsreg")
```

Demonstration
-------

``` r
# Load the package
library("bsreg")

# Estimate a Bayesian linear model using cigarette demand data
x <- blm(log(sales) ~ log(price), data = cigarettes, lags = 1)

# Check convergence via trace and density plots
plot(x)
```

References
-------

Nikolas Kuschnig (2021). Bayesian spatial econometrics and the need for software. *Working Paper*, DOI: [10.13140/RG.2.2.11269.68328](https://doi.org/10.13140/RG.2.2.11269.68328).
