<!-- README.md is generated from README.Rmd. Please edit that file -->

# Blend

> Robust Bayesian Longitudinal Regularized Semiparametric Mixed Models
<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/Blend)](https://cran.r-project.org/package=Blend)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/Blend)](https://www.r-pkg.org:443/pkg/Blend)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/Blend)](https://www.r-pkg.org:443/pkg/Blend)
<!-- badges: end -->

Our recently developed fully robust Bayesian semiparametric mixed-effect model for high-dimensional longitudinal studies with heterogeneous observations 
can be implemented through this package. This model can distinguish between time-varying interactions and constant-effect-only 
cases to avoid model misspecifications. Facilitated by spike-and-slab priors, this model leads to superior performance in estimation,
identification and statistical inference. In particular, robust Bayesian inferences in terms of valid Bayesian credible intervals on 
both parametric and nonparametric effects can be validated on finite samples. The Markov chain Monte Carlo algorithms of the proposed 
and alternative models are efficiently implemented in 'C++'.

## How to install

  - To install from github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("kunfa/Blend")

  - Released versions of Blend are available on CRAN
    [(link)](https://cran.r-project.org/package=Blend), and can be
    installed within R via

<!-- end list -->

    install.packages("Blend")

## Examples

#### Example.1 (default method: robust structural with sparsity)

    library(Blend)
    data(dat)
    
    fit = Blend(y,x,t,J,kn,degree) 
    
    u = as.numeric(fit$basis$u)
    truth = list(gamma0 = function(u) 2+2*sin(2*pi*u),
     gamma1 = function(u) 2.5*exp(2.5*u-1),
     gamma2 = function(u) 3*u^2-2*u+2,
     gamma3 = function(u) -4*u^3+3,
    gamma4 = function(u) 3-2*u,
     gamma5 = 0.8,
    gamma6 = -1.2,
    gamma7 = 0.7,
     gamma8 = -1.1,
    gamma9 = 0,
    gamma10 = 0,
     gamma11 = 0,
     gamma12 = 0,
     gamma13 = 0,
     gamma14 = 0,
     gamma15 = 0
     )
    out = Coverage(fit, truth)
    print(out)
    
    plot_Blend(fit,sparse=TRUE)
    
#### Example.2 (alternative: robust non-structural with sparsity)

    fit = Blend(y,x,t,J,kn,degree,structural=FALSE) 
    
#### Example.3 (alternative: non-robust structural with sparsity)

    fit = Blend(y,x,t,J,kn,degree, robust=FALSE)
   
#### Example.4 (alternative: non-robust non-structural without sparsity)

    fit = Blend(y,x,t,J,kn,degree, robust=FALSE, structural=FALSE, sparse = FALSE) 
    
    
## Methods

This package provides implementation for methods proposed in

  -Fan, K., Ren, J., Ma, S. and Wu, C. (2026). Robust Bayesian Variable Selection and Inference Under Misspecified Nonparametric Mixed Models in Longitudinal Studies.[Statistics in Medicine. 45(25-27), e70830](https://doi.org/10.1002/sim.70830)
