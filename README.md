<!-- README.md is generated from README.Rmd. Please edit that file -->

# Blend

> Robust Bayesian Longitudinal Regularized Semiparametric Mixed Models


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

#### Example.1 (default method)

    library(Blend)
    data(dat)
    
    fit = Blend(y,x,t,J,kn,degree) 
    fit$coefficient 
    Coverage(fit)
    plot_Blend(fit,sparse=TRUE)
#### Example.2 (alternative: robust non-structural)

    fit = Blend(y,x,t,J,kn,degree,structural=FALSE) 
    
#### Example.3 (alternative: non-robust structural)

    fit = Blend(y,x,t,J,kn,degree, robust=FALSE)
   
#### Example.4 (alternative: non-robust non-structural)

    fit = Blend(y,x,t,J,kn,degree, robust=FALSE, structural=FALSE)   
    
## Methods

This package provides implementation for methods proposed in

  -Fan, K., Ren, J., Ma, Shuangge and Wu, C. (2025). robust Bayesian Regularized Semiparametric Mixed Models in Longitudinal Studies. (submitted).
