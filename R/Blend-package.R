#' @useDynLib Blend, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

"_PACKAGE"
#' @keywords overview
#' @name Blend-package
#' @title Robust Bayesian Longitudinal Regularized Semiparametric Mixed Model
#' @aliases Blend-package
#' @description In this package, we further extend the sparse robust Bayesian mixed models to nonlinear longitudinal interactions. Specifically, the proposed Bayesian semiparametric model is robust not only to outliers and heavy‐tailed distributions of the response variable, but also to the misspecification of interaction effect in the forms other than non-linear interactions. We have developed the Gibbs sampler with the spike‐and‐slab priors to promote sparse identification of appropriate forms of main and interaction effects.
#' In addition to the default method, users can also choose different selection structures for separation of constant and varying effects or not, methods without spike--and--slab priors and non-robust methods. In total, \emph{Blend} provides 8 different methods (4 robust and 4 non-robust) under the random intercept and slope model. All the methods in this package are developed for the first time. Please read the Details below for how to configure the method used.

#' @details The user friendly, integrated interface \strong{Blend()} allows users to flexibly choose the fitting methods by specifying the following parameter:
#' \tabular{rl}{
#' robust: \tab whether to use robust methods for modelling.\cr\cr
#' structural: \tab whether to incorporate structural identification(separation of constant and varying effects) .\cr\cr
#' sparse: \tab whether to use the spike-and-slab priors to impose sparsity.
#' }
#'
#' The function Blend() returns a Blend object that contains the posterior estimates of each coefficients and other useful information for selection().
#' S3 generic functions selection() and print() are implemented for Blend objects.
#' selection() takes a Blend object and returns the variable selection results.
#'
#' @references
#' Fan, K., Ren, J., Ma, Shuangge and Wu, C. (2025+). Robust Bayesian Regularized Semiparametric Mixed Models in Longitudinal Studies. (submitted)
#'
#' Fan, K., Subedi, S., Yang, G., Lu, X., Ren, J., and Wu, C. (2024). Is Seeing Believing? A Practitioner’s Perspective on High-Dimensional Statistical Inference in Cancer Genomics Studies. Entropy, 26(9), 794.
#'
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions.
#' {\emph{Biometrics}, 79(2), 684-694 } \doi{10.1111/biom.13670}
#'
#' Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian regularized quantile varying coefficient model. Computational Statistics & Data Analysis, 187, 107808.
#'
#' Zhou, F., Lu, X., Ren, J., Fan, K., Ma, S., & Wu, C. (2022). Sparse group variable selection for gene–environment interactions in the longitudinal study. Genetic epidemiology, 46(5-6), 317-340.
#'
#' Zhou, F., Ren, J.,  Li, G., Jiang, Y., Li, X., Wang, W. and Wu, C. (2019). Penalized Variable Selection for Lipid-Environment Interactions in a Longitudinal Lipidomics Study.
#' {\emph{Genes}, 10(12), 1002} \doi{10.3390/genes10121002}
#'
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2020). roben: Robust Bayesian Variable Selection for Gene-Environment Interactions.
#' R package version 0.1.1. \url{https://CRAN.R-project.org/package=roben}
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2021). Gene–Environment Interaction: a Variable Selection Perspective.
#' {\emph{Epistasis. Methods in Molecular Biology.} 2212:191–223} \doi{10.1007/978-1-0716-0947-7_13}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y. and Wu, C. (2020) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 39: 617– 638} \doi{10.1002/sim.8434}
#'
#' Ren, J., Zhou, F., Li, X., Wu, C. and Jiang, Y. (2019) spinBayes: Semi-Parametric Gene-Environment Interaction via Bayesian Variable Selection.
#' R package version 0.1.0. \url{https://CRAN.R-project.org/package=spinBayes}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' {\emph{Statistics in Medicine}, 37:437–456} \doi{10.1002/sim.7518}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' {\emph{Statistics in Medicine}, 33(28), 4988–4998} \doi{10.1002/sim.6287}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' {\emph{Technical Report. Michigan State University.}}
#'
#' @seealso \code{\link{Blend}}
NULL
