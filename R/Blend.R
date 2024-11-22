#' fit a Bayesian longitudinal regularized semi-parametric quantile mixed model
#'
#' @keywords models
#' @param x the matrix of repeated - measured predictors (genetic factors) with intercept. Each row should be an observation vector for each measurement.
#' @param y the vector of repeated - measured response variable. The current version of mixed only supports continuous response.
#' @param t the vector of scheduled time points.
#' @param J the vector of number of repeated measurement for each subject.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B spline basis.
#' @param iterations the number of MCMC iterations.
#' @param burn.in the number of iterations for burn-in.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param quant specify different quantiles when applying robust methods.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param structural logical flag. If TRUE, the coefficient functions with varying effects and constant effects will be penalized separately.
#' @return an object of class `Blend' is returned, which is a list with component:
#' \item{posterior}{the posteriors of coefficients.}
#' \item{coefficient}{the estimated coefficients.}
#' \item{burn.in}{the total number of burn-ins.}
#' \item{iterations}{the total number of iterations.}
#'
#' @details Consider the data model described in "\code{\link{data}}":
#' \deqn{Y_{ij} = \alpha_0(t_{ij})+\sum_{k=1}^{m}\beta_{k}(t_{ij})X_{ijk}+\boldsymbol{Z^\top_{ij}}\boldsymbol{\zeta_{i}}+\epsilon_{ij}.}
#' The basis expansion and changing of basis with B splines will be done automatically:
#' \deqn{\beta_{k}(\cdot)\approx \gamma_{k1} + \sum_{u=2}^{q}{B}_{ku}(\cdot)\gamma_{ku}}
#' where \eqn{B_{ku}(\cdot)} represents B spline basis. \eqn{\gamma_{k1}} and \eqn{(\gamma_{k2}, \ldots, \gamma_{kq})^\top} correspond to the constant and varying parts of the coefficient functional, respectively.
#' q=kn+degree+1 is the number of basis functions. By default, kn=degree=2. User can change the values of kn and degree to any other positive integers.
#' When `structural=TRUE`(default), the coefficient functions with varying effects and constant effects will be penalized separately. Otherwise, the coefficient functions with varying effects and constant effects will be penalized together.
#'
#' When `sparse="TRUE"` (default), spike-and-slab priors are imposed on individual and/or group levels to identify important constant and varying effects. Otherwise, Laplacian shrinkage will be used.
#'
#' When `robust=TRUE` (default), the distribution of \eqn{\epsilon_{ij}} is defined as a Laplace distribution with density.
#'
#' \eqn{
#' f(\epsilon_{ij}|\theta,\tau) = \theta(1-\theta)\exp\left\{-\tau\rho_{\theta}(\epsilon_{ij})\right\}
#' }, (\eqn{i=1,\dots,n,j=1,\dots,J_{i} }), which leads to a Bayesian formulation of quantile regression. If `robust=FALSE`, \eqn{\epsilon_{ij}} follows a normal distribution.

#' \cr
#'
#' Please check the references for more details about the prior distributions.
#'
#' @seealso \code{\link{data}}
#'
#' @examples
#' data(dat)
#'
#' ## default method

#' fit = Blend(y,x,t,J,kn,degree)
#' fit$coefficient
#'
#' \donttest{
#' ## alternative: robust non-structural
#' fit = Blend(y,x,t,J,kn,degree, structural=FALSE)
#' fit$coefficient
#'
#' ## alternative: non-robust structural
#' fit = Blend(y,x,t,J,kn,degree, robust=FALSE)
#' fit$coefficient
#'
#' ## alternative: non-robust non-structural
#' fit = Blend(y,x,t,J,kn,degree, robust=FALSE, structural=FALSE)
#' fit$coefficient
#'
#' }
#'
#' @export

Blend <- function(y,x,t,J,kn,degree,iterations=10000, burn.in=NULL, robust=TRUE, quant=0.5, sparse="TRUE", structural=TRUE)
{
  x = as.matrix(x); kn = as.integer(kn); degree = as.integer(degree)
  m = ncol(x)-1
  u = as.matrix(stats::pnorm((t-mean(t))/stats::sd(t)))
  n = length(u)
  q = kn + degree + 1
  u.k = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots = as.numeric(stats::quantile(u, u.k))
  pi.u = splines::bs(u, knots=Knots, intercept=TRUE, degree=degree)[,1:(q)]

  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
    BI = floor(iterations/2)
  }else if(burn.in>=1){
    BI = as.integer(burn.in)
  }else{
    stop("burn.in must be a positive integer.")
  }
  if(iterations<=BI) stop("iterations must be larger than burn.in.")
  if(structural){
      pi.u = cbind(1,pi.u[,-1])
      Pi = apply(x, 2, function(x) x*pi.u)
      X1 = matrix(Pi[,1], nrow=n)             # intercept
      X2 = Pi[c(1:n), -1]                     # constant part
      X3 = matrix(Pi[-c(1:n), -1], nrow=n)    # varying part
      out = LONBVSS_SI(y,X1,X2,X3,J,q,iterations,robust,quant,sparse)
      INT = apply(out$posterior$GS.alpha[-c(1:BI),,drop=FALSE], 2, stats::median)
      CC = apply(out$posterior$GS.beta[-c(1:BI),,drop=FALSE], 2, stats::median)
      VV = apply(out$posterior$GS.eta[-c(1:BI),,drop=FALSE], 2, stats::median)
      coeff = cbind(INT, rbind(CC, matrix(VV, nrow = q-1)))
    }else{
      xx1 = as.data.frame(matrix(0, n, (m+1)*q))
      for(j in 1:(m+1)){
        last = j*q; first = last-q+1
        xx1[,first:last] = pi.u*x[,j]
      }
      xx1 = as.matrix(xx1)
      X2=xx1[,-(1:q)]
      X1=pi.u
      out = LONBVSS(y,X1,X2,J,m,q,iterations,robust,quant,sparse)
      INT = apply(out$posterior$GS.alpha[-c(1:BI),,drop=FALSE], 2, stats::median)
      VV = apply(out$posterior$GS.eta[-c(1:BI),,drop=FALSE], 2, stats::median)
      coeff = cbind(INT, matrix(VV, nrow = q))

    }
  colnames(coeff) = c("intercept",c(1:m))
  rownames(coeff) = paste("basis", 0:(q-1), sep="")
  basis = list(q=q, kn=kn, degree=degree, t=t,m = m,u=u)
  fit = list(posterior = out$posterior, coefficient=coeff, burn.in = BI, iterations=iterations)
  fit$basis = basis
  fit
}
