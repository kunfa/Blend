#' fit a robust Bayesian longitudinal regularized semi-parametric mixed model
#'
#' @keywords models
#' @param x the matrix of predictors for repeated measurements, including an intercept in the first column.
#' @param y the vector of repeated measurements of the response variable. The current version supports only continuous responses.
#' @param t the vector of observation time points corresponding to the repeated measurements of all subjects.
#' @param J the vector specifying the number of repeated measurements for each subject.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B spline basis.
#' @param iterations the number of MCMC iterations.
#' @param burn.in the number of iterations for burn-in. If NULL (default), the first half of the iterations are discarded as burn-in.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param structural logical flag. If TRUE, the coefficient functions with varying effects and constant effects will be penalized separately.
#' @return an object of class `Blend' is returned, which is a list with component:
#' \item{posterior}{the posteriors of coefficients.}
#' \item{coefficient}{the posterior estimates of the model coefficients.}
#' \item{burn.in}{the total number of burn-ins.}
#' \item{iterations}{the total number of iterations.}
#' \item{structural}{a logical value indicating whether the structural model is used.}
#' \item{basis}{a list containing the number of B-spline basis
#' functions (\code{q}), the number of interior knots (\code{kn}),
#' the degree of B-splines (\code{degree}), the original observation
#' times (\code{t}), the number of covariates excluding the intercept
#' (\code{m}), and the transformed observation times (\code{u}).}

#'
#' @details Consider the data model described in \code{\link{data}}:
#'
#' \deqn{Y_{ij} = \alpha_0(t_{ij}) +
#' \sum_{k=1}^{m}\beta_k(t_{ij})X_{ijk} +
#' Z_{ij}^\top\zeta_i + \epsilon_{ij}.}
#'
#' The coefficient functions are approximated using B-spline
#' basis expansions. When \code{structural=TRUE} (default),
#' an equivalent change of basis is performed automatically:
#'
#' \deqn{\beta_k(\cdot) \approx \gamma_{k1} +
#' \sum_{u=2}^{q} B_{ku}(\cdot)\gamma_{ku},}
#'
#' where \eqn{B_{ku}(\cdot)} denotes the B-spline basis function,
#' and \eqn{\gamma_{k1}} and
#' \eqn{(\gamma_{k2}, \ldots, \gamma_{kq})^\top}
#' represent the constant and varying components of the
#' coefficient function, respectively.
#' The number of basis functions is \eqn{q=kn+degree+1}.
#' By default, \code{kn=2} and \code{degree=2}.
#'
#' When \code{structural=TRUE} (default), constant and
#' time-varying effects are regularized separately.
#' Otherwise, each coefficient function is regularized
#' as a whole without separating constant and
#' time-varying effects.
#'
#' When \code{sparse=TRUE} (default), spike-and-slab priors
#' are imposed at the individual and/or group levels for
#' variable selection. Otherwise, Laplacian shrinkage priors
#' are used without inducing exact sparsity.
#'
#' When \code{robust=TRUE} (default), the error terms
#' \eqn{\epsilon_{ij}} are modeled using a Laplace distribution
#' with density
#'
#' \deqn{f(\epsilon_{ij}\mid\tau) =
#' \frac{\tau}{4}\exp\left\{-\frac{\tau}{2}
#' |\epsilon_{ij}|\right\},}
#'
#' where \eqn{\tau>0} and \eqn{2/\tau} is the scale parameter.
#'
#' When \code{robust=FALSE}, the error terms are modeled
#' using a Gaussian distribution.
#'
#' Posterior inference is conducted through Gibbs sampling.
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

Blend <- function(y,x,t,J,kn,degree,iterations=10000, burn.in=NULL, robust=TRUE, sparse=TRUE, structural=TRUE)
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
      out = LONBVSS_SI(y,X1,X2,X3,J,q,iterations,robust,quant=0.5,sparse)
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
      out = LONBVSS(y,X1,X2,J,m,q,iterations,robust,quant=0.5,sparse)
      INT = apply(out$posterior$GS.alpha[-c(1:BI),,drop=FALSE], 2, stats::median)
      VV = apply(out$posterior$GS.eta[-c(1:BI),,drop=FALSE], 2, stats::median)
      coeff = cbind(INT, matrix(VV, nrow = q))

    }
  colnames(coeff) = c("intercept",c(1:m))
  rownames(coeff) = paste("basis", 0:(q-1), sep="")
  basis = list(q=q, kn=kn, degree=degree, t=t,m = m,u=u)
  fit = list(
    posterior = out$posterior,
    coefficient = coeff,
    burn.in = BI,
    iterations = iterations,
    structural = structural
  )
  fit$basis = basis
  fit
}
