#' Coverage indicators for a Blend object
#'
#' Calculate coverage indicators for varying and constant effects based on posterior credible intervals.
#'
#' @param x A fitted Blend object.
#' @param truth A list of true coefficient functions and constant coefficients, including the intercept.
#' @param level The credible interval level. The default is 0.95.
#' @param n.grid The number of equally spaced grid points used to evaluate coefficient function coverage. The default is 200.
#'
#' @return A named list of binary coverage indicators, with one element for each coefficient, including the intercept.
#'
#' @seealso \code{\link{Blend}}
#'
#' @examples
#'
#' data(dat)
#' fit = Blend(y, x, t, J, kn, degree)
#'
#' u = as.numeric(fit$basis$u)
#'
#' truth = list(
#'   gamma0 = function(u) 2+2*sin(2*pi*u),
#'   gamma1 = function(u) 2.5*exp(2.5*u-1),
#'   gamma2 = function(u) 3*u^2-2*u+2,
#'   gamma3 = function(u) -4*u^3+3,
#'   gamma4 = function(u) 3-2*u,
#'   gamma5 = 0.8,
#'   gamma6 = -1.2,
#'   gamma7 = 0.7,
#'   gamma8 = -1.1,
#'   gamma9 = 0,
#'   gamma10 = 0,
#'   gamma11 = 0,
#'   gamma12 = 0,
#'   gamma13 = 0,
#'   gamma14 = 0,
#'   gamma15 = 0
#' )
#'
#' out = Coverage(fit, truth)
#' print(out)
#'
#'
#' @export

Coverage = function(x, truth, level=0.95, n.grid=200){

  ## Model settings

  q = x$basis$q
  m = x$basis$m
  structural = x$structural

  ## Input checking

  if(!is.list(truth) || length(truth) != m+1L){
    stop("truth must be a list of length m+1, including the intercept.")
  }

  if(length(level) != 1L ||
     !is.numeric(level) ||
     !is.finite(level) ||
     level <= 0 || level >= 1){
    stop("level must be between 0 and 1.")
  }

  for(j in seq_along(truth)){

    if(!is.function(truth[[j]]) &&
       !(is.numeric(truth[[j]]) &&
         length(truth[[j]]) == 1L &&
         is.finite(truth[[j]]))){
      stop("Each element of truth must be a function or a finite scalar.")
    }
  }

  ## Identify varying and constant effects

  is.varying = vapply(truth, is.function, logical(1))

  needs.grid = any(is.varying) || !structural

  ## Credible interval probabilities

  probs = c((1-level)/2, 1-(1-level)/2)

  ## Posterior samples after burn-in

  keep = seq.int(x$burn.in+1L, x$iterations)

  GS.alpha = x$posterior$GS.alpha[keep, , drop=FALSE]

  if(structural){

    GS.beta = x$posterior$GS.beta[keep, , drop=FALSE]

  }

  if(needs.grid){

    GS.eta = x$posterior$GS.eta[keep, , drop=FALSE]

  }

  ## Evaluation grid and B-spline basis

  if(needs.grid){

    if(length(n.grid) != 1L ||
       !is.numeric(n.grid) ||
       !is.finite(n.grid) ||
       n.grid < 2L ||
       n.grid != floor(n.grid)){
      stop("n.grid must be an integer greater than one.")
    }

    u = as.numeric(x$basis$u)
    kn = x$basis$kn
    degree = x$basis$degree

    ## Evaluation grid

    u.plot = seq(min(u), max(u), length.out=n.grid)

    ## Reconstruct the original B-spline basis

    u.star = seq(0, 1, length.out=kn+2L)[-c(1L, kn+2L)]

    Knots = as.numeric(stats::quantile(u, u.star))

    pi.u = splines::bs(
      u,
      knots=Knots,
      intercept=TRUE,
      degree=degree
    )

    ## Evaluate the basis at the grid points

    pi.star = stats::predict(pi.u, newx=u.plot)

    ## Structural reparameterization

    if(structural){

      pi.star = cbind(1, pi.star[, -1, drop=FALSE])

    }
  }

  ## Coefficient function coverage

  varying.coverage = function(draws, true.values){

    CI = apply(
      draws, 2, stats::quantile,
      probs=probs,
      names=FALSE
    )

    covered = true.values >= CI[1L, ] &
      true.values <= CI[2L, ]

    as.integer(all(covered))
  }

  ## Constant parameter coverage

  constant.coverage = function(draws, true.value){

    CI = stats::quantile(
      draws,
      probs=probs,
      names=FALSE
    )

    as.integer(
      true.value >= CI[1L] &
        true.value <= CI[2L]
    )
  }

  ## Extract posterior coefficient draws

  coefficient.draws = function(j){

    ## Intercept

    if(j == 0L){

      return(GS.alpha)

    }

    ## Structural model

    if(structural){

      index = seq.int(
        (j-1L)*(q-1L)+1L,
        j*(q-1L)
      )

      return(
        cbind(
          GS.beta[, j, drop=FALSE],
          GS.eta[, index, drop=FALSE]
        )
      )

    }else{

      ## Non-structural model

      index = seq.int(
        (j-1L)*q+1L,
        j*q
      )

      return(
        GS.eta[, index, drop=FALSE]
      )
    }
  }

  ## Calculate coverage

  coverage = vector("list", m+1L)

  for(j in 0:m){

    true.effect = truth[[j+1L]]

    if(is.varying[j+1L] || !structural){

      ## Coefficient function coverage

      if(is.varying[j+1L]){

        true.values = true.effect(u.plot)

      }else{

        true.values = rep(true.effect, n.grid)

      }

      if(!is.numeric(true.values) ||
         length(true.values) != n.grid ||
         any(!is.finite(true.values))){
        stop(paste("Invalid true function values for coefficient", j))
      }

      coef.draws = coefficient.draws(j)

      gamma.draws = coef.draws %*% t(pi.star)

      coverage[[j+1L]] = varying.coverage(
        gamma.draws,
        true.values
      )

    }else{

      ## Constant parameter coverage

      if(j == 0L){

        coverage[[j+1L]] = constant.coverage(
          GS.alpha[, 1L],
          true.effect
        )

      }else{

        coverage[[j+1L]] = constant.coverage(
          GS.beta[, j],
          true.effect
        )
      }
    }
  }

  ## Output

  if(!is.null(names(truth))){

    names(coverage) = names(truth)

  }else{

    names(coverage) = paste0("gamma", 0:m)

  }

  return(coverage)
}
