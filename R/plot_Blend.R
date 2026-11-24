#' Plot a Blend object
#'
#' Plot the intercept function and the identified varying effects with pointwise posterior credible intervals.
#'
#' @param x A fitted Blend object.
#' @param sparse Logical value indicating whether a sparse model is used.
#' @param prob probability for credible interval, between 0 and 1. e.g. prob=0.95 leads to 95\% credible interval.
#' @param n.grid The number of evaluation grid points. The default is 200.
#' @usage plot_Blend(x, sparse, prob=0.95, n.grid=200)
#' @return Plots of the intercept function and selected coefficient functions.
#' @seealso \code{\link{Blend}}
#'
#' @examples
#' data(dat)
#' fit = Blend(y,x,t,J,kn,degree)
#' plot_Blend(fit,sparse=TRUE)
#'
#' @export
plot_Blend = function(x, sparse, prob=0.95, n.grid=200){
  pe = LL = UL = NULL
  ## Model settings

  q = x$basis$q
  kn = x$basis$kn
  degree = x$basis$degree
  structural = x$structural

  t = x$basis$t
  u = as.numeric(x$basis$u)

  ## Credible interval probabilities

  probs = c((1-prob)/2, 0.5, 1-(1-prob)/2)

  ## Posterior samples after burn-in

  keep = seq.int(x$burn.in+1L, x$iterations)

  GS.alpha = x$posterior$GS.alpha[keep, , drop=FALSE]
  GS.eta = x$posterior$GS.eta[keep, , drop=FALSE]

  if(structural){

    GS.beta = x$posterior$GS.beta[keep, , drop=FALSE]

  }

  ## Evaluation grid

  t.plot = seq(min(t), max(t), length.out=n.grid)

  u.plot = stats::pnorm((t.plot-mean(t))/stats::sd(t))

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

  ## Variable selection

  selected = selection(x, sparse=sparse)

  ind.V = selected$indices$Varying

  ## Plot coefficient functions

  for(j in c(0, ind.V)){

    ## Extract posterior coefficient draws

    if(j == 0L){

      coef.draws = GS.alpha

    }else if(structural){

      index = seq.int(
        (j-1L)*(q-1L)+1L,
        j*(q-1L)
      )

      coef.draws = cbind(
        GS.beta[, j, drop=FALSE],
        GS.eta[, index, drop=FALSE]
      )

    }else{

      index = seq.int(
        (j-1L)*q+1L,
        j*q
      )

      coef.draws = GS.eta[, index, drop=FALSE]

    }

    ## Reconstruct posterior coefficient functions

    gamma.draws = coef.draws %*% t(pi.star)

    ## Posterior medians and credible intervals

    CI = apply(
      gamma.draws, 2, stats::quantile,
      probs=probs,
      names=FALSE
    )

    ## Plotting data

    data = data.frame(
      t=t.plot,
      pe=CI[2L, ],
      LL=CI[1L, ],
      UL=CI[3L, ]
    )

    ## Plot

    p = ggplot2::ggplot(data, ggplot2::aes(x=t)) +
      ggplot2::geom_line(
        ggplot2::aes(y=pe),
        color="gray35",
        linewidth=1
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y=LL),
        color="steelblue",
        linetype=2,
        linewidth=1,
        alpha=0.9
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y=UL),
        color="steelblue",
        linetype=2,
        linewidth=1,
        alpha=0.9
      ) +
      ggplot2::labs(
        title=as.character(j+1L),
        x="Time",
        y=bquote(~ beta[.(j)] ~ (t))
      )

    print(p)
  }

}
