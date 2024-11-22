#' plot a Blend object
#'
#' plot the identified varying effects
#'
#' @param x Blend object.
#' @param sparse sparsity.
#' @param prob probability for credible interval, between 0 and 1. e.g. prob=0.95 leads to 95\% credible interval
#' @usage plot_Blend(x, sparse, prob=0.95)
#' @return plot
#' @seealso \code{\link{Blend}}
#'
#' @examples
#' data(dat)
#' fit = Blend(y,x,t,J,kn,degree)
#' plot_Blend(fit,sparse=TRUE)
#'
#' @export
plot_Blend=function(x, sparse,prob=0.95){

  q = x$basis$q
  kn = x$basis$kn
  degree = x$basis$degree
  t = x$basis$t
  BI = x$burn.in
  s = x$basis$m
  u = x$basis$u

  adj=0.1

  varName=1:(s+1)

  xlab="Time"
  lt = (1-prob)/2; ut= 1-lt
  n = length(u)
  u.plot = seq(min(u), max(u), length.out = n)
  u.star = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots.star = as.numeric(stats::quantile(u.plot, u.star))
  pi.star = splines::bs(u.plot, knots=Knots.star, intercept=TRUE, degree=degree)#[,1:(q)]
  pi.star = cbind(1,pi.star[,-1])

  selected = selection(x,sparse=sparse)
  ind.V = selected$indices$Varying

  if(length(ind.V)>0){
    temp = matrix(rbind(x$posterior$GS.beta, matrix(x$posterior$GS.eta, ncol = s)), ncol = s*q)
  }

  for(j in c(0, ind.V)){
    if(j==0){
      coeff.mat = pi.star %*% t(x$posterior$GS.alpha[-c(1:BI),])
    }else{
      last = j*q; first = last-q+1
      coeff.mat = pi.star %*% t(temp[-c(1:BI),first:last])
    }
    pe = apply(coeff.mat, 1, stats::median)

    LL = apply(coeff.mat, 1, function(t) stats::quantile(t, lt))
    UL = apply(coeff.mat, 1, function(t) stats::quantile(t, ut))

    data <- data.frame(t,pe,LL,UL)
    p <- ggplot2::ggplot(data, ggplot2::aes(x=t)) +
      ggplot2::geom_line(ggplot2::aes(y = pe), color = "gray35", size=1) +
      ggplot2::geom_line(ggplot2::aes(y = LL), color="steelblue", linetype=2, size=1, alpha=0.9) +
      ggplot2::geom_line(ggplot2::aes(y = UL), color="steelblue", linetype=2, size=1, alpha=0.9) +
      ggplot2::labs(title=varName[j+1],
                    x= xlab,
                    y= bquote(~ beta[.(j)] ~ (t)))
    print(p)
  }
}
