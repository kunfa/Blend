#' 95\% coverage for a Blend object with structural identification
#'
#' calculate 95\% coverage for varying effects and constant effects under example data
#'
#' @param x Blend object.
#' @usage Coverage(x)
#' @return coverage
#' @seealso \code{\link{Blend}}
#'
#' @examples
#' data(dat)
#' fit = Blend(y,x,t,J,kn,degree)
#' Coverage(fit)
#'
#' @export
Coverage = function(x){
  q = x$basis$q
  kn = x$basis$kn
  degree = x$basis$degree
  t = x$basis$t
  BI = x$burn.in
  m = x$basis$m
  u = x$basis$u
  n = length(u)
  u.plot = seq(min(u), max(u), length.out = n)
  gamma0 = 2+2*sin(u.plot*2*pi)
  gamma2 = 3*u.plot^2-2*u.plot+2
  gamma1 = 2.5*exp(2.5*u.plot-1)
  gamma3= -4*u.plot^3+3
  gamma4 = 3-2*u.plot
  gamma5 = 0.8
  gamma6 = -1.2
  gamma7 = 0.7
  gamma8 = -1.1
  u.star = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots.star = as.numeric(stats::quantile(u.plot, u.star))
  pi.star = splines::bs(u.plot, knots=Knots.star, intercept=TRUE, degree=degree)[,1:(q)]
  pi.star = cbind(1,pi.star[,-1])
  c1.C=rep(0,q)
  for (i in 1:q) {
    c1.C[i]=stats::quantile(x$posterior$GS.alpha[5001:10000,i],0.5)
  }

  c2.C=rep(0,dim(x$posterior$GS.beta)[2])
  for (i in 1:dim(x$posterior$GS.beta)[2]) {
    c2.C[i]=stats::quantile(x$posterior$GS.beta[5001:10000,i],0.5)
  }

  c3.C=rep(0,dim(x$posterior$GS.eta)[2])
  for (i in 1:dim(x$posterior$GS.eta)[2]) {
    c3.C[i]=stats::quantile(x$posterior$GS.eta[5001:10000,i],0.5)
  }
  C = matrix(rbind(c2.C, matrix(c3.C, ncol = m)), ncol = m*q)
  C = matrix(C,nrow=q)
  coeffmatrix.C=as.matrix(cbind(c1.C,C))
  c1.C=rep(0,q)
  for (i in 1:q) {
    c1.C[i]=stats::quantile(x$posterior$GS.alpha[5001:10000,i],0.025)
  }

  c2.C=rep(0,dim(x$posterior$GS.beta)[2])
  for (i in 1:dim(x$posterior$GS.beta)[2]) {
    c2.C[i]=stats::quantile(x$posterior$GS.beta[5001:10000,i],0.025)
  }

  c3.C=rep(0,dim(x$posterior$GS.eta)[2])
  for (i in 1:dim(x$posterior$GS.eta)[2]) {
    c3.C[i]=stats::quantile(x$posterior$GS.eta[5001:10000,i],0.025)
  }
  C = matrix(rbind(c2.C, matrix(c3.C, ncol = m)), ncol = m*q)
  C = matrix(C,nrow=q)
  coeffmatrix.C1=as.matrix(cbind(c1.C,C))
  c1.C=rep(0,q)
  for (i in 1:q) {
    c1.C[i]=stats::quantile(x$posterior$GS.alpha[5001:10000,i],0.975)
  }

  c2.C=rep(0,dim(x$posterior$GS.beta)[2])
  for (i in 1:dim(x$posterior$GS.beta)[2]) {
    c2.C[i]=stats::quantile(x$posterior$GS.beta[5001:10000,i],0.975)
  }

  c3.C=rep(0,dim(x$posterior$GS.eta)[2])
  for (i in 1:dim(x$posterior$GS.eta)[2]) {
    c3.C[i]=stats::quantile(x$posterior$GS.eta[5001:10000,i],0.975)
  }
  C = matrix(rbind(c2.C, matrix(c3.C, ncol = m)), ncol = m*q)
  C = matrix(C,nrow=q)
  coeffmatrix.C2=as.matrix(cbind(c1.C,C))
  ind = 1:5
  gamma.var.RBGLSS =pi.star %*% coeffmatrix.C[,ind]
  gamma.var.RBGLSS1 =pi.star %*% coeffmatrix.C1[,ind]
  gamma.var.RBGLSS2 =pi.star %*% coeffmatrix.C2[,ind]

  CI1 = rep(0,n)
  for(i in 1:n){
    CI1[i] = ifelse(gamma0[i]>=gamma.var.RBGLSS1[i,1]&gamma0[i]<=gamma.var.RBGLSS2[i,1],1,0)
  }

  C1 = sum(CI1)
  C1 = ifelse(C1!=0,1,0)
  CI2 = rep(0,n)
  for(i in 1:n){
    CI2[i] = ifelse(gamma1[i]>=gamma.var.RBGLSS1[i,2]&gamma1[i]<=gamma.var.RBGLSS2[i,2],1,0)
  }
  C2 = sum(CI2)
  C2 = ifelse(C2!=0,1,0)
  CI3 = rep(0,n)
  for(i in 1:n){
    CI3[i] = ifelse(gamma2[i]>=gamma.var.RBGLSS1[i,3]&gamma2[i]<=gamma.var.RBGLSS2[i,3],1,0)
  }
  C3 = sum(CI3)
  C3 = ifelse(C3!=0,1,0)
  CI4 = rep(0,n)
  for(i in 1:n){
    CI4[i] = ifelse(gamma3[i]>=gamma.var.RBGLSS1[i,4]&gamma3[i]<=gamma.var.RBGLSS2[i,4],1,0)
  }
  C4 = sum(CI4)
  C4 = ifelse(C4!=0,1,0)

  #100%
  CI5 = rep(0,n)
  for(i in 1:n){
    CI5[i] = ifelse(gamma4[i]>=gamma.var.RBGLSS1[i,5]&gamma4[i]<=gamma.var.RBGLSS2[i,5],1,0)
  }
  C5 = sum(CI5)
  C5 = ifelse(C5!=0,1,0)

  #100%

  C6 = ifelse(gamma5>=coeffmatrix.C1[1,6]& gamma5<=coeffmatrix.C2[1,6],1,0)

  C7 = ifelse(gamma6>=coeffmatrix.C1[1,7]&gamma6<=coeffmatrix.C2[1,7],1,0)

  C8 = ifelse(gamma7>=coeffmatrix.C1[1,8]&gamma7<=coeffmatrix.C2[1,8],1,0)

  C9 = ifelse(gamma8>=coeffmatrix.C1[1,9]&gamma8<=coeffmatrix.C2[1,9],1,0)

  coverage = list(C1,C2,C3,C4,C5,C6,C7,C8,C9)
  coverage
}
