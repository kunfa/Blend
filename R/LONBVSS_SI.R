LONBVSS_SI <- function(y,X1,X2,X3,J,q,max.steps,robust,quant,sparse){
  n = nrow(X2)
  n1 =length(J)
  J2 = c()
  for(i in 1:n1){
  z_i = cbind(rep(1,J[i]),c(1:J[i]))
  c = ncol(z_i)

  ata = stats::rnorm(c,0,1)

  J1 = kronecker(t(ata),z_i)
  J2 = rbind(J2,J1)
  }
  J3 = matrix(c(1,0,0,0,0,1,0,0,0,0,0,1),nrow=4)
  C = J2%*%J3

  m = ncol(X2)
  q1 = ncol(X1)
  hatAlpha = rep(1,q1)     ## coeff for intercept
  hatBeta =rep(1,m) ## coeff for constant ## coeff for varying part
  hatEta = matrix(rep(1,m*(q-1)),ncol=m)
  hatAta = rep(1,3)
  hatBeta1 = rep(0,m)
  hatEta1 = matrix(rep(0,m*(q-1)),ncol=m)
  hatBeta1 = hatBeta1+10^-5
  hatEta1 = hatEta1+10^-5
  xi1 = (1-2*quant)/(quant*(1-quant))
  xi2 = sqrt(2/(quant*(1-quant)))
  hatTau = 1
  hatV = rep(1,n)
  hatEtaSq1 = 1
  hatEtaSq2 = 1
  hatEtaSq3 = 1
  hatEtaSq4 = 1
  hatSg1 = rep(1, m)
  hatSg2 = rep(1, m)
  hatSg3 = 1
  hatSg4 = 1
  hatPiBeta = 0.5
  hatPiEta = 0.5
  hatPi3 = 0.5
  hatPi4 = 0.5
  invSigAlpha0 = diag(10^-3, q1)
  r1=r2=a=b=sh1=sh0=1
  r3 =1
  r4 =1
  hatEta2 = rep(1,m*(q-1))
  hatInvTauSq1=rep(0.1,m)
  hatInvTauSq22=rep(0.1,m)
  hatInvTauSq3=0.1
  hatInvTauSq4=0.1
  hatLambdaSqStar1=1
  hatLambdaSqStar2=1
  hatLambdaSqStar3=1
  hatLambdaSqStar4=1
  hatSigmaSq=1
  a0=1
  b0=1.5
  aStar=1
  bStar=1.5
  a1=1
  b1=1.5
  a2=1
  b2=1.5
  alpha=0.2
  gamma=0.1
  mu0=1
  nu0=1
  muStar=1
  nuStar=1
  mu1=1
  nu1=1
  mu2=1
  nu2=1
  max.steps=10000
  invSigAlpha1 = rep(10^-3,q1)
  debugging=FALSE

  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)
  if(robust){
    fit=switch (sparse,
                "TRUE" = RBVSS_SI(y,X1,X2,X3,C,max.steps,q,hatBeta,hatEta,hatAlpha,hatAta,hatTau,hatV,hatSg1,hatSg2,hatSg3,hatSg4,invSigAlpha0,hatPiBeta,hatPiEta,hatPi3,hatPi4,hatEtaSq1,hatEtaSq2,hatEtaSq3,hatEtaSq4,xi1,xi2,r1,r2,r3,r4,a,b,sh1,sh0,progress),
                "FALSE" = RBV_SI(y,X1,X2,X3,C,max.steps,q,hatBeta1,hatEta1,hatAlpha,hatAta,hatTau,hatV,hatSg1,hatSg2,hatSg3,hatSg4,invSigAlpha0,hatPi3,hatPi4,hatEtaSq1,hatEtaSq2,hatEtaSq3,hatEtaSq4,xi1,xi2,r1,r2,r3,r4,a,b,sh1,sh0,progress)
    )
  }else{
    fit=switch (sparse,
                "TRUE" = BVSS_SI(y,X1,X2,X3,C,m,q,max.steps,hatAlpha,hatBeta,hatAta,hatEta2,invSigAlpha1,hatInvTauSq1,hatInvTauSq22,hatInvTauSq3,hatInvTauSq4,hatPiBeta,hatPiEta,hatPi3,hatPi4,hatLambdaSqStar1,hatLambdaSqStar2,hatLambdaSqStar3,hatLambdaSqStar4,hatSigmaSq,a0,b0,aStar,bStar,a1,b1,a2,b2,alpha,gamma,mu0,muStar,mu1,mu2,nu0,nuStar,nu1,nu2,progress) ,
                "FALSE" = BV_SI(y,X1,X2,X3,C,m,q,max.steps,hatAlpha,hatBeta,hatEta2,hatAta,invSigAlpha1,hatInvTauSq1,hatInvTauSq22,hatInvTauSq3,hatInvTauSq4,hatPi3,hatPi4,hatLambdaSqStar1,hatLambdaSqStar2,hatLambdaSqStar3,hatLambdaSqStar4,hatSigmaSq,a0,b0,aStar,bStar,a1,b1,a2,b2,alpha,gamma,mu1,mu2,nu1,nu2,progress)
    )
  }
 if(sparse){
   GS.phi = fit$GS.tRsRs
 }else{
   GS.phi = NULL
 }
  out = list(posterior =list(GS.alpha = fit$GS.alpha,GS.phi = GS.phi,GS.beta = fit$GS.beta,GS.eta = fit$GS.eta))
  out

}
