LONBVSS <- function(y,X1,X2,J,m,q,max.steps,robust,quant,sparse){
  n = nrow(X2)
  J2 = c()
  n1 =length(J)
  for(i in 1:n1){
    z_i = cbind(rep(1,J[i]),c(1:J[i]))
    c = ncol(z_i)

    ata = stats::rnorm(c,0,1)

    J1 = kronecker(t(ata),z_i)
    J2 = rbind(J2,J1)
  }
  J3 = matrix(c(1,0,0,0,0,1,0,0,0,0,0,1),nrow=4)
  C = J2%*%J3
  q1 = ncol(X1)
  hatAlpha = rep(1,q1)

  hatBeta = matrix(rep(1,m*q),ncol=m)
  hatBeta1 = matrix(rep(0,m*q),ncol=m)
  hatAta = rep(1,3)
  hatBeta1 = hatBeta1+10^-5

  xi1 = (1-2*quant)/(quant*(1-quant))
  xi2 = sqrt(2/(quant*(1-quant)))
  hatTau = 1
  hatV = rep(1,n) # rgamma(n, shape=1, rate=hatTau)
  hatEtaSq1 = 1
  hatEtaSq2 =1
  hatEtaSq3 = 1
  hatSg1 = rep(1, m)
  hatSg2=1
  hatSg3=1
  hatPiBeta = 0.5
  hatPi1=0.5
  hatPi2=0.5
  invSigAlpha0 = diag(10^-3, q1)
  r1=r2=a=b=sh1=sh0=1
  r3 =1
  r4 =1
  hatBeta2 = rep(1,m*q)
  hatInvTauSq1=rep(0.1,m)
  hatInvTauSq2=0.1
  hatInvTauSq3=0.1
  hatLambdaSqStar1=1
  hatLambdaSq1=1
  hatLambdaSq2=1
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
  debugging=FALSE

  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)
  if(robust){
    fit=switch (sparse,
                "TRUE" = RBVSS(X2,y,X1,C,m,q,max.steps,hatAlpha,hatBeta,hatAta,hatTau,hatV,hatSg1,hatSg2,hatSg3,invSigAlpha0,hatPiBeta,hatPi1,hatPi2,hatEtaSq1,hatEtaSq2,hatEtaSq3,xi1,xi2,r1,r2,r3,a,b,sh1,sh0,progress),
                "FALSE" = RBV(X2,y,X1,C,m,q,max.steps,hatAlpha,hatBeta1,hatAta,hatTau,hatV,hatSg1,hatSg2,hatSg3,invSigAlpha0,hatPi1,hatPi2,hatEtaSq1,hatEtaSq2,hatEtaSq3,xi1,xi2,r1,r2,r3,a,b,sh1,sh0,progress)
    )
  }else{
    fit=switch (sparse,
                "TRUE" = BVSS(y,X1,X2,C,m,q,max.steps,hatAlpha,hatBeta2,hatAta,hatInvTauSq1,hatInvTauSq2,hatInvTauSq3,invSigAlpha0,hatPiBeta,hatPi1,hatPi2,hatLambdaSqStar1,hatLambdaSq1,hatLambdaSq2,hatSigmaSq,aStar,bStar,a1,b1,a2,b2,alpha,gamma,mu0,nu0,mu1,nu1,mu2,nu2,progress),
                "FALSE" = BV(X2,y,X1,C,m,q,max.steps,hatBeta2,hatAlpha,hatAta,hatInvTauSq1,hatInvTauSq2,hatInvTauSq3,invSigAlpha0,hatPi1,hatPi2,hatLambdaSqStar1,hatLambdaSq1,hatLambdaSq2,hatSigmaSq,aStar,bStar,a1,b1,a2,b2,alpha,gamma,mu1,nu1,mu2,nu2,progress)
    )
  }
  if(sparse){
    GS.phi = fit$GS.tRsRs
  }else{
    GS.phi = NULL
  }
  out = list( posterior = list(GS.alpha = fit$GS.alpha,GS.phi = GS.phi,GS.beta = NULL,GS.eta = fit$GS.beta))
  out

}
