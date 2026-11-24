Selection_Sparse = function(obj, burn.in=obj$burn.in){

  ## Model settings

  BI = if(is.null(burn.in)) 0 else burn.in

  structural = obj$structural
  m = obj$basis$m

  ## Posterior samples

  GS.beta = obj$posterior$GS.beta
  GS.phi = obj$posterior$GS.phi

  ## Remove burn-in

  if(BI>0){

    if(structural){
      GS.beta = GS.beta[-seq_len(BI), , drop=FALSE]
    }

    GS.phi = GS.phi[-seq_len(BI), , drop=FALSE]
  }

  ## Posterior inclusion probabilities

  if(structural){

    Selectbeta = colMeans(GS.beta != 0)
    Selecteta = colMeans(GS.phi)

  }else{

    Selectbeta = rep(0, m)
    Selecteta = colMeans(GS.phi)

  }

  ## Median probability model

  MPM.V = which(Selecteta > 0.5)
  MPM.C = setdiff(which(Selectbeta > 0.5), MPM.V)

  ## Summary

  numb = matrix(
    c(length(MPM.C), length(MPM.V)),
    ncol=1,
    dimnames=list(c("Constant effect", "Varying effect"), "#")
  )

  ## Selected variable indices

  Main = if(length(MPM.C)>0){
    stats::setNames(MPM.C, MPM.C)
  }else{
    NULL
  }

  Varying = if(length(MPM.V)>0){
    stats::setNames(MPM.V, MPM.V)
  }else{
    NULL
  }

  ## Output

  sel = list(Constant=Main, Varying=Varying)

  method = "Median Probability Model (MPM)"

  out = list(method=method, indices=sel, summary=numb)

  out
}
