Selection_Nonsparse = function(obj, burn.in=obj$burn.in, prob=0.95){

  ## Model settings

  BI = if(is.null(burn.in)) 0 else burn.in

  q = obj$basis$q
  m = obj$basis$m
  structural = obj$structural

  ## Posterior samples

  GS.beta = obj$posterior$GS.beta
  GS.eta = obj$posterior$GS.eta

  ## Remove burn-in

  if(BI>0){

    if(structural){
      GS.beta = GS.beta[-seq_len(BI), , drop=FALSE]
    }

    GS.eta = GS.eta[-seq_len(BI), , drop=FALSE]
  }

  ## Variable selection

  if(structural){

    Selectbeta = Selection.CI(GS.beta, 1, prob)
    Selecteta = Selection.CI(GS.eta, q-1, prob)

  }else{

    Selectbeta = rep(0, m)
    Selecteta = Selection.CI(GS.eta, q, prob)

  }

  ## Identify constant and varying effects

  CIM.V = which(Selecteta > 0)
  CIM.C = setdiff(which(Selectbeta > 0), CIM.V)

  ## Summary

  numb = matrix(
    c(length(CIM.C), length(CIM.V)),
    ncol=1,
    dimnames=list(c("Constant effect", "Varying effect"), "#")
  )

  ## Selected variable indices

  Main = if(length(CIM.C)>0){
    stats::setNames(CIM.C, CIM.C)
  }else{
    NULL
  }

  Varying = if(length(CIM.V)>0){
    stats::setNames(CIM.V, CIM.V)
  }else{
    NULL
  }

  ## Output

  sel = list(Constant=Main, Varying=Varying)

  method = paste0(prob*100, "% credible interval")

  out = list(method=method, indices=sel, summary=numb)

  out
}
