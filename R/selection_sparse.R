Selection_Sparse=function(obj, burn.in=obj$burn.in){
  BI = ifelse(is.null(burn.in), 0, burn.in)
  max_BI = obj$iterations - BI
  GS.beta = obj$posterior$GS.beta
  GS.eta = obj$posterior$GS.eta
  GS.phi = obj$posterior$GS.phi
  m = obj$basis$m
  if(BI>0){
    GS.beta = GS.beta[-c(1:BI),]
    GS.eta = GS.eta[-(1:BI),]
    GS.phi = GS.phi[-c(1:BI),]
  }
  Selectbeta = if(is.null(GS.beta)){ 0 }else{ apply(GS.beta, 2, function(t) sum(t!=0))}
  Selecteta = apply(GS.phi, 2, sum)

  MPM.V = which(Selecteta > max_BI/2)
  MPM.C = setdiff(which(Selectbeta > max_BI/2), MPM.V)

  numb = matrix(c(length(MPM.C), length(MPM.V)), ncol=1,
                dimnames=list(c("Constant effect", "Varying effect"), "#"))

  Var.names = 1:m
  if(length(MPM.C)>0){
    Main = MPM.C
    names(Main) = Var.names[MPM.C]
  }else{
    Main = NULL
  }

  if(length(MPM.V)>0){
    Varying = MPM.V
    names(Varying) = Var.names[MPM.V]
  }else{
    Varying = NULL
  }


  sel = list(Constant=Main, Varying=Varying)
  method = paste("Median Probability Model (MPM)", sep = "")

  out = list(method=method, indices=sel, summary=numb)
  out
}
