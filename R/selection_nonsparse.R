Selection_Nonsparse=function(obj, burn.in=obj$burn.in, prob=0.95){
  BI = ifelse(is.null(burn.in), 0, burn.in)
  GS.beta = obj$posterior$GS.beta
  GS.eta = obj$posterior$GS.eta
  q = obj$basis$q
  m = obj$basis$m
  if(BI>0){
    GS.beta = GS.beta[-c(1:BI),]
    GS.eta = GS.eta[-(1:BI),]
  }

  if(is.null(GS.beta)){
    Selectbeta = 0
    Selecteta = Selection.CI(GS.eta, q, prob)
  }else{
    Selectbeta = Selection.CI(GS.beta, 1, prob)
    Selecteta = Selection.CI(GS.eta, max(1, q-1), prob)
  }

  CIM.V = which(Selecteta > 0)
  CIM.C = setdiff(which(Selectbeta > 0), CIM.V)

  numb = matrix(c(length(CIM.C), length(CIM.V)), ncol=1,
                  dimnames=list(c("Constant effect", "Varying effect"), "#"))

                Var.names = 1:m
                if(length(CIM.C)>0){
                  Main = CIM.C
                  names(Main) = Var.names[CIM.C]
                }else{
                  Main = NULL
                }

                if(length(CIM.V)>0){
                  Varying = CIM.V
                  names(Varying) = Var.names[CIM.V]
                }else{
                  Varying = NULL
                }

                sel = list(Constant=Main, Varying=Varying)

                method = paste(prob*100,"% credible interval", sep = "")
                out = list(method=method, indices=sel, summary=numb)
                out
}
