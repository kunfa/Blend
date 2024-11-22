Selection.CI <- function(GS.r, L, level){
  lt = (1-level)/2; ut= 1-lt
  limits = apply(GS.r, 2, stats::quantile, probs = c(lt, ut))
  temp = matrix(abs(sign(limits[1,]) + sign(limits[2,]))==2, nrow = L)
  (apply(temp, 2, sum)>=1)*1
}
