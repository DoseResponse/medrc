vcov.metadrc <- function(object, ...){
  return(object$vcov)  
}

df.residual.metadrc <- function(object, ...){
  df.resid <- object$rma$k.eff - object$rma$p.eff
  return(df.resid)
}