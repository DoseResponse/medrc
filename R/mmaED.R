mmaED <- function(object, ..., respLev, ic = c("AIC", "BIC"), interval = c("none", "buckland", "kang"), level = 0.95, marginal=FALSE, nGQ=5, rfinterval=c(0, 1000)){
  interval <- match.arg(interval)
  ic <- match.arg(ic)
  lllist <- list(object, ...)
  ismedrc <- sapply(lllist, function(x) inherits(x, "medrc") | inherits(x, "glsdrc"))
  mllist <- lllist[ismedrc]  
  Call <- match.call()
  Call$respLev <- NULL
  Call$ic <- NULL
  Call$interval <- NULL
  Call$level <- NULL
  Call$marginal <- NULL
  Call$nGQ <- NULL
  Call$rfinterval <- NULL
  mnames <- as.character(Call[-1L])[ismedrc]
  # number of curves
  if (any(ncol(mllist[[1]]$parmMat) != sapply(mllist, function(x) ncol(x$parmMat)))) stop("Number of curves are not the same for all models!")
  # model weights
  if (ic == "AIC") inc <- sapply(mllist, AIC)
  if (ic == "BIC") inc <- sapply(mllist, BIC)
  d <- inc-min(inc)
  wi <- exp(-0.5*d)/sum(exp(-0.5*d))
  names(wi) <- mnames
  
  # calculate ED per model
  if (marginal == TRUE){
    edlist <- lapply(mllist, function(x) EDmarg(x, respLev=respLev, display=FALSE, nGQ=nGQ, rfinterval=rfinterval))
  } else {
    edlist <- lapply(mllist, function(x) ED(x, respLev=respLev, display=FALSE))
  }
  edm <- rbind(sapply(edlist, function(x) x[,1]))
  colnames(edm) <- mnames
  
  # model-averaged ED
  edma <- apply(t(t(edm) * wi), 1, sum)
  
  if (interval[1] == "none"){
    retMat <- cbind("Estimate" = edma)
  }  
  
  if (interval[1] == "buckland"){
    edse <- rbind(sapply(edlist, function(x) x[,2]))
    sema <- apply(sqrt(edse^2 + (edm - apply(edm, 1, mean))^2), 1, function(x) sum(x*wi))
    quant <- qnorm(1 - (1 - level)/2) * sema
    retMat <- as.matrix(cbind(edma, sema, edma - quant, edma + quant))
    colnames(retMat) <- c("Estimate", "SE", "Lower", "Upper")
  }
  
  if (interval[1] == "kang"){
    if (marginal == TRUE){
      edclist <- lapply(mllist, function(x) EDmarg(x, respLev=respLev, interval="delta", display=FALSE, nGQ=nGQ, rfinterval=rfinterval))
    } else {
      edclist <- lapply(mllist, function(x) ED(x, respLev=respLev, interval="delta", display=FALSE))
    }
    ll <- rbind(sapply(edclist, function(x) x[,3]))
    ul <- rbind(sapply(edclist, function(x) x[,4]))
    retMat <- as.matrix(cbind(edma, apply(ll, 1, function(x) sum(x*wi)), apply(ul, 1, function(x) sum(x*wi))))
    colnames(retMat) <- c("Estimate", "Lower", "Upper")
  }  
  
  out <- list()
  out$retMat <- retMat
  out$weights <- wi
  out$EDi <- edm
  if (interval[1] == "buckland") {out$SEi <- edse}  
  return(out)
}

