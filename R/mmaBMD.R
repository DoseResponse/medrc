#' BMD estimation for averaged medrc or glsdrc models
#' 
#' Estimates the benchmark dose for averaged fixed effects of several medrc or glsdrc objects
#' 
#' @param object an medrc object
#' @param ... further medrc objects
#' @param respLev a numeric vector containing the benchmark response levels
#' @param ic the information criterion used for model averaging
#' @param interval character string specifying the type of confidence intervals to be supplied.
#' @param level confidence level, default at 0.95
#' @param bmd benchmark dose estimation (smallest dose resulting in a probability of an abnormal response)
#' @param background probability of an abnormal response
#' @param marginal logical. If TRUE, marginalized BMD estimates are obtained. See BMDmarg for further information.
#' @param nGQ integer. Specifies the number nof nodes for Gauss-Hermite quadrature.
#' @param rfinterval numeric vector. Interval for root finding (uniroot) to search for ED values.
#' 
#' @keywords htest

mmaBMD <- function(object, ..., respLev, ic = c("AIC", "BIC"), interval = c("none", "buckland", "kang"), level = 0.95, bmd = c("additional", "extra"), background = 0.05, marginal=FALSE, nGQ=5, rfinterval=c(0, 1000)){
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
  Call$background <- NULL
  Call$bmd <- NULL
  Call$dmlist <- NULL
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
  
  # adjusting respLev for BMD
  if (bmd[1] == "extra") respLev <- respLev * (1-background)  
  lenRL <- length(respLev)
  respLevMat <- lapply(1:lenRL, function(i) sapply(1:length(mllist),function(x){
    objectFit <- mllist[[x]]
    cobj <- objectFit$parmMat
    prnames <- rownames(cobj)
    cvals <- cobj[prnames == "c", ]
    if (length(cvals) == 0) {cvals <- 0}  # setting c if missing (not estimated)
    dvals <- cobj[prnames == "d", ]
    if (length(cvals) == 1 & length(dvals) > 1) cvals <- rep(cvals, length(dvals))
    if (length(dvals) == 1 & length(cvals) > 1) dvals <- rep(dvals, length(cvals))
    if (any(cvals > dvals)) {
      tempd <- apply(cbind(cvals, dvals), 1, function(x) sort(x))
      cvals <- tempd[1,]
      dvals <- tempd[2,]
    }
    if (class(objectFit)[1] == "medrc"){
      varcorr <- VarCorr(objectFit)
      return(100 * (qnorm(1 - background) - qnorm(1-(background+respLev[i]/100))) * as.numeric(varcorr[attr(varcorr, "dimnames")[[1]] %in% "Residual", 2]) / (dvals - cvals))
    }
    if (class(objectFit)[1] == "glsdrc"){
      return(100 * (qnorm(1 - background) - qnorm(1-(background+respLev[i]/100))) * objectFit$fit$sigma / (dvals - cvals))     
    }
  }))
  
  # ED estimation
  edl <- lapply(1:length(mllist), function(i){
    mobj <- mllist[[i]]
    parmMat <- mobj$parmMat
    if (ncol(parmMat) > 1){
      if (marginal == TRUE){
        lapply(1:length(respLev), function(rl) sapply(1:ncol(parmMat), function(j)  EDmarg(mobj, respLev=respLevMat[[rl]][j,i], clevel=colnames(parmMat)[j], display=FALSE, nGQ=nGQ, rfinterval=rfinterval)))
      } else {
        lapply(1:length(respLev), function(rl) sapply(1:ncol(parmMat), function(j)  ED(mobj, respLev=respLevMat[[rl]][j,i], clevel=colnames(parmMat)[j], display=FALSE)))
      }
    } else {
      if (marginal == TRUE){
        lapply(1:length(respLev), function(rl) t(EDmarg(mobj, respLev=respLevMat[[rl]][i], display=FALSE, nGQ=nGQ, rfinterval=rfinterval)))  
      } else {
        lapply(1:length(respLev), function(rl) t(ED(mobj, respLev=respLevMat[[rl]][i], display=FALSE)))  
      }
    }    
  })
  edm <- sapply(1:length(edl), function(i) sapply(edl[[i]], function(x) x[1,]))
  edse <- sapply(1:length(edl), function(i) sapply(edl[[i]], function(x) x[2,]))
  
  # model-averaged ED
  edma <- apply(t(t(edm) * wi), 1, sum)
  
  # rownames
  parmMat <- mllist[[1]]$parmMat
  if (ncol(parmMat) > 1){
    rcnames <- paste(colnames(parmMat), rep(respLev, each=ncol(parmMat)), sep=":")
  } else {
    rcnames <- as.character(respLev)
  }
  
  if (interval[1] == "none"){
    retMat <- cbind("Estimate" = edma)
  }  
  
  if (interval[1] == "buckland"){
    sema <- apply(sqrt(edse^2 + (edm - apply(edm, 1, mean))^2), 1, function(x) sum(x*wi))
    quant <- qnorm(1 - (1 - level)/2) * sema
    retMat <- as.matrix(cbind(edma, sema, edma - quant, edma + quant))
    colnames(retMat) <- c("Estimate", "SE", "Lower", "Upper")
  }
  
  if (interval[1] == "kang"){
    edlci <- lapply(1:length(mllist), function(i){
      mobj <- mllist[[i]]
      parmMat <- mobj$parmMat
      if (ncol(parmMat) > 1){
        if (marginal == TRUE){
          lapply(1:length(respLev), function(rl) sapply(1:ncol(parmMat), function(j)  EDmarg(mobj, respLev=respLevMat[[rl]][j,i], clevel=colnames(parmMat)[j], display=FALSE, interval="delta", nGQ=nGQ, rfinterval=rfinterval)))
        } else {
          lapply(1:length(respLev), function(rl) sapply(1:ncol(parmMat), function(j)  ED(mobj, respLev=respLevMat[[rl]][j,i], clevel=colnames(parmMat)[j], display=FALSE, interval="delta")))
        }
      } else {
        if (marginal == TRUE){
          lapply(1:length(respLev), function(rl) t(ED(mobj, respLev=respLevMat[[rl]][i], display=FALSE, interval="delta", nGQ=nGQ, rfinterval=rfinterval)))   
        } else {
          lapply(1:length(respLev), function(rl) t(ED(mobj, respLev=respLevMat[[rl]][i], display=FALSE, interval="delta")))   
        }
      }    
    })
    ll <- sapply(1:length(edlci), function(i) sapply(edlci[[i]], function(x) x[3,]))
    ul <- sapply(1:length(edlci), function(i) sapply(edlci[[i]], function(x) x[4,]))
    retMat <- as.matrix(cbind(edma, apply(ll, 1, function(x) sum(x*wi)), apply(ul, 1, function(x) sum(x*wi))))
    colnames(retMat) <- c("Estimate", "Lower", "Upper")
  }  
  
  rownames(retMat) <- rcnames
  rownames(edm) <- rcnames
  colnames(edm) <- mnames
  rlm <- sapply(1:length(mllist), function(i) sapply(respLevMat, function(x) rbind(x)[,i]))
  rownames(rlm) <- rcnames
  colnames(rlm) <- mnames
  
  out <- list()
  out$effectiveRespLev <- rlm
  out$retMat <- retMat
  out$weights <- wi
  out$EDi <- edm
  if (interval[1] == "buckland") {out$SEi <- edse}  
  return(out)  
  
}
