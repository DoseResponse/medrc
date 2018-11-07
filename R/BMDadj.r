#' Adjusted ED response levels for BMD estimation based on medrc or glsdrc models
#' 
#' Calculates adjusted response levels for estimation of the BMD
#' 
#' @param object an medrc object
#' @param respLev a numeric vector containing the benchmark response levels
#' @param bmd benchmark dose estimation (smallest dose resulting in a probability of an abnormal response)
#' @param background probability of an abnormal response
#' 
#' @keywords htest

BMDadjresp <- function (object, respLev, bmd = c("additional", "extra"), background = 0.05){ 
  bmd <- match.arg(bmd)
  if (bmd[1] == "extra") 
    respLev <- respLev * (1 - background)
  lenRL <- length(respLev)
  sapply(1:lenRL, function(i) {
    cobj <- object$parmMat
    prnames <- rownames(cobj)
    if (class(object)[1] == "drc") {
      prnames <- unique(object$parNames[[2]])
      rownames(cobj) <- prnames
    }
    cvals <- cobj[prnames == "c", ]
    if (length(cvals) == 0) {
      cvals <- 0
    }
    dvals <- cobj[prnames == "d", ]
    if (length(cvals) == 1 & length(dvals) > 1) 
      cvals <- rep(cvals, length(dvals))
    if (length(dvals) == 1 & length(cvals) > 1) 
      dvals <- rep(dvals, length(cvals))
    if (any(cvals > dvals)) {
      tempd <- apply(cbind(cvals, dvals), 1, function(x) sort(x))
      cvals <- tempd[1, ]
      dvals <- tempd[2, ]
    }
    if (class(object)[1] == "medrc") {
      varcorr <- VarCorr(object)
      return(100 * (qnorm(1 - background) - qnorm(1 - (background + respLev[i]/100))) * as.numeric(varcorr[attr(varcorr, "dimnames")[[1]] %in% "Residual", 2])/(dvals - cvals))
    }
    if (class(object)[1] == "glsdrc") {
      return(100 * (qnorm(1 - background) - qnorm(1 - (background + respLev[i]/100))) * object$fit$sigma/(dvals - cvals))
    }
    if (class(object)[1] == "drc") {
      return(100 * (qnorm(1 - background) - qnorm(1 - (background + respLev[i]/100))) * summary(object)$rseMat[1, 1]/(dvals - cvals))
    }
  })
}