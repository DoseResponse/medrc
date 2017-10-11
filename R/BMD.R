BMD <- function(object, respLev, interval = c("none", "delta", "fls", "tfls"), level = 0.95, bmd = c("additional", "extra"), background = 0.05, display=TRUE){
  interval <- match.arg(interval)
  bmd <- match.arg(bmd)
  # adjusting respLev for BMD
  if (bmd[1] == "extra") respLev <- respLev * (1-background)  
  lenRL <- length(respLev)
  adjrespLev <- sapply(1:lenRL, function(i){ 
    cobj <- object$parmMat
    prnames <- rownames(cobj)
    if (class(object)[1] == "drc"){
      prnames <- unique(object$parNames[[2]])
      rownames(cobj) <- prnames
    }
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
    if (class(object)[1] == "medrc"){
      varcorr <- VarCorr(object)
      return(100 * (qnorm(1 - background) - qnorm(1-(background+respLev[i]/100))) * as.numeric(varcorr[attr(varcorr, "dimnames")[[1]] %in% "Residual", 2]) / (dvals - cvals))
    }
    if (class(object)[1] == "glsdrc"){
      return(100 * (qnorm(1 - background) - qnorm(1-(background+respLev[i]/100))) * object$fit$sigma / (dvals - cvals))     
    }
    if (class(object)[1] == "drc"){
      return(100 * (qnorm(1 - background) - qnorm(1-(background+respLev[i]/100))) * summary(object)$rseMat[1,1] / (dvals - cvals))     
    }
  }) 
  intLabel <- NULL
  if (interval == "delta") intLabel <- "Delta method"
  if (interval == "tfls") intLabel <- "To and from log scale"
  if (interval == "fls") intLabel <- "Back-transformed from log scale"
  if (ncol(object$parmMat) > 1){ 
    cid <- colnames(object$parmMat)
    bmdrn <- rownames(ED(object, respLev=adjrespLev, display=FALSE)[[1]])
    el <- lapply(1:length(cid), function(ci) ED(object, respLev=adjrespLev[ci,], clevel=cid[ci], interval=interval, level=level, display=FALSE)[[1]])
    bmdest <- matrix(simplify2array(el), ncol=ncol(el[[1]]))
    rownames(bmdest) <- as.vector(sapply(el, rownames))
    colnames(bmdest) <- colnames(el[[1]])    
  } else {
    bmdrn <- rownames(ED(object, respLev=respLev, display=FALSE)[[1]])
    bmdest <- ED(object, respLev=adjrespLev, interval=interval, level=level, display=FALSE)[[1]]
    rownames(bmdest) <- bmdrn
  }
  eval(parse(text="drc:::resPrint(bmdest, 'Estimated benchmark doses', interval, intLabel, display = display)"))
}