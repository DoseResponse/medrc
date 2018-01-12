#' Estimating marginal effective doses for medrc objects
#' 
#' \code{ED} estimates marginal effective doses (ECp/EDp/ICp) for given reponse levels, conditional on the estimated variance components.
#' 
#' @param object an object of class 'medrc'.
#' @param respLev a numeric vector containing the response levels.
#' @param interval character string specifying the type of confidence intervals to be supplied. The default is "none". Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution). Use "fls" for from logarithm scale based confidence intervals (in case the parameter in the model is log(ED50) as for the \code{\link{llogistic2}}) models. The only alternative for model-robust fits is using inverse regression.
#' @param clevel character string specifying the curve id in case on estimates for a specific curve or compound is requested. By default estimates are shown for all curves.
#' @param level numeric. The level for the confidence intervals. The default is 0.95.
#' @param reference character string. Is the upper limit or the control level the reference?
#' @param type character string. Whether the specified response levels are absolute or relative (default).
#' @param nGQ integer. Specifies the number nof nodes for Gauss-Hermite quadrature.
#' @param rfinterval numeric vector. Interval for root finding (uniroot) to search for ED values. 
#' @param lref numeric value specifying the lower limit to serve reference.
#' @param uref numeric value specifying the upper limit to serve reference (eg. 100\%).
#' @param bound logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful for hormesis models.
#' @param display logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).
#' @param logBase numeric. The base of the logarithm in case logarithm transformed dose values are used.
#' @param ... additional arguments
#' 
#' @return A matrix with two or more columns, containing the estimates and the corresponding estimated standard errors and possibly lower and upper confidence limits.
#' 
#' @keywords htest

EDmarg <- function (object, respLev, interval = c("none", "delta", "fls", "tfls"), clevel=NULL, level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), type = c("relative", "absolute"), nGQ=5, rfinterval=c(0, 1000), lref, uref, bound = TRUE, display = TRUE, logBase = NULL, ...){
  interval <- match.arg(interval)
  reference <- match.arg(reference)
  type <- match.arg(type)
  if ((type == "relative") && (bound)) {
    if (any(respLev <= 0 | respLev >= 100)) {
      stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
    }
  }
  EDlist <- object$fct$edfct
  if (is.null(EDlist)) {
    stop("ED values cannot be calculated")
  }
  indexMat <- object$indexMat
  parmMat <- object$parmMat
  
  # construct integration grid
  pnames <- rownames(parmMat)
  
  # Pinheiro function
  varRan <- function(object, level = 1){
    sigE <- object$sig^2
    sigE*pdMatrix(object$modelStruct$reStruct)[[level]]
  }

  # extract variance components
  ranefs <- ranef(object$fit)
  nvc <- length(object$fit$modelStruct$reStruct)
  vclist <- lapply(1:nvc, function(i) varRan(object$fit, level=i))
  strspl <- lapply(vclist, function(x) strsplit(names(diag(x)), ".", fixed=TRUE))
  nrnl <- lapply(strspl, function(x) sapply(x, function(x) x[1]))
  
  ###
  gq <- gauss.quad.prob(n=nGQ, dist="normal", sigma=1)  
  
  igwlist <- lapply(1:nvc, function(v){
    nrn <- nrnl[[v]]
    nv <- length(diag(vclist[[v]]))
    
    eg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$nodes", nv), collapse=","), ")", sep="")))
    weg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$weights", nv), collapse=","), ")", sep="")))
    w <- apply(weg, 1, function(x) prod(x))
    
    cfvv <- chol(vclist[[v]])
    z <- as.matrix(eg) %*% cfvv
    
    intgrid <- matrix(0, ncol=length(pnames), nrow=nrow(z))
    wn <- sapply(1:length(nrn), function(i) which(pnames == nrn[i]))
    intgrid[,wn] <- as.matrix(z)
    colnames(intgrid) <- pnames
    return(list(intgrid, w))
  })
  ###   
  iglist <- lapply(igwlist, function(x) x[[1]])
  ni <- sapply(iglist, nrow)
  vllist <- lapply(1:length(ni), function(i) apply(iglist[[i]], 2, function(ir){
    if (i == 1) return(rep(ir, each=prod(ni[(i+1):length(ni)])))
    if (i > 1 & i < length(ni)) return(rep(rep(ir, each=prod(ni[(i+1):length(ni)])), times=prod(ni[1:(i-1)])))
    if (i == length(ni)) return(rep(ir, times=prod(ni[1:(i-1)])))    
  }))
  arr <- simplify2array(vllist)
  intgrid <- apply(arr, c(1,2), sum)
  
  wlist <- lapply(igwlist, function(x) x[[2]])
  wmat <- sapply(1:length(ni), function(i){
    if (i == 1) return(rep(wlist[[i]], each=prod(ni[(i+1):length(ni)])))
    if (i > 1 & i < length(ni)) return(rep(rep(wlist[[i]], each=prod(ni[(i+1):length(ni)])), times=prod(ni[1:(i-1)])))
    if (i == length(ni)) return(rep(wlist[[i]], times=prod(ni[1:(i-1)])))    
  })
  w <- apply(wmat, 1, prod)
  
  strParm0 <- sort(colnames(parmMat))
  curveNames <- colnames(parmMat)
  options(warn = -1)
  if (any(is.na(as.numeric(curveNames)))) {
    curveOrder <- order(curveNames)
  } else {
    curveOrder <- 1:length(curveNames)
  }
  options(warn = 0)
  strParm0 <- curveNames[curveOrder]
  indexMat <- indexMat[, curveOrder, drop = FALSE]
  parmMat <- parmMat[, curveOrder, drop = FALSE]
  strParm <- strParm0
  vcMat <- vcov(object)
  ncolIM <- ncol(indexMat)
  indexVec <- 1:ncolIM
  lenPV <- length(respLev)
  noRows <- ncolIM * lenPV
  dimNames <- rep("", noRows)
  EDmat <- matrix(0, noRows, 2)
  oriMat <- matrix(0, noRows, 2)
  if (identical(length(unique(strParm)), 1)) {
    strParm[indexVec] <- rep("", ncolIM)
  } else {
    strParm <- paste(strParm, ":", sep = "")
  }
  
  # ED estimation by root finding based on marginal prediction
  EDlistm <- function(parmChosen, respLev, reference=reference, type=type, intgrid=intgrid, intweights=w, rfinterval=rfinterval){
    parm <- object$fct$fixed
    parm[is.na(parm)] <- parmChosen
    parm2 <- object$fct$fixed
    p <- 100-eval(parse(text="drc:::EDhelper(parm, respLev, reference = reference, type = type)"))
    cip <- if ("c" %in% colnames(intgrid)) intgrid[,"c"] else rep(0, nrow(intgrid))
    dip <- if ("d" %in% colnames(intgrid)) intgrid[,"d"] else rep(0, nrow(intgrid))
    mint <- function(d, parmChosen, intgrid, intweights, object, cip, dip) sapply(d, function(dx) sum(na.omit(intweights * sapply(1:nrow(intgrid), function(x){
      pc <- if (is.na(parm2[2])) parmChosen["c"] else parm[2]
      pd <- if (is.na(parm2[3])) parmChosen["d"] else parm[3]
      object$fct$fct(dx, rbind(parmChosen + intgrid[x,])) - ((pc + cip[x]) + ((pd + dip[x]) - (pc + cip[x])) * (p/100)) 
    }))))    
    myenv <- new.env()
    assign("object", object, envir = myenv)
    assign("parmChosen", parmChosen, envir = myenv) 
    assign("rfinterval", rfinterval, envir = myenv) 
    assign("intgrid", intgrid, envir = myenv) 
    assign("intweights", intweights, envir = myenv)   
    assign("cip", cip, envir = myenv)
    assign("dip", dip, envir = myenv)
    ede <- suppressWarnings(numericDeriv(quote(uniroot(mint, interval = rfinterval, parmChosen = parmChosen, intgrid = intgrid, intweights = intweights, object = object, cip=cip, dip=dip)$root), "parmChosen", myenv))
    out <- list()
    out[[1]] <- ede
    out[[2]] <- as.vector(attr(ede, "gradient"))
    return(out)
  }
  
  lenIV <- length(indexVec)
  dEDmat <- matrix(0, lenPV * lenIV, nrow(vcMat))
  rowIndex <- 1
  for (i in indexVec) {
    parmChosen <- parmMat[, i]
    parmInd <- indexMat[, i]
    varCov <- vcMat[parmInd, parmInd]
    if ((is.null(clevel)) || (strParm0[i] %in% clevel)){
      for (j in 1:lenPV) {
        EDeval <- EDlistm(parmChosen, respLev[j], reference = reference, type = type, intgrid=intgrid, intweights=w, rfinterval=rfinterval)
        EDval <- EDeval[[1]]
        dEDval <- EDeval[[2]]
        dEDmat[(i-1)*lenPV + j, parmInd] <- dEDval 
        oriMat[rowIndex, 1] <- EDval
        oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
        if (!is.null(logBase)) {
          EDval <- logBase^(EDval)
          dEDval <- EDval * log(logBase) * dEDval
        }
        EDmat[rowIndex, 1] <- EDval
        EDmat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
        dimNames[rowIndex] <- paste(strParm[i], respLev[j], sep = "")
        rowIndex <- rowIndex + 1
      }    
    } else {
      rowsToRemove <- rowIndex:(rowIndex + lenPV - 1)
      EDmat <- EDmat[-rowsToRemove, , drop = FALSE]
      dimNames <- dimNames[-rowsToRemove]
    }
  } 
  
  colNames <- c("Estimate", "Std. Error")
  if (interval == "delta") {
    intMat <- eval(parse(text="drc:::confint.basic(EDmat, level, object$type, df.residual(object), FALSE)"))
    intLabel <- "Delta method"
  }
  if (interval == "tfls") {
    intMat <- eval(parse(text="exp(drc:::confint.basic(matrix(c(log(oriMat[, 1]), oriMat[, 2]/oriMat[, 1]), ncol = 2), level, object$type, df.residual(object), FALSE))"))
    intLabel <- "To and from log scale"
  }
  if (interval == "fls") {
    if (is.null(logBase)) {
      logBase <- exp(1)
      EDmat[, 1] <- exp(EDmat[, 1])
    }
    intMat <- eval(parse(text="logBase^(drc:::confint.basic(oriMat, level, object$type, df.residual(object), FALSE))"))
    intLabel <- "Back-transformed from log scale"
    EDmat <- EDmat[, -2, drop = FALSE]
    colNames <- colNames[-2]
  }
  if (identical(interval, "none")) {
    intLabel <- NULL
  } else {
    EDmat <- as.matrix(cbind(EDmat, intMat))
    colNames <- c(colNames, "Lower", "Upper")
  }
  dimnames(EDmat) <- list(dimNames, colNames)
  eval(parse(text="drc:::resPrint(EDmat, 'Estimated effective doses', interval, intLabel, display = display)"))
  invisible(list(EDmat, EDmultcomp = parm(EDmat[, 1], (dEDmat %*% vcMat %*% t(dEDmat))[1:nrow(EDmat), 1:nrow(EDmat),drop=FALSE])))
}