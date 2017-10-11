BMDmarg <- function (object, respLev, interval = c("none", "delta", "fls", "tfls"), clevel=NULL, level = ifelse(!(interval == "none"), 0.95, NULL), bmd = c("additional", "extra"), background = 0.05, nGQ=5, rfinterval=c(0, 1000), display = TRUE, ...){
  interval <- match.arg(interval)
  if (bmd[1] == "extra") respLev <- respLev * (1-background)  
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
  EDlistm <- function(parmChosen, respLev, intgrid=intgrid, intweights=w, rfinterval=rfinterval){
    parm <- object$fct$fixed
    parm[is.na(parm)] <- parmChosen
    parm2 <- object$fct$fixed
    
    varcorr <- VarCorr(object)
    vcr <- as.numeric(varcorr[attr(varcorr, "dimnames")[[1]] %in% "Residual", 2])
    adj1 <- 100 * (qnorm(1 - background) - qnorm(1-(background+respLev/100))) * vcr 
        
    cip <- if ("c" %in% colnames(intgrid)) intgrid[,"c"] else rep(0, nrow(intgrid))
    dip <- if ("d" %in% colnames(intgrid)) intgrid[,"d"] else rep(0, nrow(intgrid))
    mint <- function(d, parmChosen, intgrid, intweights, object, cip, dip, adj1){
      sapply(d, function(dx){
        sum(na.omit(intweights * sapply(1:nrow(intgrid), function(x){           
          pc <- if (is.na(parm2[2])) parmChosen["c"] else parm[2]
          pd <- if (is.na(parm2[3])) parmChosen["d"] else parm[3]
          p <- (adj1 / abs((pd + dip[x]) - (pc + cip[x])))
          if (parm[1] > 0) p <- 100 - p
          tval <- (pc + cip[x]) + ((pd + dip[x]) - (pc + cip[x])) * (p/100)
          object$fct$fct(dx, rbind(parmChosen + intgrid[x,])) - tval
        }) ))
      })  
    }
    myenv <- new.env()
    assign("object", object, envir = myenv)
    assign("parmChosen", parmChosen, envir = myenv) 
    assign("rfinterval", rfinterval, envir = myenv) 
    assign("intgrid", intgrid, envir = myenv) 
    assign("intweights", intweights, envir = myenv)   
    assign("cip", cip, envir = myenv)
    assign("dip", dip, envir = myenv)
    assign("adj1", adj1, envir = myenv)
    ede <- suppressWarnings(numericDeriv(quote(uniroot(mint, interval = rfinterval, parmChosen = parmChosen, intgrid = intgrid, intweights = intweights, object = object, cip=cip, dip=dip, adj1=adj1)$root), "parmChosen", myenv))
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
        EDeval <- EDlistm(parmChosen, respLev[j], intgrid=intgrid, intweights=w, rfinterval=rfinterval)
        EDval <- EDeval[[1]]
        dEDval <- EDeval[[2]]
        dEDmat[(i-1)*lenPV + j, parmInd] <- dEDval 
        oriMat[rowIndex, 1] <- EDval
        oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
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
    logBase <- exp(1)
    EDmat[, 1] <- exp(EDmat[, 1])
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