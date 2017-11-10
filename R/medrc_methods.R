print.medrc <- function (x, ..., digits = max(3, getOption("digits") - 3)){
  object <- x
  classList <- class(object)
  cat(paste("\n", "A 'medrc' model.", "\n", sep = ""))
  cat("\nCall:\n", paste(deparse(object$call), collapse="\n"), "\n\n", sep = "")
  if (length(coef(object)) > 0) {
    cat("Coefficients:\n")
    print.default(format(coef(object), digits = digits), print.gap = 2, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(object)
}




residuals.medrc <- function(object, ...){
  resid(object$fit, type="n")
}

predict.medrc <- function(object, ..., newdata=NULL, level=NULL, type=c("conditional", "marginal"), nGQ=5){
  fct <- object$fct
  if (is.null(newdata)) newdata <- object$data
  if (is.null(level)) level <- object$fit$dims$Q
  if (type[1] == "conditional"){
    return(predict(object$fit, newdata=newdata, level=level))
  }
  if (type[1] == "marginal"){
    ## curveid parmMat matching...
    parmMat <- object$parmMat
    if (ncol(parmMat) == 1){
      idi <- rep(1, nrow(newdata))
    } else {
      nci <- as.character(object$curveid[3])
      if (!nci %in% colnames(newdata)) stop(paste("Cannot find the factor ", nci, " in newdata!", sep=""))
      cl <- colnames(parmMat)
      if (!all(unique(newdata[,nci]) %in% cl)) stop(paste("Cannot find the level(s) ", unique(newdata[,nci])[!(unique(newdata[,nci]) %in% cl)], " of factor ", nci, " in newdata!", sep=""))
      idi <- sapply(newdata[,nci], function(x) which(x == cl))       
    }
    
    dname <- as.character(object$form[3])
    dose <- newdata[,dname]
    
    # integration grid conditional on estimated variance components
    pnames <- rownames(parmMat)
    vc <- VarCorr(object$fit)
    nv <- nrow(vc)
    if (nv == 2){
      nvc <- names(ranef(object$fit))
    } else{
      nvc <- rownames(vc)[-nv] 
    }
    ###
    # single vc per parameter
    rnspl <- strsplit(nvc, ".", fixed=TRUE)
    nrn <- sapply(rnspl, function(x) x[1])
    rest <- sapply(rnspl, function(x) x[2])
    
    std <- as.numeric(vc[-nv,2])
    cormat <- diag(length(std))
    if (ncol(vc) > 2){    
      cormat[upper.tri(cormat)] <- cormat[lower.tri(cormat)] <- na.omit(as.numeric(vc[-c(1,nv),-c(1,2)]))
    }
    gq <- gauss.quad.prob(n=nGQ, dist="normal", sigma=1)
    
    eg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$nodes", nv-1), collapse=","), ")", sep="")))
    weg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$weights", nv-1), collapse=","), ")", sep="")))
    w <- apply(weg, 1, function(x) prod(x))
    
    ee <- eigen(cormat)
    A <- ee$vectors %*% diag(sqrt(ee$values))
    z <- data.frame(t(std*(A %*% t(as.matrix(eg)))))
    
    intgrid <- matrix(0, ncol=length(pnames), nrow=nrow(z))
    wn <- sapply(1:length(nrn), function(i) which(pnames == nrn[i]))
    intgrid[,wn] <- as.matrix(z)
    colnames(intgrid) <- pnames
    mpred <- function(d, parmChosen, intgrid, w) sum(na.omit(w*apply(intgrid, 1, function(x) object$fct$fct(d, rbind(parmChosen + x)))))
    
    
    # marginal predictions by Gauss-Hermite quadrature
    mp <- sapply(1:nrow(newdata), function(i){
      mpred(dose[i], parmMat[,idi[i]], intgrid, w)
    })
    return(mp)                 
  }  
}


ranef.medrc <- function(object, ...){
  ranef(object$fit, ...)
}

VarCorr.medrc <- function(x, sigma = 1, ...){
  VarCorr(x$fit, sigma=sigma)
}

vcov.medrc <- function(object, ...){
  #est <- fixef(object$fit)
  #stdFixed <- sqrt(diag(as.matrix(object$fit$varFix)))
  #std <- sqrt(object$fit$dims$N/(object$fit$dims$N - length(stdFixed))) * stdFixed
  #cr <- array(t(object$fit$varFix/stdFixed)/stdFixed, dim(object$fit$varFix), list(names(est), names(est)))
  #return(cor2cov(round(cr,10), std))
  object$vc
}

summary.medrc <- function(object, ...){
  summary(object$fit)
}


df.residual.medrc <- function(object, ...){
  # need to define a better residual.df !!!!!
  object$fit$dims$N - length(coefficients(object)) - nrow(VarCorr(object$fit))
}


AIC.medrc <- function(object, ..., k = 2){
  objlist <- list(object, ...)
  if (length(objlist) > 1){
    ismedrc <- sapply(objlist, function(x) inherits(x, "medrc"))
    medrclist <- objlist[ismedrc]  
    Call <- match.call()
    Call$k <- NULL
    names(medrclist) <- as.character(Call[-1L])[ismedrc]
    nlmelist <- lapply(medrclist, function(x) x$fit)
    ftext <- paste("AIC(", paste(paste("nlmelist[['", names(nlmelist), "']]", sep=""), collapse=","), ", k=", k, ")", sep="")
    AICtab <- eval(parse(text=ftext))
    rownames(AICtab) <- names(nlmelist)
  } else {
    AICtab <- AIC(object$fit, k=k)
  }
  return(AICtab)
}

BIC.medrc <- function(object, ...){
  objlist <- list(object, ...)
  if (length(objlist) > 1){
    ismedrc <- sapply(objlist, function(x) inherits(x, "medrc"))
    medrclist <- objlist[ismedrc]  
    Call <- match.call()
    Call$k <- NULL
    names(medrclist) <- as.character(Call[-1L])[ismedrc]
    nlmelist <- lapply(medrclist, function(x) x$fit)
    ftext <- paste("BIC(", paste(paste("nlmelist[['", names(nlmelist), "']]", sep=""), collapse=","), ")", sep="")
    BICtab <- eval(parse(text=ftext))
    rownames(BICtab) <- names(nlmelist)
  } else {
    BICtab <- BIC(object$fit)
  }
  return(BICtab)
}

logLik.medrc <- function(object, REML = FALSE, ...){
  logLik(object$fit, REML=REML, ...) 
}



