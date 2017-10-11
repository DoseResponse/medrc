medrm <-
function(form, curveid=NULL, data, fct, random, correlation=NULL, weights=NULL, control=NULL, start=NULL, REML=FALSE){
  callDetail <- match.call()
  # rewrite drc fct into nls selfstart function
  # drcfunction() is defined in global environment...
  makehelpfunction(fct)

  # one-sided curveid formula
  if (!is.null(curveid) & length(curveid) < 3){
    levelnames <- levels(data[, as.character(curveid)[2]])
    cid1s <- curveid
    curveid <- NULL    
  } else {
    levelnames <- NULL
    cid1s <- NULL
  }  
  
  # dose-response formula
  mform <- deparse(as.formula(paste(as.character(form)[2], as.character(form)[1], "drcfunction(", as.character(form)[3], ", ", paste(fct$names, collapse=", "), ")", sep="")))  
  
  # fixed effects formula
  # paste list for each curveid:parameter combination
  if (is.null(curveid)){
    mfixed <- deparse(as.formula(paste(paste(fct$names, collapse=" + "), "~ 1")))
  } else {
    cin <- all.vars(curveid[[2]])
    pnl <- sapply(fct$names, function(x) if (x %in% cin) paste(x, "~", as.character(curveid)[3], "-1") else paste(x, "~ 1"))
    mfixed <- paste("list(", paste(pnl, collapse=", "), ")", sep="")
  }
  
  # starting values for fixed effects (drc fct selfstart slot)
  if (is.null(start)){       
    ini <- findfixedstartvals(form, data, as.character(curveid)[3], fct, fid=!fct$names %in% cin, mform)
#     if (randomstart == TRUE){
#       rini <- findrandomstartvals(form, data, fct, random)
#       ini <- list(fixed=ini, random=rini)
#     }
  } else {
    ini <- start
  }     

  # nlme call
  if (REML == TRUE) remlmethod <- "'REML'" else remlmethod <- "'ML'"
  nlmecall <- paste("nlme(", mform, ", fixed=", mfixed, ", data=data, random=", deparse(callDetail$random), ", start=ini, correlation=correlation, weights=weights, control=control, method=", remlmethod, ")", sep="")
  fmmixed <- eval(parse(text=nlmecall))
  
  
  ### output
  out <- list()
  out$call <- callDetail
  out$form <- form
  out$curveid <- curveid
  out$rform <- random
  out$fit <- fmmixed
  out$fct <- fct
  out$data <- data
  if (is.null(curveid)){ 
    if (is.null(levelnames)){
      out$coefficients <- fixef(fmmixed)
      out$parmMat <- cbind(fixef(fmmixed))
      colnames(out$parmMat) <- "1"
      out$vc <- vcov(fmmixed)
      out$indexMat <- cbind(1:length(fixef(fmmixed)))
    } else {
      pnames <- fct$names      
      ocoef <- fixef(fmmixed)
      flev <- length(levelnames)
      rcoefs <- rep(ocoef, each=flev)
      cnames <- paste(rep(pnames, each=flev), levelnames, sep=".")
      names(rcoefs) <- cnames
      out$coefficients <- rcoefs      
      cf <- matrix(rep(ocoef, each=flev), ncol=flev, byrow=TRUE)
      colnames(cf) <- levelnames
      rownames(cf) <- pnames
      out$parmMat <- cf
      imat <- matrix(1, nrow=flev, ncol=length(ocoef))
      out$imat <- imat
      ii <- as.vector(t(apply(imat, 1, function(x) (1:ncol(imat))[as.logical(x)])))
      indmat <- t(sapply(ii, function(i){
        x <- numeric(length=ncol(imat))
        x[i] <- 1
        return(x)
      }))
      out$vc <- indmat %*% vcov(fmmixed) %*% t(indmat)
      out$indexMat <- matrix(seq(1, flev*length(pnames), by=1), ncol=flev, byrow=TRUE)         
    }} else {
    pnames <- fct$names
    cin <- all.vars(curveid[[2]])
    cp <- pnames %in% cin
    lev <- levels(data[,as.character(curveid)[3]])
    flev <- length(levels(data[,as.character(curveid)[3]]))        
    imat <- t(ldply(lapply(cp, function(i) if (i == TRUE) diag(flev) else rep(1, flev)), rbind))
    out$imat <- imat
    cf <- apply(imat, 1, function(ind) fixef(fmmixed)[as.logical(ind)])
    colnames(cf) <- lev
    rownames(cf) <- pnames
    out$parmMat <- cf
    ii <- as.vector(t(apply(imat, 1, function(x) (1:ncol(imat))[as.logical(x)])))
    indmat <- t(sapply(ii, function(i){
      x <- numeric(length=ncol(imat))
      x[i] <- 1
      return(x)
    }))
    out$vc <- indmat %*% vcov(fmmixed) %*% t(indmat)
    cvec <- as.vector(indmat %*% fixef(fmmixed))
    cnames <- paste(rep(pnames, each=flev), lev, sep=".")
    names(cvec) <- cnames
    out$coefficients <- cvec
    out$indexMat <- matrix(seq(1, flev*length(pnames), by=1), ncol=flev, byrow=TRUE)    
  }
  out$type = "continuous"
  pll <- attr(logLik(fmmixed), "df")
  nll <- nobs(logLik(fmmixed))
  out$mselect <- c("AIC"=AIC(fmmixed), "AICc"=AIC(fmmixed) + ((2*pll*(pll+1))/(nll-pll-1)), "BIC"=BIC(fmmixed), "logLik"=logLik(fmmixed), "df"=pll)
  out$start <- ini
  if (!is.null(cid1s)) out$curveid <- cid1s
  class(out) <- c("medrc", "drc")
  return(out)
}


