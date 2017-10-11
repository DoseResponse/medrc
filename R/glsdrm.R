glsdrm <-
function(form, curveid=NULL, data, fct, correlation=NULL, weights=NULL, control=NULL, start=NULL){
  callDetail <- match.call()
  # rewrite drc fct into nls selfstart function
  # drcfunction() is defined in global environment...
  makehelpfunction(fct)

  # one-sided curveid formula
  if (!is.null(curveid) & length(curveid) < 3){
    levelnames <- levels(data[, as.character(curveid)[2]])
    curveid <- NULL
  } else {
    levelnames <- NULL
  }
  
  # dose-response formula
  mform <- deparse(as.formula(paste(as.character(form)[2], as.character(form)[1], "drcfunction(", as.character(form)[3], ", ", paste(fct$names, collapse=", "), ")", sep="")))  
  
  # paste list for each curveid:parameter combination
  if (is.null(curveid)){
    mfixed <- deparse(as.formula(paste(paste(fct$names, collapse=" + "), "~ 1")))
  } else {
    stsp <- strsplit(as.character(curveid)[2], "+")[[1]]
    cin <- stsp[!stsp %in% c(" ", "+")]
    pnl <- sapply(fct$names, function(x) if (x %in% cin) paste(x, "~", as.character(curveid)[3], "-1") else paste(x, "~ 1"))
    mfixed <- paste("list(", paste(pnl, collapse=", "), ")", sep="")
  }
  
  # starting values for fixed effects (drc fct selfstart slot)
  if (is.null(start)){       
    ini <- findfixedstartvalsglsdrm(form, data, as.character(curveid)[3], fct, fid=!fct$names %in% cin, mform)
  } else {
    ini <- start
  }     

  # gnls call
  if (is.null(ini)){
    gnlscall <- paste("gnls(", mform, ", params=", mfixed, ", data=data, correlation=correlation, weights=weights, control=control)", sep="")
  } else {
    gnlscall <- paste("gnls(", mform, ", params=", mfixed, ", data=data, start=ini, correlation=correlation, weights=weights, control=control)", sep="")
  }
  fmgls <- eval(parse(text=gnlscall))
  
  ### output
  out <- list()
  out$call <- callDetail
  out$form <- form
  out$curveid <- curveid
  out$fit <- fmgls
  out$fct <- fct
  out$data <- data    
  if (is.null(curveid)){  
    if (is.null(levelnames)){
      out$coefficients <- coefficients(fmgls)
      out$parmMat <- cbind(coefficients(fmgls))
      colnames(out$parmMat) <- "1"
      out$vc <- vcov(fmgls)
      out$indexMat <- cbind(1:length(coefficients(fmgls)))
    } else {
      pnames <- fct$names      
      ocoef <- coefficients(fmgls)
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
      out$vc <- indmat %*% vcov(fmgls) %*% t(indmat)
      out$indexMat <- matrix(seq(1, flev*length(pnames), by=1), ncol=flev, byrow=TRUE)   
    }} else {
    pnames <- fct$names
    stsp <- strsplit(as.character(curveid)[2], "+")[[1]]
    cin <- stsp[!stsp %in% c(" ", "+")]
    cp <- pnames %in% cin
    lev <- levels(data[,as.character(curveid)[3]])
    flev <- length(levels(data[,as.character(curveid)[3]]))        
    imat <- t(ldply(lapply(cp, function(i) if (i == TRUE) diag(flev) else rep(1, flev)), rbind))
    out$imat <- imat
    cf <- apply(imat, 1, function(ind) coefficients(fmgls)[as.logical(ind)])
    colnames(cf) <- lev
    rownames(cf) <- pnames
    out$parmMat <- cf
    ii <- as.vector(t(apply(imat, 1, function(x) (1:ncol(imat))[as.logical(x)])))
    indmat <- t(sapply(ii, function(i){
      x <- numeric(length=ncol(imat))
      x[i] <- 1
      return(x)
    }))
    out$vc <- indmat %*% vcov(fmgls) %*% t(indmat)
    cvec <- as.vector(indmat %*% coefficients(fmgls))
    cnames <- paste(rep(pnames, each=flev), lev, sep=".")
    names(cvec) <- cnames
    out$coefficients <- cvec
    out$indexMat <- matrix(seq(1, flev*length(pnames), by=1), ncol=flev, byrow=TRUE)    
  }  
  out$type = "continuous"
  pll <- attr(logLik(fmgls), "df")
  nll <- nobs(logLik(fmgls))
  out$mselect <- c("AIC"=AIC(fmgls), "BIC"=BIC(fmgls), "logLik"=logLik(fmgls), "df"=pll)
  out$start <- ini
  class(out) <- c("glsdrc", "drc")
  return(out)
}


