#' Mixed effects dose-response curves
#' 
#' Implements drc nonlinear functions into the nlme framework for mixed effects dose-response modelling.
#' 
#' An application of medrm is shown on the help pages of data \code{\link{broccoli}}. EDx and selectivity indices can be calculated with functions \code{\link[drc]{ED}} and \code{\link[drc]{EDcomp}}. Model-averaged ED can be computed by function \code{\link{mmaED}}.
#' Marginal EDx and corresponding function are implemented in functions \code{\link{EDmarg}}, using numerical integration conditional on the etimated variance components.
#' 
#' @param form Formula describing the dose-response relationship
#' @param curveid Formula with parameter names on the left hand side (divided by +) and a column name in data, denoting a factor, to estimate separate parameters per factor-level. If NULL only fixed effects for a single curve will be estimated.
#' @param data a data.frame object
#' @param fct a function of the drc package
#' @param random a list or one-sided formula describing the random effects
#' @param correlation additional corClasses object
#' @param weights additional varClasses object
#' @param control list with nlme control arguments
#' @param start optional list with initial values for the fixed components. If NULL the initial values will be found automatically.
#' @param REML logical value. If TRUE the model is fit by maximizing the restricted log-likelihood, if FALSE log-likelihood is maximized.
#' 
#' @return An object of class medrc
#' 
#' @seealso \code{\link[drc]{drm}}, \code{\link[nlme]{nlme}}
#' 
#' @keywords models


medrm <-
function(form, curveid=NULL, data, fct, random, correlation=NULL, weights=NULL, control=NULL, start=NULL, REML=FALSE){
  callDetail <- match.call()
  fname <- fct$name
  ffixed <- fct$fixed

  mefct <- switch(fname, 
                  L.2 = "meL.2",
                  L.3 = "meL.3",
                  L.4 = "meL.4",
                  L.5 = "meL.5",
                  LL.2 = "meLL.2",
                  LL.3 = "meLL.3",
                  LL.4 = "meLL.4",
                  LL.5 = "meLL.5",
                  W1.2 = "meW1.2",
                  W1.3 = "meW1.3",
                  W1.4 = "meW1.4",
                  W2.2 = "meW2.2",
                  W2.3 = "meW2.3",
                  W2.4 = "meW2.4",
                  LN.2 = "meLN.2",
                  LN.3 = "meLN.3",
                  LN.4 = "meLN.4",
                  AR.3 = "meAR.3",
                  AR.2 = "meAR.2",
                  EXD.3 = "meEXD.3",
                  EXD.2 = "meEXD.2")
  
  mepnames <- switch(fname, 
                     L.2 = c("b", "e"),
                     L.3 = c("b", "d", "e"),
                     L.4 = c("b", "c", "d", "e"),
                     L.5 = c("b", "c", "d", "e", "f"),
                     LL.2 = c("b", "e"),
                     LL.3 = c("b", "d", "e"),
                     LL.4 = c("b", "c", "d", "e"),
                     LL.5 = c("b", "c", "d", "e", "f"),
                     W1.2 = c("b", "e"),
                     W1.3 = c("b", "d", "e"),
                     W1.4 = c("b", "c", "d", "e"),
                     W2.2 = c("b", "e"),
                     W2.3 = c("b", "d", "e"),
                     W2.4 = c("b", "c", "d", "e"),
                     LN.4 = c("b", "c", "d", "e"),
                     LN.3 = c("b", "d", "e"),
                     LN.2 = c("b", "e"),
                     AR.3 = c("c", "d", "e"),
                     AR.2 = c("d", "e"),
                     EXD.3 = c("c", "d", "e"),
                     EXD.2 = c("d", "e"))
  
  ppos <- switch(fname, 
                 L.2 = c(1, 4),
                 L.3 = c(1, 3, 4),
                 L.4 = 1:4,
                 L.5 = 1:5,
                 LL.2 = c(1, 4),
                 LL.3 = c(1, 3, 4),
                 LL.4 = 1:4,
                 LL.5 = 1:5,
                 W1.2 = c(1, 4),
                 W1.3 = c(1, 3, 4),
                 W1.4 = 1:4,
                 W2.2 = c(1, 4),
                 W2.3 = c(1, 3, 4),
                 W2.4 = 1:4,
                 LN.4 = 1:4,
                 LN.3 = c(1, 3, 4),
                 LN.2 = c(1, 4),
                 AR.3 = c(2, 3, 4),
                 AR.2 = c(3, 4),
                 EXD.3 = c(2, 3, 4),
                 EXD.2 = c(3, 4))
  
  pest <- mepnames[is.na(ffixed[ppos])]
  pfixed <- paste(mepnames[!is.na(ffixed[ppos])], "=", ffixed[ppos][!is.na(ffixed[ppos])], sep="")
  
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
  mform <- deparse(as.formula(paste(as.character(form)[2], as.character(form)[1], mefct, "(", as.character(form)[3], ", ", paste(pest, collapse=", "), if (all(pfixed != "=")) paste(",", paste(pfixed, collapse=", ")), ")", sep="")))  
  
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


