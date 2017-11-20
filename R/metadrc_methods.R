vcov.metadrc <- function(object, ...){
  return(object$vcov)  
}

df.residual.metadrc <- function(object, ...){
  df.resid <- object$rma$k.eff - object$rma$p.eff
  return(df.resid)
}

coef.metadrc <- function (object, ...){
  cdat <- data.frame(object$coefmat)
  names(cdat) <- object$parNames[[1]]
  return(cdat)
}

print.metadrc <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nTwo-stage meta-analysis dose-response model\nModel fitted: ", x$fct$text, "\n",  sep = "")
  
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}


summary.metadrc <- function(object, digits = max(3, getOption("digits") - 3), ...){
  cat("\nTwo-stage meta-analysis dose-response model\nModel fitted: ", object$fct$text, "\n",  sep = "")
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  x <- object$rma
  
  cat("Variance estimates:\n")
  mng <- max(nchar(x$g.names))
  right <- TRUE
  
  tau2 <- formatC(x$tau2, digits=digits, format="f")
  tau <- formatC(sqrt(x$tau2), digits=digits, format="f")
  rho <- formatC(x$rho, digits=digits, format="f")
  
  
  if (is.element(x$struct[1], c("CS","AR","ID"))) {
    
    vc <- cbind(tau2, tau)
    vc <- rbind(vc, c(rho, ""))
    colnames(vc) <- c("estim", "sqrt")
    rownames(vc) <- c("tau^2    ", "rho")
    if (x$struct[1] == "ID")
      vc <- vc[1,,drop=FALSE]
    print(vc, quote=FALSE, right=right, print.gap=2)
    
  }
  
  if (is.element(x$struct[1], c("HCS","HAR","DIAG"))) {
    
    vc <- cbind(tau2, tau)
    vc <- rbind(vc, c(rho, ""))
    colnames(vc) <- c("estim", "sqrt")
    if (length(x$tau2) == 1) {
      rownames(vc) <- c("tau^2   ", "rho")
    } else {
      rownames(vc) <- c(paste("tau^2.", seq_along(x$tau2), "  ", sep=""), "rho")
    }
    if (x$struct[1] == "DIAG")
      vc <- vc[seq_along(tau2),,drop=FALSE]
    print(vc, quote=FALSE, right=right, print.gap=2)
    
  }
  
  if (is.element(x$struct[1], c("UN","UNHO"))) {
    
    if (x$struct[1] == "UN") {
      vc <- cbind(tau2, tau)
    }
    colnames(vc) <- c("estim", "sqrt")
    if (length(x$g.levels.k) == 1) {
      rownames(vc) <- c("tau^2")
    } else {
      rownames(vc) <- paste("tau^2.", seq_along(x$g.levels.k), "  ", sep="")
    }
    print(vc, quote=FALSE, right=right, print.gap=2)
    cat("\n")
    
    if (length(x$rho) == 1) {
      G <- matrix(NA_real_, nrow=2, ncol=2)
    } else {
      G <- matrix(NA_real_, nrow=x$g.nlevels.f[1], ncol=x$g.nlevels.f[1])
    }
    G[upper.tri(G)] <- rho
    G[lower.tri(G)] <- t(G)[lower.tri(G)]
    diag(G) <- 1
    #G[upper.tri(G)] <- ""
    
    
    vc <- G
    colnames(vc) <- paste("rho.", abbreviate(x$g.levels.f[[1]]), sep="")
    rownames(vc) <- x$g.levels.f[[1]]
    print(vc, quote=FALSE, right=right, print.gap=2)
    
  }
  
  cat("\n")
  
  cat("\nCoefficients:\n")
  res.table <- cbind(estimate=c(x$beta), se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
  rownames(res.table) <- object$parNames[[1]]
  if (object$type == "continuous"){
    colnames(res.table)[c(1,2,3,4)] <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  } else {
    colnames(res.table)[c(1,2,3,4)] <- c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
  }
  printCoefmat(res.table[,c(1,2,3,4)], has.Pvalue = TRUE)
}

confint.metadrc <- function(object, parm, level=0.95, fixed=TRUE, random=FALSE, ...){
  confint.rma.mv(object$rma, fixed=fixed, random=random, ...)  
}