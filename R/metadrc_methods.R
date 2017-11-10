vcov.metadrc <- function(object, ...){
  return(object$vcov)  
}

df.residual.metadrc <- function(object, ...){
  df.resid <- object$rma$k.eff - object$rma$p.eff
  return(df.resid)
}

coef.metadrc <- function (object, ...){
  return(object$coefficients)
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