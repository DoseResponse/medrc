#' Asymptotic regression model
#' 
#' Providing the mean function and the corresponding self starter function for the asymptotic regression model.
#' 
#' The asymptotic regression model is a three-parameter model with mean function:
#' \deqn{ f(x) = c + (d-c)(1-\exp(-x/e))}
#' The parameter \eqn{c} is the lower limit (at \eqn{x=0}), the parameter \eqn{d} is the upper limit and the parameter \eqn{e>0} is determining the steepness of the increase as \eqn{x}.
#' 
#' @param x numeric vector specifying the dose
#' @param c lower limit at at \eqn{x=0}
#' @param d upper limit
#' @param e \eqn{e>0} is determining the steepness of the increase as \eqn{x}
#' 
#' @keywords models
#' 
#' @rdname meAR
meAR.2 <- function(x, d, e){
  .value <- d * (1 - exp(-exp(log(x) - log(e))))
  .actualArgs <- c("d", "e")
  t2 <- exp(log(x) - log(e))
  t3 <- exp(-t2)
  .grad <- cbind(1 - t3, 
                 -d * xexpx(x/e, 1)/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meAR
meAR.3 <- function(x, c, d, e){
  .value <- c + (d - c) * (1 - exp(-exp(log(x) - log(e))))
  .actualArgs <- c("c", "d", "e")
  t1 <- d - c
  t2 <- exp(log(x) - log(e))
  t3 <- exp(-t2)
  .grad <- cbind(1 - (1 - t3), 
                 1 - t3, 
                 -t1 * xexpx(x/e, 1)/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}