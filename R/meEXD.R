#' Exponential decay model
#' 
#' Exponential decay model with or without a nonzero lower limit.
#' 
#' The exponential decay model is a three-parameter model with mean function:
#' \deqn{f(x) = c + (d-c)(\exp(-x/e))}
#' The parameter init is the upper limit (attained at \eqn{x=0}), the parameter plateau is the lower limit reached for x going to infinity and the parameter \eqn{e>0} is determining the steepness of the decay. The curve is monotonously decreasing in \eqn{x}.
#' 
#' @param x numeric dose vector
#' @param c lower limit
#' @param d upper limit
#' @param e steepness
#' 
#' @references Organisation for Economic Co-operation and Development (OECD) (2006) \emph{Current approaches in the statistical analysis of ecotoxicity data: A guidance to application - annexes}, Paris: OECD (p. 80).
#' 
#' @keywords models
#' 
#' @rdname meEXD
meEXD.2 <- function(x, d, e){
  .value <- d * exp(-exp(log(x) - log(e)))
  .actualArgs <- c("d", "e")
  t2 <- exp(log(x) - log(e))
  t3 <- exp(-t2)
  .grad <- cbind(t3, 
                 d * divAtInf(t2, exp(t2))/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meEXD
meEXD.3 <- function(x, c, d, e){
  .value <- c + (d - c) * exp(-exp(log(x) - log(e)))
  .actualArgs <- c("c", "d", "e")
  t1 <- d - c
  t2 <- exp(log(x) - log(e))
  t3 <- exp(-t2)
  .grad <- cbind(1 - t3, 
                 t3, 
                 t1 * divAtInf(t2, exp(t2))/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}
