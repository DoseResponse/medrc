#' Four-parameter Weibull functions
#' 
#' 'W1.X' and 'W2.X' provide the X-parameter Weibull functions.
#' 
#' As pointed out in Seber and Wild (1989), there exist two different parameterisations of the Weibull model. They do not yield the same fitted curve for a given dataset.
#' The four-parameter Weibull type 1 model ('weibull1') is
#' \deqn{ f(x) = c + (d-c) \exp(-\exp(b(\log(x)-\log(e)))).}
#' Thw four-parameter Weibull type 2 model ('weibull2') is
#' \deqn{ f(x) = c + (d-c) (1 - \exp(-\exp(b(\log(x)-\log(e))))).} 
#' Both four-parameter model functions are asymmetric with inflection point at the dose equal \eqn{e}.
#' 
#' @param x numeric dose vector.
#' @param b steepness
#' @param c lower limit
#' @param d upper limit
#' @param e ED50
#' 
#' @references 
#' Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 330--331).
#' Ritz, C (2009) Towards a unified approach to dose-response modeling in ecotoxicology \emph{To appear in Environ Toxicol Chem}.
#' 
#' @keywords models
#' 
#' @rdname meW
meW1.2 <- function(x, b, e){
  .value <- exp(-exp(b * (log(x) - log(e))))
  .actualArgs <- c("b", "e")
  t2 <- exp(b * (log(x) - log(e)))
  .grad <- cbind(-1 * divAtInf(xlogx(x/e, b), exp(t2)), 
                 divAtInf(t2, exp(t2)) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meW
meW1.3 <- function(x, b, d, e){
  .value <- d * exp(-exp(b * (log(x) - log(e))))
  .actualArgs <- c("b", "d", "e")
  t2 <- exp(b * (log(x) - log(e)))
  t3 <- exp(-t2)
  .grad <- cbind(-d * divAtInf(xlogx(x/e, b), exp(t2)), 
                 t3, 
                 d * divAtInf(t2, exp(t2)) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meW
meW1.4 <- function(x, b, c, d, e){
  .value <- c + (d - c) * exp(-exp(b * (log(x) - log(e))))
  .actualArgs <- c("b", "c", "d", "e")
  t1 <- d - c
  t2 <- exp(b * (log(x) - log(e)))
  t3 <- exp(-t2)
  .grad <- cbind(-t1 * divAtInf(xlogx(x/e, b), exp(t2)), 
                 1 - t3, 
                 t3, 
                 t1 * divAtInf(t2, exp(t2)) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meW
meW2.2 <- function(x, b, e){
  .value <- (1 - exp(-exp(b * (log(x) - log(e)))))
  .actualArgs <- c("b", "e")
  t2 <- exp(b * (log(x) - log(e)))
  t3 <- exp(-t2)
  .grad <- cbind(xexplogx(x/e, b), 
                 -1 * xexpx(x/e, b) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meW
meW2.3 <- function(x, b, d, e){
  .value <- d * (1 - exp(-exp(b * (log(x) - log(e)))))
  .actualArgs <- c("b", "d", "e")
  t2 <- exp(b * (log(x) - log(e)))
  t3 <- exp(-t2)
  .grad <- cbind(d * xexplogx(x/e, b), 
                 1 - t3, 
                 -d * xexpx(x/e, b) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meW
meW2.4 <- function(x, b, c, d, e){
  .value <- c + (d - c) * (1 - exp(-exp(b * (log(x) - log(e)))))
  .actualArgs <- c("b", "c", "d", "e")
  t1 <- d - c
  t2 <- exp(b * (log(x) - log(e)))
  t3 <- exp(-t2)
  .grad <- cbind(t1 * xexplogx(x/e, b), 
                 1 - (1 - t3), 
                 1 - t3, 
                 -t1 * xexpx(x/e, b) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}