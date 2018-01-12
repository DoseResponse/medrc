#' The log-logistic function
#' 
#' 'meLL.X' provides the X-parameter log-logistic function.
#' 
#' The five-parameter logistic function is given by the expression 
#' \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}} 
#' or in another parameterisation 
#' \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-e)))^f}} 
#' The function is asymmetric for \eqn{f} different from 1.
#' 
#' @param x numeric dose vector.
#' @param b steepness
#' @param c lower limit
#' @param d upper limit
#' @param e ED50
#' @param f asymmetry
#' 
#' @references Finney, D. J. (1979) Bioassay and the Practise of Statistical Inference, \emph{Int. Statist. Rev.}, \bold{47}, 1--12.
#' 
#' @keywords models
#' 
#' @rdname meLL
meLL.2 <- function(x, b, e){
  .value <- 1/(1 + exp(b*log(x/e)))
  .actualArgs <- c("b", "e")
  t2 <- exp(b*(log(x) - log(e)))
  t5 <- (1 + t2)
  .grad <- cbind(-1 * xlogx(x/e, b, 2), 
                 divAtInf(t2, (1 + t2)^2) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meLL
meLL.3 <- function(x, b, d, e){
  .value <- d/(1 + exp(b*log(x/e)))
  .actualArgs <- c("b", "d", "e")
  t2 <- exp(b*(log(x) - log(e)))
  t5 <- (1 + t2)
  .grad <- cbind(-d * xlogx(x/e, b, 2), 
                 1/t5,
                 d * divAtInf(t2, (1 + t2)^2) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meLL
meLL.4 <- function(x, b, c, d, e){
  .value <- c + (d-c)/(1 + exp(b*log(x/e)))
  .actualArgs <- c("b", "c", "d", "e")
  t1 <- d - c
  t2 <- exp(b*(log(x) - log(e)))
  t5 <- (1 + t2)
  .grad <- cbind(-t1 * xlogx(x/e, b, 2), 
                 1-1/t5, 
                 1/t5,
                 t1 * divAtInf(t2, (1 + t2)^2) * b/e)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meLL
meLL.5 <- function(x, b, c, d, e, f){
  .value <- c + (d-c)/(1 + exp(b*log(x/e))^f)
  .actualArgs <- c("b", "c", "d", "e", "f")
  t1 <- d - c
  t2 <- exp(b*(log(x) - log(e)))
  t5 <- (1 + t2)^f
  .grad <- cbind(-t1 * xlogx(x/e, b, f+1) * f, 
                 1-1/t5, 
                 1/t5,
                 t1 * f * divAtInf(t2, (1 + t2)^(f+1)) * b/e,
                 -t1 * divAtInf(log(1 + t2), t5))
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}