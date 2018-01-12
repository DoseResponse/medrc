#' The logistic model
#' 
#' The general asymmetric five-parameter logistic model for describing dose-response relationships.
#' 
#' The default arguments yields the five-parameter logistic mean function given by the expression \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f}} The model is sometimes referred to as the Boltzmann model.
#' 
#' @param x numeric dose vector.
#' @param b steepness
#' @param c lower limit
#' @param d upper limit
#' @param e ED50
#' @param f asymmetry
#' 
#' @keywords models
#' 
#' @rdname meL
meL.2 <- function(x, b, e){
  .value <- 1/((1 + exp(b * (x - e))))
  .actualArgs <- c("b", "e")
  t2 <- exp(b*(x - e))
  t4 <- ((1 + t2)^(-2))
  t5 <- (1 + t2)
  .grad <- cbind(-1 * t2 * t4 * (x - e), 
                 t2 * t4 * b)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meL
meL.3 <- function(x, b, d, e){
  .value <- d/((1 + exp(b * (x - e))))
  .actualArgs <- c("b", "d", "e")
  t2 <- exp(b*(x - e))
  t4 <- ((1 + t2)^(-2))
  t5 <- (1 + t2)
  .grad <- cbind(-d * t2 * t4 * (x - e), 
                 1/t5, 
                 d * t2 * t4 * b)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meL
meL.4 <- function(x, b, c, d, e){
  .value <- c + (d - c)/((1 + exp(b * (x - e))))
  .actualArgs <- c("b", "c", "d", "e")
  t1 <- d - c
  t2 <- exp(b*(x - e))
  t4 <- ((1 + t2)^(-2))
  t5 <- (1 + t2)
  .grad <- cbind(-t1 * t2 * t4 * (x - e), 
                 1 - 1/t5, 
                 1/t5, 
                 t1 * t2 * t4 * b)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}

#' @rdname meL
meL.5 <- function(x, b, c, d, e, f){
  .value <- c + (d - c)/((1 + exp(b * (x - e)))^f)
  .actualArgs <- c("b", "c", "d", "e", "f")
  t1 <- d - c
  t2 <- exp(b*(x - e))
  t4 <- f * ((1 + t2)^(-1 * f - 1))
  t5 <- (1 + t2)^f
  .grad <- cbind(-t1 * t2 * t4 * (x - e), 
                 1 - 1/t5, 
                 1/t5, 
                 t1 * t2 * t4 * b,
                 -t1 * log(1 + t2)/t5)
  dimnames(.grad) <- list(NULL, .actualArgs)
  attr(.value, "gradient") <- .grad
  return(.value)
}