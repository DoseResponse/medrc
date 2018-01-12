#' Log-normal dose-response model
#' 
#' \code{LN.X} and the accompanying convenience functions provide a general framework for specifying the mean function of the decreasing or incresing log-normal dose-response model.
#' 
#' For the case where log(ED50), denoted \eqn{e} in the equation below, is a parameter in the model, the mean function is:
#' \deqn{f(x) = c + (d-c)(\Phi(b(\log(x)-e)))}
#' and the mean function is:
#' \deqn{f(x) = c + (d-c)(\Phi(b(\log(x)-\log(e))))}  
#' in case ED50, which is also denoted \eqn{e}, is a parameter in the model. If the former model is fitted any estimated ED values will need to be back-transformed subsequently in order to obtain effective doses on the original scale.
#' The mean functions above yield the same models as those described by Bruce and Versteeg (1992), but using a different parameterisation (among other things the natural logarithm is used).
#' For the case \eqn{c=0} and \eqn{d=1}, the log-normal model reduces the classic probit model (Finney, 1971) with log dose as explanatory variable (mostly used for quantal data). This special case is available through the convenience function \code{meLN.2}.
#' The case \eqn{c=0} is available as the function \code{meLN.3}. The full four-parameter model is available through \code{meLN.4}.
#' 
#' @param x numeric dose vector
#' @param b steepness
#' @param c lower limit
#' @param d upper limit
#' @param e ED50
#' 
#' @references 
#' Finney, D. J. (1971) \emph{Probit analysis}, London: Cambridge University Press.
#' Bruce, R. D. and Versteeg, D. J. (1992) A statistical procedure for modeling continuous toxicity data, \emph{Environ. Toxicol. Chem.}, \bold{11}, 1485--1494.
#' 
#' @keywords models
#' 
#' @rdname meLN
meLN.2 <- function(x, b, e){
  .expr4 <- log(x) - log(e)
  .expr5 <- b * .expr4
  .expr9 <- dnorm(.expr5)
  .value <- pnorm(.expr5)
  .grad <- array(0, c(length(.value), 2L), list(NULL, c("b", "e")))
  tempVec <- .expr9 * .expr4
  tempVec[!is.finite(tempVec)] <- 0
  .grad[, "b"] <- tempVec
  .grad[, "e"] <- -(.expr9 * (b * (1/e)))
  attr(.value, "gradient") <- .grad
  .value
}

#' @rdname meLN
meLN.3 <- function(x, b, d, e){
  .expr1 <- d
  .expr4 <- log(x) - log(e)
  .expr5 <- b * .expr4
  .expr6 <- pnorm(.expr5)
  .expr9 <- dnorm(.expr5)
  .value <- .expr1 * .expr6
  .grad <- array(0, c(length(.value), 3L), list(NULL, c("b", "d", "e")))
  tempVec <- .expr9 * .expr4
  tempVec[!is.finite(tempVec)] <- 0
  .grad[, "b"] <- .expr1 * tempVec
  .grad[, "d"] <- .expr6
  .grad[, "e"] <- -(.expr1 * (.expr9 * (b * (1/e))))
  attr(.value, "gradient") <- .grad
  .value
}

#' @rdname meLN
meLN.4 <- function(x, b, c, d, e){
  .expr1 <- d - c
  .expr4 <- log(x) - log(e)
  .expr5 <- b * .expr4
  .expr6 <- pnorm(.expr5)
  .expr9 <- dnorm(.expr5)
  .value <- c + .expr1 * .expr6
  .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
  tempVec <- .expr9 * .expr4
  tempVec[!is.finite(tempVec)] <- 0
  .grad[, "b"] <- .expr1 * tempVec
  .grad[, "c"] <- 1 - .expr6
  .grad[, "d"] <- .expr6
  .grad[, "e"] <- -(.expr1 * (.expr9 * (b * (1/e))))
  attr(.value, "gradient") <- .grad
  .value
}

#' @rdname meLN
meLLN.4 <- function(x, b, c, d, e){
  .expr1 <- d - c
  .expr3 <- log(x) - e
  .expr4 <- b * .expr3
  .expr5 <- pnorm(.expr4)
  .expr8 <- dnorm(.expr4)
  .value <- c + .expr1 * .expr5
  .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
  tempVec <- .expr8 * .expr3
  tempVec[!is.finite(tempVec)] <- 0
  .grad[, "b"] <- .expr1 * tempVec
  .grad[, "c"] <- 1 - .expr5
  .grad[, "d"] <- .expr5
  .grad[, "e"] <- -(.expr1 * (.expr8 * b))
  attr(.value, "gradient") <- .grad
  .value
}