## logistic
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



## log-logistic
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


## Weibull
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


## log normal
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


## Asymptotic regression
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


## exponential decay
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

### fplogistic
#FPL.4 <- function(x, b, c, d, e, p1, p2) {
#  .expr1 <- d - c
#  .expr3 <- log(x + 1)
#  .expr4 <- .expr3^p1
#  .expr6 <- .expr3^p2
#  .expr9 <- exp(b * .expr4 + e * .expr6)
#  .expr10 <- 1 + .expr9
#  .expr15 <- .expr10^2
#  .expr18 <- 1/.expr10
#  .value <- c + .expr1/.expr10
#  .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
#  .grad[, "b"] <- -(.expr1 * (.expr9 * .expr4)/.expr15)
#  .grad[, "c"] <- 1 - .expr18
#  .grad[, "d"] <- .expr18
#  .grad[, "e"] <- -(.expr1 * (.expr9 * .expr6)/.expr15)
#  attr(.value, "gradient") <- .grad
#  .value
#}
