divAtInf <- function(x, y){
  retVec <- x/y
  retVec[(!is.finite(y))] <- 0
  retVec
}

xexplogx <- function(x, p){
  lv <- (x < 1e-12)
  nlv <- !lv
  rv <- rep(0, length(x))
  xlv <- x[lv]
  rv[lv] <- 0
  xnlv <- x[nlv]
  rv[nlv] <- log(xnlv) * (xnlv^p[nlv]) * exp(-(xnlv^p[nlv]))
  rv
}


xexpx <- function(x, p){
  lv <- (x < 1e-12)
  nlv <- !lv
  rv <- rep(0, length(x))
  xlv <- x[lv]
  rv[lv] <- 0
  xnlv <- x[nlv]
  rv[nlv] <- (xnlv^p[nlv]) * exp(-(xnlv^p[nlv]))
  rv
}

xlogx <- function(x, p, f = 0){
  lv <- (x < 1e-12)
  nlv <- !lv
  rv <- rep(0, length(x))
  xPowerp <- x^p
  ratioVec <- divAtInf(xPowerp, (1 + xPowerp)^f)
  xlv <- x[lv]
  rv[lv] <- log(xlv^ratioVec[lv])
  xnlv <- x[nlv]
  rv[nlv] <- ratioVec[nlv] * log(xnlv)
  rv
}