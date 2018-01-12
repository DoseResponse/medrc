#' ED and BMD estimation on log scale
#' 
#' Calculate BMD estimates on log transformed scale by the dmlist argument in mmaED.
#' 
#' @param respLev response level
#' 
#' @keywords models
#' 
#' @rdname logED
LL4.loged <- function(respLev){
    paste(log(respLev/(1-respLev)), "/b+log(e)")
}
#' @rdname logED
LN4.loged <- function(respLev){
    paste(qnorm(1-respLev), "/b+log(e)")
}
#' @rdname logED
W14.loged <- function(respLev){
    paste(log(-log(1-respLev)), "/b+log(e)")
}
#' @rdname logED
W24.loged <- function(respLev){
    paste(log(-log(respLev)), "/b+log(e)")
}
