LL4.loged <- function(respLev){
    paste(log(respLev/(1-respLev)), "/b+log(e)")
}

LN4.loged <- function(respLev){
    paste(qnorm(1-respLev), "/b+log(e)")
}

W14.loged <- function(respLev){
    paste(log(-log(1-respLev)), "/b+log(e)")
}

W24.loged <- function(respLev){
    paste(log(-log(respLev)), "/b+log(e)")
}
