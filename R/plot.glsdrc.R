plot.glsdrc <- function(x, ..., ndose=25, logx=FALSE){
  fct <- x$fct
  
  yname <- as.character(x$form)[2]
  xname <- as.character(x$form)[3]
  if (is.null(x$curveid)){
    pform <- x$form
  } else {
    fname <- as.character(x$curveid)[3]
    pform <- paste(yname, "~", xname, "+", fname)
  }
  mf <- model.frame(pform, data=x$data)
  
  if (logx == TRUE){
    if (min(mf[,2]) == 0){
      m0 <- mean(unique(mf[,2])[order(unique(mf[,2]))][1:2])
    } else {
      m0 <- min(mf[,2]) 
    }    
    dr <- exp(seq(log(m0), log(max(mf[,2])), length=ndose))
    mf[mf[,2] == 0, 2] <- m0
  } else {
    dr <- seq(min(mf[,2]), max(mf[,2]), length=ndose)
  } 
  
  if (is.null(x$curveid)){
    predictions <- x$fct$fct(dr, rbind(coefficients(x$fit)))
    pdat <- data.frame(predictions, dose=dr)   
    if (logx == TRUE){
      eval(parse(text=paste("ggplot(mf, aes(x=",xname,", y=",yname,")) + coord_trans(x='log') + geom_point(alpha=0.3) + geom_line(data=pdat, aes(x=dose, y=predictions), colour='blue3')", sep="")))
    } else {
      eval(parse(text=paste("ggplot(mf, aes(x=",xname,", y=",yname,")) + geom_point(alpha=0.3) + geom_line(data=pdat, aes(x=dose, y=predictions), colour='blue3')", sep="")))
    }
  } else {
    flev <- length(levels(mf[,3])) 
    cf <- matrix(coefficients(x), ncol=flev, byrow=TRUE)
    predictions <- stack(data.frame(apply(cf, 2, function(para) x$fct$fct(dr, rbind(para)))))$values
    pdat <- data.frame(predictions, dose=rep(dr, times=ncol(cf)), curve=rep(levels(mf[,3]), each=ndose))
    if (logx == TRUE){
      eval(parse(text=paste("ggplot(mf, aes(x=",xname,", y=",yname,", colour=", fname,")) + coord_trans(x='log') + geom_point(alpha=0.3) + geom_line(data=pdat, aes(x=dose, y=predictions, colour=curve))", sep="")))
    } else {
      eval(parse(text=paste("ggplot(mf, aes(x=",xname,", y=",yname,", colour=", fname,")) + geom_point(alpha=0.3) + geom_line(data=pdat, aes(x=dose, y=predictions, colour=curve))", sep="")))
    }
  }    
}

