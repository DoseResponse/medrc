mmplot <- function(x, ..., ndose=25, logx=FALSE){
  lllist <- list(x, ...)
  ismedrc <- sapply(lllist, function(x) inherits(x, "medrc") | inherits(x, "glsdrc"))
  mllist <- lllist[ismedrc]  
  Call <- match.call()
  Call$ndose <- NULL
  Call$logx <- NULL
  mnames <- as.character(Call[-1L])[ismedrc]

  yname <- as.character(x$form)[2]
  xname <- as.character(x$form)[3]
  if (is.null(x$curveid)){
    pform <- x$form
  } else {
    if (length(as.character(x$curveid)) == 2){
      fname <- as.character(x$curveid)[2]
    } else {
      fname <- as.character(x$curveid)[3]
    }
    pform <-  paste(yname, '~', xname, '+', fname)
  }
  mf <- model.frame(pform, data=x$data)
  
  if (logx == TRUE){
    if (min(mf[,2]) == 0){
      m0 <- mean(unique(mf[,2])[order(unique(mf[,2]))][1:2])
    } else {
      m0 <- min(mf[,2]) 
    }
    dr <- exp(seq(log(m0), log(max(mf[,2])), length=ndose))
  } else {
    dr <- seq(min(mf[,2]), max(mf[,2]), length=ndose)
  } 
   
  datlist <- lapply(lllist, function(x, ndose, logx, dr, mf){
    fct <- x$fct    
    if (is.null(x$curveid)){
      predictions <- x$fct$fct(dr, rbind(fixef(x$fit)))
      pdat <- data.frame(predictions, dose=dr) 
    } else {
      flev <- length(levels(mf[,3])) 
      cf <- matrix(coefficients(x), ncol=flev, byrow=TRUE)
      predictions <- stack(data.frame(apply(cf, 2, function(para) x$fct$fct(dr, rbind(para)))))$values
      pdat <- data.frame(predictions, dose=rep(dr, times=ncol(cf)), curve=rep(levels(mf[,3]), each=ndose))
    } 
    return(pdat)
  }, ndose=ndose, logx=logx, dr=dr, mf=mf)  
  
  pdat <- ldply(datlist)  
  pdat$model <- as.factor(rep(mnames, each=nrow(datlist[[1]])))
  
  if (is.null(x$curveid)){
    if (logx == TRUE){
      eval(parse(text=paste("ggplot(pdat, aes(x=dose, y=predictions, colour=model)) + coord_trans(x='log') + geom_line()", sep="")))
    } else {
      eval(parse(text=paste("ggplot(pdat, aes(x=dose, y=predictions, colour=model)) + geom_line()", sep="")))
    }
  } else {
    if (logx == TRUE){
      eval(parse(text=paste("ggplot(data=pdat, aes(x=dose, y=predictions, linetype=curve, colour=model)) + coord_trans(x='log') + geom_line()", sep="")))
    } else {
      eval(parse(text=paste("ggplot(data=pdat, aes(x=dose, y=predictions, linetype=curve, colour=model)) + geom_line()", sep="")))
    }
  }  
}
