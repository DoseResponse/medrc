findfixedstartvalsglsdrm <- function(form, data, cid, fct, fid, mform){
  if (is.na(cid)){
    #mf <- model.frame(form, data)
    #stv <- fct$ssfct(mf)
    #plist <- lapply(stv, function(x) x)
    #names(plist) <- fct$names
    #return(plist)
    return(NULL)
  } else {
    cform <- as.formula(paste(as.character(form)[2], as.character(form)[1], as.character(form)[3], "+", cid))
    mf <- model.frame(cform, data)
    spl <- split(mf, mf[,3])
    fullp <- coefficients(drm(form, data=mf[,1:2], fct=fct))
    pmat <- t(sapply(spl, function(x){
      trycoef <- try(coefficients(drm(form, data=x, fct=fct)), silent=TRUE)  
      if (class(trycoef)[1] == "try-error") return(fullp) else return(trycoef)
    }))    
    plist <- lapply(1:length(fid), function(i){
      if (fid[i]) fullp[i] else pmat[,i]   
    }) 
    names(plist) <- fct$names
    return(plist)
  }
}
