findfixedstartvals <- function(form, data, cid, fct, fid, mform){
  if (is.na(cid)){
    mf <- model.frame(form, data)
    fullp <- coefficients(drm(form, data=mf[,1:2], fct=fct))
    return(fullp)
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
    return(unlist(plist))
  }
}



findrandomstartvals <- function(form, data, fct, random, curveid){
  if (class(random)[1] == "list"){
    rid <- names(random)
    rpars <- lapply(random, function(x){
      spr <- strsplit(deparse(x[[2]]), "+")[[1]]
      return(spr[!(spr %in% c(" ", "+"))])
    })
  }
  if (class(random)[1] == "formula"){
    rid <- as.character(random[[3]][[3]])
    rid <- rid[rid != "/"]
    spr <- strsplit(deparse(random[[2]]), "+")[[1]]
    rpars <- list()
    for (i in 1:length(rid)){
      rpars[[i]] <- spr[!(spr %in% c(" ", "+"))]
    }
  }
  cid <- as.character(curveid)[3]
  if (is.na(cid)){
    cid <- "curveid"
    data$curveid <- as.factor(1)
  }   
  rform <- as.formula(paste(as.character(form)[2], as.character(form)[1], as.character(form)[3], "+", cid, "+", paste(rid, collapse="+")))
  mf <- model.frame(rform, data)
  
  rstart <- lapply(1:length(rid), function(i){    
    sf1 <- apply(mf[,3:(i+2), drop=FALSE], 1, function(xx) paste(xx, collapse="/"))
    sf2 <- apply(mf[,4:(i+3), drop=FALSE], 1, function(xx) paste(xx, collapse="/"))
    spl <- split(mf, sf2)
    rmat <- lapply(spl, function(x){ 
      aid <- unique(apply(x[,3:(i+2), drop=FALSE], 1, function(xn) paste(xn, collapse="/")))
      submf <- mf[sf1 %in% aid,,drop=FALSE]
      mff <- droplevels(mf[sf1 %in% aid,cid])
      spl3 <- split(submf, mff)
      pfc <- sapply(spl3, function(mfx) coefficients(drm(form, data=mfx, fct=fct)))
      spl4 <- split(x, droplevels(x[,cid]))
      co <- sapply(spl4, function(mmmm) coefficients(drm(form, data=mmmm, fct=fct)))
      return(apply(co-pfc, 1, function(p) p[which(min(abs(p)) == abs(p))]))
    })
    rmat <- t(simplify2array(rmat))
    colnames(rmat) <- fct$names
    rownames(rmat) <- names(spl)    
    rmat <- rmat[,rpars[[i]], drop=FALSE]
    return(rmat)
  })  
  names(rstart) <- rid
  return(rstart)  
}


