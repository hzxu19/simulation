
.n.g.M <- function(.z,.catData,.ML,.lML,.g){
  tmp <- rep(NA,.lML)
  #for(ml in 1:.lML) tmp[ml] <- sum(.z[.catData==.ML[ml],.g])
  tmp <- sapply(1:.lML,function(i) sum(.z[.catData==.ML[i],.g]))
  return(tmp)
}
n.g.M <- cmpfun(.n.g.M)
#n.g.M(z,.catData,1); sum(n.g.M(z,.catData,1))
##---------------------
.p.g.M <- function(.z,.catData,.N,.ML,.lML,.g,.n){
 tmp <- rep(NA,.lML)
 #for(ml in 1:.lML) tmp[ml] <- sum(.z[,.g]*(.catData==.ML[ml]))/sum(.z[,.g])
 tmp <- sapply(1:.lML,function(i) sum(.z[,.g]*(.catData==.ML[i]))/sum(.z[,.g]))
 return(tmp)
}
p.g.M <- cmpfun(.p.g.M)
#p.g.M(z,.catData,1); p.g.M(z,.catData,2)
##---------------------
.ybar.g.M <- function(.z,.normData,.catData,.N,.D,.ML,.lML,.g,.n){
  jnk <- .z[,.g]*.normData
  tmp <- matrix(NA,.lML,.D)
  rownames(tmp) <- .ML
  for(ml in 1:.lML){
   for(d in 1:.D)  tmp[ml,d] <- ifelse(.n[ml]>0,
                     apply(jnk[.catData==.ML[ml],d,drop=FALSE],2,sum)/.n[ml],
                     NA)
  }
 return(tmp)
}
ybar.g.M <- cmpfun(.ybar.g.M)
#ybar.g.M(.z,.normData,useCat[[1]],.N,.D,.MML[[1]],.lMML[1],1,n.g.M(.z,useCat[[1]],.MML[[1]],.lMML[m],1))
##---------------------
.gam.GM <- function(.z,.normData,.catDataCol,.N,.D,.M,.MML,.lMML,.G,.which.pri,.hyp){
 pro <- apply(.z,2,sum)/.N
 tmp <- array(NA,dim=c(max(.lMML),.M,.G)) ##MAY NEED TO CHANGE BACK TO 0
 if(is.null(.which.pri)){
   dnm <- rep(0,.M)
   num <- 0
 }else{
   dnm <- .lMML
   num <- 1
 }
    	for(g in 1:.G){
        for(m in 1:.M){
           for(ml in 1:.lMML[m]) tmp[ml,m,g] <- {sum(.z[.catDataCol[[m]]==.MML[[m]][ml],g])+num}/ {(.N*pro[g])+dnm[m]}
      }}
 	 #   }else{
   #  	 for(g in 1:.G){
   #      for(m in 1:.M){
   #       for(ml in 1:.lMML[m]) tmp[ml,m,g] <- sum(.z[.catDataCol[[m]]==.MML[[m]][ml],g],1)/ (.N*pro[g]+.lMML[m])
   #   }}}
     return(tmp)
}
gam.GM <- cmpfun(.gam.GM)
#gam.GM(z.init,.normData,mD$catData$Orig,2,NULL,kp)
##---------------------
##---------------------
##NOTE: Line "num1[is.na(num1)] <- 0"
## In some instances of a CJ model, it is possible for there to be
## no observations in for a Jth level in group G
## This line replaces any NA values with a 0 when this occurs
## so that that level is simply replace by the prior for mu
## If no prior is used, the NA remains and the model may not be able to be estimated
##---------------
.mix.mu.GM <- function(.z,.normData,.catDataCol,.catMap,.N,.D,.M,.MML,.lMML,.G,.which.pri,.hyp){
  baseMean <- array(NA,dim=c(.D,.G))
  baseN    <- array(NA,dim=c(1,.G))
  jnk   <- tmp <- array(NA,dim=c(.D,max(.lMML),.M,.G))
  jnkN  <- array(NA,dim=c(1,max(.lMML),.M,.G))
  nR    <- nrow(.catMap)
  tmp2  <- array(NA,dim=c(nR,.D,.G))
  
  for(g in 1:.G){
    ng <- n.g(.z,g)
    baseMean[,g] <- ybar.g(.z,.normData,.D,.N,g,ng)   
    for(m in 1:.M){
      n.gm <- n.g.M(.z,.catDataCol[[m]],.MML[[m]],.lMML[m],g)  
      jnk[1:.D,1:.lMML[m],m,g] <- t(ybar.g.M(.z,.normData,.catDataCol[[m]],.N,.D,.MML[[m]],.lMML[m],g,n.gm))
    }}
  
  if(is.null(.which.pri)){
    tmp <- jnk
  }else{ 
    for(g in 1:.G){
      for(m in 1:.M){
        n.gm <- n.g.M(.z,.catDataCol[[m]],.MML[[m]],.lMML[m],g) 
        for(ml in 1:.lMML[m]){
          num1 <- n.gm[ml]*jnk[,ml,m,g]
          num1[is.na(num1)] <- 0 
          num2 <- .hyp$kp.h*.hyp$mu.h
          tmp[,ml,m,g] <- {num1+num2} / {n.gm[ml] + .hyp$kp.h}
        }}}}
  if(.M==1){
    for(g in 1:.G) tmp2[,1:.D,g] <- tmp[1:.D,,.M,g]
  }else{
    .tmpDelta <- array(NA,dim=c(.D,.M,nrow(.catMap),.G));
    for(g in 1:.G){
      for(r in 1:nR){
        for(m in 1:.M) .tmpDelta[,m,r,g] <- tmp[1:.D,.catMap[r,m],m,g] - baseMean[,g]
      }}
    for(g in 1:.G) tmp2[,,g] <- matrix(rep(baseMean[,g],nR),ncol=.D,byrow=TRUE)
    for(g in 1:.G) tmp2[,,g] <- tmp2[,,g]+t(apply(.tmpDelta,c(1,3,4),sum)[,,g])# + baseMean[,g]
  }
  return(list(mean=tmp,meanMod=tmp2))
}
mix.mu.GM <- cmpfun(.mix.mu.GM)
#mix.mu.k(z,.normData,.catData,2,NULL,.01)
#mix.mu.k(z,.normData,.catData,2,"MOD",.01)
##---------------------
##.catData can be either of the "Orig" or "Repr" form
##.catMap contains l* rows and M columns which identify
##  all combination of levels from M multinomial variables
.mix.M.stepI <- function(.modelName,.normData,.catData,.catMap,.z,.N,.D,.M,.MML,.lMML,.G,.which.pri,.hyp){
  tmp.mult <- gam.GM(.z,.normData,.catData,.N,.D,.M,.MML,.lMML,.G,.which.pri,.hyp)
  tmp.m    <- my.m.step(.modelName,.z,.normData,.D,.N,.G,.which.pri,.hyp)
  return(list(pro=tmp.m$par$pro,mean=tmp.m$par$mean,variance=tmp.m$par$variance,mult=tmp.mult))
}
##---------------------
.mix.M.stepC <- function(.modelName,.normData,.catData,.catMap,.z,.N,.D,.M,.MML,.lMML,.G,.which.pri,.hyp){
  tmp.mult<- gam.GM(.z,.normData,.catData,.N,.D,.M,.MML,.lMML,.G,.which.pri,.hyp)
  tmp.mu  <- mix.mu.GM(.z,.normData,.catData,.catMap,.N,.D,.M,.MML,.lMML,.G,.which.pri,.hyp)
  tmp.var <- my.m.step(.modelName,.z,.normData,.D,.N,.G,.which.pri,.hyp)
  return(list(pro=tmp.var$par$pro,mean=tmp.mu$mean,meanMod=tmp.mu$meanMod,variance=tmp.var$par$variance,mult=tmp.mult))
}
mix.M.stepI <- cmpfun(.mix.M.stepI)
mix.M.stepC <- cmpfun(.mix.M.stepC)
##---------------------
#.par <- mix.M.stepI("EII",mD$norm,mD$catData$Orig,.z,2,which.pri=NULL,kp)
#.par <- mix.M.stepC("E",.normData,.catDataAll$Orig,.catDataAll$catMap,.z,.G,which.pri=NULL,kp)
#.parI$mean;.parJ$mean
 #for(m in .mclustModel) print(mix.M.step2(mod,.normData,.catDataCol,z.init,.G,NULL,.01)$pro)
#par$var$sigma
#m.EII(z,.normData,.G,.which.pri,.kp)
##------------------------------------
.useDF.I <- function(.modelName,.D,.G,.M,.lMML){
  multDF <- .G*{sum(.lMML)-.M}
  meanDF <- .G*.D
  proDF  <- .G-1
  varDF  <- nVarParams(.modelName, .D, .G) 
  tmp    <- proDF + multDF + meanDF + varDF
return(tmp)
}
.useDF.C <- function(.modelName,.D,.G,.M,.lMML){
  multDF <- .G*(sum(.lMML)-.M)
  meanDF <- .G*.D*sum(.lMML)
  proDF  <- .G-1
  varDF  <- nVarParams(.modelName, .D, .G) 
  tmp    <- proDF + multDF + meanDF + varDF
return(tmp)
}
useDF.I <- cmpfun(.useDF.I)
useDF.C <- cmpfun(.useDF.C)
##------------------------------------
.iDens <- function(g,.N,.D,.M,.MML,.lMML,.normData,.catDataCol,.par,.norm){
  .tmpC  <- matrix(0,.N,.M)
  .tmpND <- rep(NA,.N)
  if(.D==1){
    .tmpND <- eval(parse(text=paste(.norm,
                                       "(.normData,.par$mean[,",g,"],sqrt(.par$variance$sigmasq[",g,"]))",sep="")))
  }else{
    .tmpND <- eval(parse(text=paste(.norm,
                                       "(.normData,.par$mean[,",g,"],.par$variance$sigma[,,",g,"])",sep="")))
  } 
  for(m in 1:.M)  .tmpC[,m] <- ddiscrete(.catDataCol[[m]],prob=.par$mult[1:.lMML[m],m,g],value=.MML[[m]])
  .pro.F <- .par$pro[g] * {.tmpND*apply(.tmpC,1,prod)} 
  return(.pro.F)
} 
iDens <- cmpfun(.iDens)
##-----------------------------------------------
.mix.E.stepI <- function(.modelName,.normData,.catDataCol,.catMap,.N,.D,.M,.MML,.lMML,.G,.par,.dfFunc,.useRows){
  norm <- ifelse(.D==1,"dnorm","dmvnorm")
  if(.D==1 & .modelName=="E") .par$variance$sigmasq <- rep(.par$variance$sigmasq,.G)
  
  pro.F <- sapply(1:.G, function(i) iDens(i,.N,.D,.M,.MML,.lMML,.normData,.catDataCol,.par,norm))  
  #pro.F[which(is.na(pro.F))] <- 0 ##Should be unnecessary now, remove.
  if(!any(apply(pro.F,1,function(x) sum(x,na.rm=TRUE))==0) | any(is.na(pro.F))){
    pro.F1 <- apply(pro.F,1,sum)
    llik   <- sum(log(pro.F1))
    z      <- pro.F/pro.F1
    class  <- map(z)
    npar   <- .dfFunc(.modelName,.D,.G,.M,.lMML)
    bic    <- 2*llik - npar*log(.N)
   ##Add ERROR FLAG
    eFlag <- ifelse(length(na.omit(unique(class))) != .G,"Error - Unique classes NOT EQUAL to G",NA)
  }else{
    llik   <- -Inf
    z      <- matrix(NA,.N,.G)
    class  <- rep(NA,.N)
    npar   <- NA
    bic    <- NA
    eFlag  <- "0 density"
  }
  return(list(z=z,class=class,loglik=llik,bic=bic,npar=npar,par=.par,eFlag=eFlag))#,Cat=.tmpCD,Nom=.tmpND))
}
mix.E.stepI <- cmpfun(.mix.E.stepI)
##------------------------------------
.cDens <- function(g,.N,.D,.M,.MML,.lMML,.normData,.catDataCol,.catMap,.par,.useRows,.norm){
  .tmpCD <- array(0,dim=c(.N,.M));  
  .tmpND <- rep(NA,.N)
  for(r in 1:nrow(.catMap)){
    useRows <- which(.useRows==r)
    if(.D==1){
      .tmpND[useRows] <- eval(parse(text=paste(.norm,
                                                "(.normData[useRows,],.par$meanMod[",r,",,",g,"],sqrt(.par$variance$sigmasq[",g,"]))",sep="")))
    }else{
      .tmpND[useRows] <- eval(parse(text=paste(.norm,
                                                "(.normData[useRows,],.par$meanMod[",r,",,",g,"],.par$variance$sigma[,,",g,"])",sep="")))
  }}
  for(m in 1:.M) .tmpCD[,m] <- ddiscrete(.catDataCol[[m]],prob=.par$mult[1:.lMML[m],m,g],value=.MML[[m]])
  .pro.F <- .par$pro[g] * {.tmpND * apply(.tmpCD,1,prod)}
  return(.pro.F)
}
cDens <- cmpfun(.cDens)
##------------------------------------
.mix.E.stepC <- function(.modelName,.normData,.catDataCol,.catMap,.N,.D,.M,.MML,.lMML,.G,.par,.dfFunc,.useRows){
  norm <- ifelse(.D==1,"dnorm","dmvnorm")
  if(.D==1 & .modelName=="E") .par$variance$sigmasq <- rep(.par$variance$sigmasq,.G)
  pro.F  <- sapply(1:.G, function(g) cDens(g,.N,.D,.M,.MML,.lMML,.normData,.catDataCol,.catMap,.par,.useRows,norm))
  if(!any(apply(pro.F,1,function(x) sum(x,na.rm=TRUE))==0) | any(is.na(pro.F))){
    pro.F1 <- apply(pro.F,1,sum)
    llik   <- sum(log(pro.F1))
    z      <- pro.F/pro.F1
    class  <- map(z)
    npar   <- .dfFunc(.modelName,.D,.G,.M,.lMML)
    bic    <- 2*llik - npar*log(.N)
  ##Add ERROR FLAG
    eFlag <- ifelse(length(na.omit(unique(class))) != .G,"Error - Unique classes NOT EQUAL to G",NA)
  }else{
    llik   <- -Inf
    z      <- matrix(NA,.N,.G)
    class  <- rep(NA,.N)
    npar   <- NA
    bic    <- NA
    eFlag  <- "0 density"
  }
  return(list(z=z,class=class,loglik=llik,bic=bic,npar=npar,par=.par,eFlag=eFlag))
  
}
mix.E.stepC <- cmpfun(.mix.E.stepC)
##------------------------------------
.mix.me <- function(.FullModelName,.normData,.catDataAll,.z,.G,.which.pri,.kp,.tol,.max.iter,.useRows){
  .N <- dim(.normData)[1]
  .D <- dim(.normData)[2]

  if(substr(.FullModelName,5,5)=="I"){ m.func <- mix.M.stepI; e.func <- mix.E.stepI; dfFunc <- useDF.I
   }else{                              m.func <- mix.M.stepC; e.func <- mix.E.stepC; dfFunc <- useDF.C }#"R" 
  
  if(substr(.FullModelName,6,6)=="I"){ useCat <- .catDataAll$Orig
   }else{                              useCat <- .catDataAll$Repr}#"J"
  
  if(!is.null(.which.pri)){  .hyp <- Hyper(.normData,.D,.G,.kp,.tol,.max.iter)
  }else{                     .hyp <- list(kp=.kp,tol=.tol,max.iter=.max.iter) }
  
  M <- dim(useCat)[2]
  ML <- vector("list",M) 
  lMML <- rep(NA,M)
  for(m in 1:M){
    ML[[m]] <- levels(useCat[[m]])
    lMML[m] <- length(ML[[m]])
  }
  if(!is.null(M)){
  .end <- ifelse(.D==1,1,3)
  .modelName <- unique(substr(.FullModelName,1,.end))
  m.init <- m.func(.modelName,.normData,useCat,.catDataAll$catMap,.z,.N,.D,M,ML,lMML,.G,.which.pri,.hyp)
  e.init <- e.func(.modelName,.normData,useCat,.catDataAll$catMap,.N,.D,M,ML,lMML,.G,m.init,dfFunc,.useRows)
 
  ##Begin iteration set-up
  iter <- 1
  .m <- .e <- vector("list",.max.iter)
  
  if(is.na(e.init$eFlag)){
   if(.G==1 & length(e.init$z)==0) e.init$z <- matrix(1,.N,1)
    .m[[iter]] <- m.func(.modelName,.normData,useCat,.catDataAll$catMap,e.init$z,.N,.D,M,ML,lMML,.G,.which.pri,.hyp)
     if(!is.na(.m[[iter]]$variance$sigma[1])){
       .e[[iter]] <- e.func(.modelName,.normData,useCat,.catDataAll$catMap,.N,.D,M,ML,lMML,.G,.m[[iter]],dfFunc,.useRows)
       tst        <- abs(.e[[iter]]$loglik - e.init$loglik)
     }else{ .e[[iter]] <- list(loglik=NA,bic=NA,z=NA,eFlag="ERROR - Variance Problems")}
   }else{ .e[[iter]] <- list(eFlag=e.init$eFlag) }
 
  if(is.na(.e[[iter]]$eFlag)){
   while(tst > .tol & iter < .max.iter & is.na(.e[[iter]]$eFlag)){
   iter <- iter+1

   if(.G==1 & length(.e[[iter-1]]$z) ==0) .e[[iter-1]]$z <- matrix(1,.N,1)
   .m[[iter]]   <- m.func(.modelName,.normData,useCat,.catDataAll$catMap,.e[[iter-1]]$z,.N,.D,M,ML,lMML,.G,.which.pri,.hyp)
    if(!is.na(.m[[iter]]$variance$sigma[1])){
     .e[[iter]] <- e.func(.modelName,.normData,useCat,.catDataAll$catMap,.N,.D,M,ML,lMML,.G,.m[[iter]],dfFunc,.useRows)
     tst        <- abs(.e[[iter]]$loglik - .e[[iter-1]]$loglik)
    }else{ .e[[iter]] <- list(loglik=NA,bic=NA,z=NA,eFlag="ERROR - Variance Problems")}
  }}else{  .e[[iter+1]] <- list(loglik=NA,bic=NA,z=NA,eFlag=.e[[iter]]$eFlag)}
 
  res <- .e[[iter]]
  #z <- res$z
  if(is.na(res$eFlag)){
    if(length(na.omit(unique(res$class))) != .G){
      res$loglik <- -Inf
      res$bic <- NA
    }
  }else{
    res$loglik <- -Inf
    res$bic    <- NA
  }
  return(list(z=res$z,class=res$class,loglik=res$loglik,bic=res$bic,npar=res$npar,par=res$par,iter=iter,eFlag=res$eFlag))
 }else{
  cat("ERROR: catData is 1D -- MAKE DATA FRAME \n")
}}
mix.me <- cmpfun(.mix.me)
##------------------------------------
#load("d1_150.rdata")
#.varTypes = varType1; 
#.normData=useData$norm; .catDataAll=useData$catData;.g1=2;.g2=6;.AllModelNames=tryMods 
#a <- mixClust(useData$norm,useData$catData,2,4,tryMods[1:4],"MOD",.kp,.tol,.max.iter)
#.which.pri="MOD";.kp=kp;.tol=tol;.max.iter=max.iter
.mixClust <- function(.normData,.catDataAll,.g1,.g2,.AllModelNames,.which.pri,.kp,.tol,.max.iter){
 if(nchar(.AllModelNames[1])==6){
  .D <- dim(.normData)[2]
  .M <- dim(.catDataAll$Orig)[2]
 if(.D==1){
  MmodelName <- NULL
  Cmodels <- unique(substr(.AllModelNames,4,6))
  if(any(substr(.AllModelNames,1,1) %in% c("E",.mclustModel[c(1,3,7)])))  MmodelName <- c(MmodelName,"EXX")
  if(any(substr(.AllModelNames,1,1) %in% c("V",.mclustModel[-c(1,3,7)]))) MmodelName <- c(MmodelName,"VXX")
  MmodelName <- paste(rep(MmodelName,each=length(Cmodels)),Cmodels,sep="")
  if(.M==1){    if(any(substr(.AllModelNames,6,6)=="J")) MmodelName <- MmodelName[-which(substr(MmodelName,6,6)=="J")] }	 
 }else{
  if(.M==1){
  MmodelName=.AllModelNames
  if(any(substr(.AllModelNames,6,6)=="J")) MmodelName <- .AllModelNames[-which(substr(.AllModelNames,6,6)=="J")]
   }else{ MmodelName <- .AllModelNames }}	
 
 MmodelName <- na.omit(MmodelName)
 if(length(MmodelName)>0){
  if(any(substr(MmodelName,5,5)=="R")){
   if(.M==1){  useRows <- as.numeric(.catDataAll$Orig[[1]])
    }else{     useRows <- as.numeric(.catDataAll$Repr[[1]]) }
  }else{       useRows <- NULL  }
  
  BIC <- matrix(NA,.g2-.g1+1,length(MmodelName))
  dimnames(BIC) <- list(paste(.g1:.g2),MmodelName)
  ##Initialization
  if(.D==1){ tmphc <- hc(modelName="E",.normData)
   }else{    tmphc <- hc(modelName="VII",.normData)   }
  tmp.cl <- hclass(tmphc,c(1:.g2))
 
  tmpMod <- vector("list",length(MmodelName))
  names(tmpMod) <- MmodelName
  tmpG <- vector("list",.g2)
  for(g in .g1:.g2){
    z.init <- unmap(tmp.cl[,g], groups=1:max(tmp.cl[,g]))
    tmpG[[g]] <- sapply(MmodelName, function(.mod) 	mix.me(.mod,.normData,.catDataAll,z.init,g,.which.pri,.kp,.tol,.max.iter,useRows),simplify=FALSE )
    BIC[paste(g),] <- sapply(MmodelName, function(.mod) tmpG[[g]][[.mod]]$bic)
  }

 if(all(is.na(BIC))){
   cat(paste("No mixClust model:",paste(MmodelName,collapse=","),"with",paste(.g1,":",.g2,sep=""),"clusters could be fit \n"))
   return(mc=list(bic=NA))
 }else{
  maxBIC <- which(BIC==max(BIC,na.rm=TRUE),arr.ind=TRUE)
  BIC <- matrix(as.numeric(BIC),dim(BIC))
  dimnames(BIC) <- list(paste(.g1:.g2),MmodelName)
  maxG   <- as.numeric(rownames(BIC)[maxBIC[1]])
  maxMod <- colnames(BIC)[maxBIC[2]]
  mc 	    <- tmpG[[maxG]][[maxMod]] #mix.me(maxMod, .normData,.catDataAll,z.init,maxG,.which.pri,.kp,.tol,.max.iter)
  mc$BIC   <- BIC
  mc$G     <- maxG
  mc$model <- maxMod
  cat(paste("best mixClust model:",maxMod,"with",maxG,"clusters \n"))
  return(mc)
 }}else{cat("No valid models available to fit \n")}
 }else{   cat("Invalid Model Name: Add .II,.RI,.IJ,.RJ for categorical \n") }
}
mixClust <- cmpfun(.mixClust)
