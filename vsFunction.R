##----------------------
##             Data Prep: newCat, makeData
##----------------------
##Function creates a single categorical variable from
## multiple variables
## Can be used when >1 cat variables have same role
## MUST be used when >1 cat variables have same role in presence
##  of normal variables
.newCat <- function(.catData,.N,.M){
  if(.M>1){
   catData = matrix(NA,.N,.M)
   for(m in 1:.M) catData[,m] <- factor(.catData[,m])
   levelTable <- eval(parse(text=paste("table(",
	               paste("catData[,",1:.M,"]",sep="",collapse=","),")")))
   names(dimnames(levelTable)) <- paste("d",1:.M)
   expTable <- expand.table(levelTable)
   uniq <- unique(expTable)
   nR = nrow(uniq)
   
   .newVar = NULL
   for(i in 1:nR){
	  .tmp = NULL
    for(c in 1:ncol(expTable)) .tmp <- c(.tmp,
	               paste("catData[,",c,"]==uniq[",i,",",c,"]",sep=""))
    chkTmp <- paste(.tmp,collapse=" & ")
    who <- eval(parse(text=paste("which(",chkTmp,")",sep="")))
    .newVar[who] <- i
   }
   newVar = matrix(.newVar,.N,1)
   rownames(uniq) <- paste(1:nR)
  }else{
   .newVar = factor(.catData)
   newVar = .catData
   uniq = data.frame(matrix(levels(.newVar),length(levels(.newVar)),1))
   rownames(uniq) = levels(newVar)
  }
   return(list(newCat=as.factor(newVar),catMap=uniq))
 }
newCat <- cmpfun(.newCat)
##---------------------------
.makeData <- function(.propVar,.baseVar,.InOrOut,.varTypes,.fullData){
  N = dim(.fullData)[1]
  useVar = .baseVar
  if(.InOrOut=="Incl")   useVar = c(.propVar,useVar)
  useVarTF = 1:length(.fullData) %in% useVar
  useTypes = .varTypes[useVarTF]
  nN = sum(useTypes=="N")
  nC = sum(useTypes=="C")

  if(nN>0){
   tmpNorm = matrix(NA,N,nN)
   for(i in 1:nN) tmpNorm[,i] <- .fullData[[which(useVarTF & .varTypes=="N")[i]]]
  }else{   tmpNorm = NULL   }
  if(nC>0){
   tmpCat1 = list(); tmpMat = matrix(0,N,nC)
   for(i in 1:nC) tmpMat[,i] <- tmpCat1[[i]] <- factor(as.numeric(factor(.fullData[[which(useVarTF & .varTypes=="C")[i]]])))
   tmpCat <- newCat(tmpMat,N,nC)
  }else{
   tmpCat = list(newCat=NULL,catMap=NULL)
   tmpCat1 = tmpCat
  }
    return(list(norm=tmpNorm,catData=list(Orig=data.frame(tmpCat1),Repr=data.frame(tmpCat$newCat),catMap=tmpCat$catMap),useTypes = useTypes)) 
}
makeData <- cmpfun(.makeData)      #save time  (byte code compiler)
#rm(.propVar);rm(.baseVar);rm(.InOrOut)
#md = makeData(NULL,c(1,5),"Excl",.varTypes,d[[1]])
##-------------------------------
##              Cluster Functions: m2m, whichClust
##-------------------------------
## Transform results from mixmodCluster
##  to be of same form as results from:
##  Mclust, Mclustpri (mine), mixClust (mine)
m2m <- function(mx){
 res = list(modelname=NULL,n=NULL,d=NULL,
            d=NULL,parameters=list(proportions=NULL,
			mean=NULL,variance=NULL),
            bic=NULL,loglik=NULL,class=NULL,z=NULL)
 res$modelname <- mx['bestResult']['model']
 res$n <- length(mx['bestResult']['partition'])
 res$d <- ncol(mx['bestResult']['parameters']['mean'])
 res$parameters$proportions <- mx['bestResult']['parameters']['proportions']
 res$parameters$mean <- t(mx['bestResult']['parameters']['mean'])
 res$parameters$variance <- mx['bestResult']['parameters']['variance']
 res$bic <- -mx['bestResult']['criterionValue']
 res$loglik <- mx['bestResult']['likelihood']
 res$class <- mx['bestResult']['partition']
 res$z <- mx['bestResult']['proba']
return(res)
}
##-------------------------------
##whichClust: Function
## Function Dependencies:
##   useData (user Function),
##   Rmixmod (need R>= 2.15.1)
##   and Mclust (potentially), not necessarily?
## Variable Dependencies:
##   propVar  (variable proposed for Inclusion/Exclusion)
##   remnVar  (Remaining variables currently INCLUDED)
##   fullData (data frame of complete data, regardless of inclusion/exclusion)
##   G        (Current number of groups)
##   varTypes ("C" - Categorical, "N" - Continuous)
##   FullModelName ("Mixed-Data" model name - Create a 'matching' function for
##                names in Normal / Discrete settings)
##   which.pri, kp, tol, max.iter 
#.propVar = var$propVar; .remnVar = var$remnVar; .step="Incl"; .FullModelName=tryMods
#.propVar=NULL;.remnVar=1:7;.step="Excl"; .modelName="EII"
#whichClust(NULL,.clustVar,"Excl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
#whichClust(.propVar,.remnVar,.step,.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
#.propVar=NULL; .remnVar =c(1,2,5);.step="Excl";.fullData=d[d[,8]==nS,.startVar];.g1=g1.g2=g2;
#varTypes;.mod=tryMods1;.mc=TRUE
##-------------------------------------------------
.whichClust <- function(.propVar,.remnVar,.step,.fullData,.g1,.g2,.varTypes,
                       .FullModelName,.which.pri,.kp,.tol,.max.iter,.mc){
  clustFlag <- NA
  useData = makeData(.propVar,.remnVar,.step,.varTypes,.fullData)
  D = dim(useData$norm)[2]
  if(all(useData$useTypes=="N")){				##-- Beg: Normal
   if(.mc){								              ##-- Beg: Mclust Procedure    
	  .tmpClust <- Mclustpri(useData$norm,.g1,.g2,unique(substr(.FullModelName,1,3)),.which.pri,.kp,.tol,.max.iter)	      # ?						  			     			   
   }else{ ##-- Use Mixmod	  		  		 	##-- End: Mclust Procedure, Beg: Mixmod, Procedure
    if(D==1){			  		  		  		  	##-- Beg: uni/multivariate
    	.NmodelName <- ifelse(substr(.FullModelName,1,3)=="EII","Gaussian_pk_L_I","Gaussian_pk_Lk_I") ##Simple restriction
    }else{
    	.NmodelName <- .mixmodModel[which(.mclustModel==substr(.FullModelName,1,3))]
    }			  					  					  		##-- End: uni/multivariate
    mods = new("GaussianModel",listModels=c(.NmodelName))
    .tmpClust <- m2m(mixmodCluster(useData$norm,.g1:.g2,models=mods))      #Rmixmod ???????
	  if(length(.tmpClust$modelname)==0){	##-- Beg: Error
   	 .tmpClust <- list(bic=NA); 
	   clustFlag <- "ERROR: Estimation Problems"
   }}		  					  					  			##-- End: Error, End: Mixmod Procedure
  }else{				  					  					##-- End: Normal, Beg: Other1
    if(all(c("C","N") %in% useData$useTypes)){##-- Beg: Mixed
     .tmpClust <- mixClust(useData$norm,useData$catData,.g1,.g2,.FullModelName,.which.pri,.kp,.tol,.max.iter)
    }else{							   					       ##-- End: Error, Mixed, Beg: Other2
      if(all(useData$useTypes=="C")){				 ##-- Beg: Cat
       clustFlag <- "ERROR: Only Categorical variables"
       .tmpClust = list(bic=NA)
  }}}		  				  				  				  	   ##-- End: Cat, Other1, Other2
  .tmpClust$clustFlag = clustFlag
  .tmpClust$data = useData
return(.tmpClust)
}
whichClust <- cmpfun(.whichClust)
##-----------------------------------------------------------
##-------------------------------
##-------------------------------
## REGRESSION FUNCTIONS: catBIC, regFunc, stepReg
##-------------------------------
##-------------------------------
.regFunc <- function(.regDep,.regInd,.varTypes,.fullData,.N){
  dF = as.data.frame(.fullData)
  .typeof.regDep = .varTypes[.regDep]
  if(.typeof.regDep=="N") useReg <- c("lm")
  if(.typeof.regDep=="C") useReg <- c("multinom")
  if(length(.regInd)>0){
   n.indTypes = which(.varTypes[.regInd]=="N")
   c.indTypes = which(.varTypes[.regInd]=="C")
   cov = NULL
   if(length(n.indTypes)>0) cov = paste("useDF[[",.regInd[n.indTypes],"]]",collapse="+",sep="")
   if(length(c.indTypes)>0) cov = paste(c(cov,paste("factor(useDF[[",.regInd[c.indTypes],"]])",collapse="+",sep="")),collapse="+",sep="")
   form = paste("useDF[[",.regDep,"]] ~", cov ,sep="")
   useDF = dF
  }else{
   form = paste("useDF[,1] ~ 1",sep="")
   useDF = as.data.frame(matrix(dF[[.regDep]],.N,1))
  }
  reg = eval(parse(text=paste(useReg,"(formula(",form,"),data=useDF)",sep="")))
  bic <- BIC(reg)
  return(list(useDep=.regDep,useInd=.regInd,reg=reg,bic=-bic))
}
regFunc <- cmpfun(.regFunc)
##----------------------
.step.regFunc <- function(.regDep,.regInd,.varTypes,.fullData,.N,.clustModel){
  dF = as.data.frame(.fullData)
  .typeof.regDep = .varTypes[.regDep]
  if(.typeof.regDep=="N"){
    useReg <- c("lm")
    if(substr(.clustModel,5,5)=="I") .regInd = .regInd[-which(.varTypes[.regInd]=="C")]
  } 
  if(.typeof.regDep=="C"){
    useReg <- c("multinom")
    if(substr(.clustModel,5,5)=="I") .regInd = NULL ##Else, .regInd stays the same
  }  
  if(length(.regInd)>0){
    n.indTypes = which(.varTypes[.regInd]=="N")
    c.indTypes = which(.varTypes[.regInd]=="C")
    cov = NULL
    if(length(n.indTypes)>0) cov = paste("useDF[[",.regInd[n.indTypes],"]]",collapse="+",sep="")
    if(length(c.indTypes)>0) cov = paste(c(cov,paste("factor(useDF[[",.regInd[c.indTypes],"]])",collapse="+",sep="")),collapse="+",sep="")
    form = paste("useDF[[",.regDep,"]] ~", cov ,sep="")
    useDF = dF
  }else{
    form = paste("useDF[,1] ~ 1",sep="")
    useDF = as.data.frame(matrix(dF[[.regDep]],.N,1))
  }
  reg = eval(parse(text=paste(useReg,"(formula(",form,"),data=useDF)",sep="")))
  bic <- BIC(reg)
  return(list(useDep=.regDep,useInd=.regInd,reg=reg,bic=-bic))
}
step.regFunc <- cmpfun(.step.regFunc)
##----------------------
excludeCriteria.r = "any(changeBIC.r<0)"
includeCriteria.r = "any(changeBIC.r>0)"
.stepReg <- function(.regDep,.regInd,.varTypes,.fullData,.N){
 dec= ""
 useInd <- .regInd
 useVarTF.r <- rep(TRUE,length(.regInd))
 vExcl = vIncl = changBIC.r = NULL

 ##Regress Dependent variable on all possible independent variables
 allBIC <- regFunc(.regDep,useInd,.varTypes,.fullData,.N)
 tmpBIC <- NULL
 ##For j in 1:length(independent vars), remove jth var, regress dependent var on remaining dependent vars
 step <- "Excl"
 for(j in 1:length(useInd[useVarTF.r]))   tmpBIC[j] <- regFunc(.regDep,useInd[useVarTF.r][-j],.varTypes,.fullData,.N)$bic
 ##Calculate BIC difference
 changeBIC.r <- allBIC$bic - tmpBIC
 crit <- ifelse(step=="Excl",excludeCriteria.r,includeCriteria.r)
 ##Begin stepwise with exclusion step:
 if(eval(parse(text=crit))){
    stepDir <- ifelse(step=="Excl","min","max")
    changeVar.r = eval(parse(text=paste("which(useVarTF.r)[which.",stepDir,"(changeBIC.r)]",sep="")))
    useVarTF.r[changeVar.r] = !useVarTF.r[changeVar.r]
    exit <- ifelse(sum(useVarTF.r==TRUE)==0,TRUE,FALSE)
 }else{ exit <- TRUE }
 ##Enter stepwise: Begin with exclusion step
 step <- "Excl" ;
 while(exit==FALSE){ 
  ##Regress Dep on remaining Independent variables
  allBIC <- regFunc(.regDep,useInd[useVarTF.r],.varTypes,.fullData,.N)
  tmpBIC <- changeBIC.r <- NULL
   ##For j in 1:length(remainingindependent vars),
   ## remove jth var, regress dependent var on remaining dependent vars
  if(step=="Excl"){ 
	 for(j in 1:length(useInd[useVarTF.r]))   tmpBIC[j] <- regFunc(.regDep,useInd[useVarTF.r][-j],.varTypes,.fullData,.N)$bic
  }else{
	 for(j in 1:length(useInd[!useVarTF.r]))  tmpBIC[j] <- regFunc(.regDep,c(useInd[useVarTF.r],useInd[!useVarTF.r][j]),.varTypes,.fullData,.N)$bic
  }
  if(step=="Excl") sign = c(1,-1)
  if(step=="Incl") sign = c(-1,1)	
  changeBIC.r <- sign[1]*allBIC$bic + sign[2]*tmpBIC; changeBIC.r
  crit <- ifelse(step=="Excl",excludeCriteria.r,includeCriteria.r)

 if(eval(parse(text=crit))){                  ##Begin: Criteria is met
    stepDir <- ifelse(step=="Excl","min","max")
    lastChangeVar.r = changeVar.r
    .usevar <- ifelse(step=="Excl","useVarTF.r","!useVarTF.r")
    changeVar.r = eval(parse(text=paste("which(",.usevar,")[which.",stepDir,"(changeBIC.r)]",sep="")))
    useVarTF.r[changeVar.r] = !useVarTF.r[changeVar.r]
    
    tst.r = lastChangeVar.r==changeVar.r
    if(!is.na(tst.r) & tst.r){	exit <- TRUE  
     }else{		                  step <- ifelse(step=="Excl","Incl","Excl")   }
    dec = "accept"
  }else{ ##If exclusion criteria is not met (i.e. ChangeVar.r = NA) -- go to Inclusion OR exit
         ##If some variables for inclusion (!useVarTF) exist, go to inclusion step
    if(dec!="reject"){                        ##Rejected two steps in a row
     if(step=="Excl"){                        ##Rejected two steps, current step=Excl
      if(sum(!useVarTF.r)>0){	  step <- ifelse(step=="Excl","Incl","Excl")
       }else{               	  exit <- TRUE }
     }else{                                   ##Else, current step=="Incl"
      if(sum(useVarTF.r)>0){    step <- ifelse(step=="Excl","Incl","Excl")
       }else{                   exit <- TRUE }
     }                                        ##End step check.
    }else{                      exit <- TRUE }##Last step was rejected 	  
  }                                           ##End Else exclusion criteria NOT Met
  dec = "reject"   
 }                                            ##End Stepwise Regression "While" 
 use.r = regFunc(.regDep,useInd[useVarTF.r],.varTypes,.fullData,.N)
 .ind = rep(FALSE,length(.fullData))
 .ind[useInd[useVarTF.r]] <- TRUE
 return(list(useDep=.regDep,useInd=useInd[useVarTF.r],useReg=use.r$reg,bic=use.r$bic))
}
stepReg <- cmpfun(.stepReg)
##----------------------
##  search is used to locate the independent variables used in regression
##  For distinguishing between variables S and var in R
.search <- function(regString.vec){
who=NULL  
  for(i in 1:length(regString.vec)){
    p1 <- regexpr("\\[\\[",regString.vec[i])
    p2 <- regexpr("\\]\\]",regString.vec[i])
   .tmp <- substr(regString.vec[i],p1,p2)
    who <- c(who,as.numeric(.tmp))
  }
return(who)
}
search <- cmpfun(.search)
##----------------------
##----------------------
## U, W functions
##----------------------
##----------------------
##----------------------
##-------------------------------
.norLik.U <- function(.useDepVars,.useRoles.S,.varTypes,.fullData,.N,.struct){
 nN = sum(.varTypes[.useDepVars]=="N")
 nD = sum(.useRoles.S=="R",na.rm=TRUE)
 
 yN <- muN <- rsN <- matrix(NA,.N,nN)
lik <- rep(NA,.N)
 uReg <- vector("list",nN)
 npar <- ifelse(.struct=="EXX",nN*(1+nD)+1,
		ifelse(.struct=="DXX",nN*(1+nD)+ nN,
						nN*(1+nD)+(nN*(nN+1))/2))
 for(d in 1:nN){
   uReg[[d]] <- regFunc(.useDepVars[d],which(.useRoles.S=="R"),.varTypes,.fullData,.N)
   yN[,d] <- .fullData[[.useDepVars[d]]]
   muN[,d] <- uReg[[d]]$reg$fitted.values
   rsN[,d] <- uReg[[d]]$reg$resid
 }
 sigN <- var(rsN)
 if(.struct=="EXX") sigN <- mean(apply(rsN,2,var))*diag(nN)
 if(.struct=="DXX") sigN <- sigN*diag(nN)
 
 for(i in 1:.N) lik[i] <- dmvnorm(yN[i,],muN[i,],matrix(sigN,nN,nN))
 bic <- -2*sum(log(lik))+npar*log(.N)
 return(list(bic=-bic,reg=uReg))
}
norLik.U <- cmpfun(.norLik.U)
##-------------------------------
.catLik.U <- function(.useDepVars,.useRoles.S,.varTypes,.fullData,.N,.struct){
 nC <- sum(.varTypes[.useDepVars]=="C")
 if(nC >1 & substr(.struct,2,2)=="V"){
   llik <- rep(NA,.N)
   bic <- npar <- NA
   jntCat <- .fullData[.useDepVars]
   d = length(.fullData)+1
   .tmpData = .fullData
   .tmpData[[d]] <- newCat(jntCat,.N,nC)$newCat
   uReg <- regFunc(d,which(.useRoles.S=="R"),c(.varTypes,"C"),.tmpData,.N)
   bic = uReg$bic
 }else{ ##.struct="XEX"
   bic <- 0
   uReg <- vector("list",nC)
   for(iC in 1:nC){
   uReg[[iC]] <- regFunc(.useDepVars[iC],which(.useRoles.S=="R"),.varTypes,.fullData,.N)
   bic <- bic + uReg[[iC]]$bic
   }}
 return(list(bic=bic,reg=uReg))
}
catLik.U <- cmpfun(.catLik.U)
##----------------------Role of U------------------------
.uLik <- function(.useRoles,.useRoles.S,.varTypes,.fullData,.N,.uStruct){
  whichVars = which(.useRoles=="U")
  useDepTypes= .varTypes[whichVars]
  nN = sum(useDepTypes=="N")
  nC = sum(useDepTypes=="C")
  u.EXX = u.DXX = u.VXX = u.XEX = u.XVX = NA;
  uRoleN = uRoleC=list()
  if(any(useDepTypes=="N")){                 ## Normal U variables
    for(.struct in .uStruct$N)  assign(paste("u.",.struct,sep=""),
                                       norLik.U(whichVars[.varTypes[whichVars]=="N"],.useRoles.S,.varTypes,.fullData,.N,.struct))
    uRoleN <- list("EXX"=u.EXX,"DXX"=u.DXX,"VXX"=u.VXX)
  }
  if(any(useDepTypes=="C")){                 ## Multinomial U variables
    for(.struct in .uStruct$C)   assign(paste("u.",.struct,sep=""),
                                        catLik.U(whichVars[.varTypes[whichVars]=="C"],.useRoles.S,.varTypes,.fullData,.N,.struct))
    uRoleC <- list("XEX"=u.XEX,"XVX"=u.XVX)
  }
  return(uRole=list(uRoleN=uRoleN,uRoleC=uRoleC))
}
uLik <- cmpfun(.uLik)
##----------------------

##----------------------
.norLik.W <- function(.useVars,.varTypes,.data,.fullData,.N,.struct){
 ##Prep
  D <- dim(.data$norm)[2]
 ##Run
   sig <- as.matrix(var(.data$norm))
   mu  <- apply(.data$norm,2,mean)
   if(substr(.struct,1,1)=="E") sig <- diag(D)*sig   
   if(substr(.struct,1,1)=="D") sig <- diag(D)*mean(diag(sig))
   d <- dmvnorm(.data$norm,mu,sig)
  
   if(substr(.struct,1,1)=="E" | D ==1) use.mN  <- "EII"
   if(substr(.struct,1,1)=="D") 	      use.mN  <- "EEI"
   if(substr(.struct,1,1)=="V") 	      use.mN  <- "EEE"
   llik <- sum(log(d))
   bic <- -2*llik + useDF.I(use.mN,D,1,1,1)*log(.N)
   return(list(loglik=llik, bic=-bic,par=list(mu=mu,sig=sig)))
}
norLik.W <- cmpfun(.norLik.W)
##----------------------
##Could use Rmixmod with G=1 to use different structures
## Right now is equivalent to the most 'flexible' model for discrete data
.catLik.W <- function(.useVars,.whichVars,.varTypes,.data,.fullData,.N,.struct){
 nC <- sum(.varTypes[.useVars]=="C")
 if(nC> 1 & substr(.struct,2,2)=="V"){
   wReg <- multinom(.data$catData$Repr[[1]]~1,data=.data$catData$Repr)
   llik <- logLik(wReg)
   bic <- BIC(wReg)
 }else{ ##.struct="xEx"
   wReg <- vector("list",nC)
   bic <- llik <- 0
   for(iC in 1:nC){
     wReg[[iC]] <- multinom(.data$catData$Orig[[iC]]~1,data=.data$catData$Orig)
     llik <- llik + logLik(wReg[[iC]])[1]
     bic <- bic + BIC(wReg[[iC]])[1]
 } }
return(list(loglik=llik,bic=-bic))
}
catLik.W <- cmpfun(.catLik.W)
##----------------------Role of W-----------------------
.wLik <- function(.useRoles,.varTypes,.fullData,.N,.wStruct){
  .whichVars = which(.useRoles=="W")
  useTypes=.varTypes[.whichVars]
  w.EXX = w.DXX = w.VXX = w.XEX = w.XVX = NA;
  wRoleN = wRoleC = list()
  .data = makeData(NULL,.whichVars,"Excl",.varTypes,.fullData)
  
  if(any(useTypes=="N")){          ## Normal W variables  
    for(w in .wStruct$N) assign(paste("w.",w,sep=""),
                                norLik.W(.whichVars[.varTypes[.whichVars]=="N"],.varTypes,.data,.fullData,.N,w))   
    wRoleN <- list("EXX"=w.EXX,"DXX"=w.DXX,"VXX"=w.VXX)
  }
  if(any(useTypes=="C")){          ## Multinomial W variables  
    for(w in .wStruct$C)  assign(paste("w.",w,sep=""),
                                 catLik.W(.whichVars[.varTypes[.whichVars]=="C"],.whichVars,.varTypes,.data,.fullData,.N,w))
    wRoleC <- list("XEX"=w.XEX,"XVX"=w.XVX)
  }
  return(wRole=list(wRoleN=wRoleN,wRoleC=wRoleC))
}
wLik <- cmpfun(.wLik)
##-------------------------------
.makeRoles <- function(.useRoles,.varTypes,.baseStruct){
 .uNom <- sum(.varTypes[.useRoles=="U"]=="N")
 .uCat <- min(sum(.varTypes[.useRoles=="U"]=="C"),2)
 .wNom <- sum(.varTypes[.useRoles=="W"]=="N")
 .wCat <- min(sum(.varTypes[.useRoles=="W"]=="C"),2)

 if(.uNom >1){ 
  .un = 2:4
 }else{ .un <- ifelse(.uNom==0,1,2) }
 if(.uCat >1){ 
  .uc = 2:3
 }else{ .uc <- ifelse(.uCat==0,1,2) }
 uStruct <- list()
 uStruct$N <- .baseStruct$N[.un]
 uStruct$C <- .baseStruct$C[.uc]
 if(.wNom >1){ 
  .wn = 2:4
 }else{ .wn <- ifelse(.wNom==0,1,2) }
 if(.wCat>1){ 
  .wc = 2:3
 }else{ .wc <- ifelse(.wCat==0,1,2) }
 wStruct <- list()
 wStruct$N <- .baseStruct$N[.wn]
 wStruct$C <- .baseStruct$C[.wc]
return(list(uStruct=uStruct,wStruct=wStruct))
}
makeRoles <- cmpfun(.makeRoles)
##-------------------------------
##-------------------------------
## STEPWISE FUNCTIONS
##-------------------------------
##-------------------------------
.makeVar.S <- function(i,.varRoles,.step){
  if(.step=="Excl"){
   propVar = which(.varRoles=="S")[i]
   remnVar = which(.varRoles=="S")[-i]
  }else{
   propVar = which(.varRoles=="Sc")[i]
   remnVar = which(.varRoles=="S")
 }
  return(list(propVar=propVar,remnVar=remnVar))
}
makeVar.S <- cmpfun(.makeVar.S)
##---------------------------
.makeRS = function(makeVar,var,i){
 tmp = rep(NA,length(makeVar))
 tmp[var$remnVar] = "R"
return(tmp)
}
makeRS <- cmpfun(.makeRS)
##-----------------------------------------------------
##-----------------------------------------------------
## For calculation of FULL BIC as in Maugis paper
##-----------------------------------------------------
##-----------------------------------------------------
## From S, Sc, determine S, R, U, W
.makeSRUW <- function(.varRoles,.fullData,.N,.varTypes){
  nV = length(.varTypes)
  if(sum(.varRoles=="S")<nV){
  	useRoles   = .varRoles
  	useRoles.S = rep(list(rep(NA,length(.varRoles))),sum(.varRoles=="Sc"))
  	tmpS <- which(.varRoles=="S")
  	tmpU <- which(.varRoles=="Sc")
  	tmpR = vector("list",length(tmpU))
  	for(i in 1:length(tmpU)){
    		tmpR[[i]] = whichReg(tmpU[i],tmpS,.varTypes,.fullData,.N)
    		if(length(tmpR[[i]]$useInd)>0){ 
      	 useRoles[tmpR[[i]]$useDep] <- "U"
      	 useRoles.S[[i]][tmpR[[i]]$useInd] <- "R"
      	 if(length(tmpR[[i]]$useInd) < length(tmpS)){
          S = tmpS[-tmpR[[i]]$useInd]
        	useRoles.S[[i]][S] <- "S"
      	 }}else{
      	 useRoles[tmpR[[i]]$useDep] <- "W"
      	 useRoles.S[[i]][tmpS] <- "S"
  	} 	}
  	useR <- rep(NA,nV)
  	for(i in 1:length(tmpU)) useR[which(useRoles.S[[i]]=="R")] <- "R"
  }else{
	 useRoles = .varRoles
	 useR = rep(NA,nV)
  }   
return(list(useRoles=useRoles,useR=useR))
}
makeSRUW <- cmpfun(.makeSRUW)
##---------------------
## From S, Sc, determine S, R, U, W
.makeSW <- function(.varRoles,.fullData,.N,.varTypes){
  nV = length(.varTypes)
  useRoles = .varRoles
  useRoles[.varRoles=="Sc"] = "W" 
  useR = rep(NA,nV)   
return(list(useRoles=useRoles,useR=useR))
}
makeSW <- cmpfun(.makeSW)
##---------------------
##Calculate Overall BIC: Sbic, U~Rbic, Wbic
.fullBIC <- function(.allRoles,.SRRoles,.g1,.g2,.mM,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc){
  nU = sum(.allRoles=="U")
  nW = sum(.allRoles=="W")
  allBIC <- matrix(NA,5,3)   ##allBIC[1,1] = BIC.s, allBIC[1:3,2,3] = Normal BIC.U,BIC.w, allBIC[4:5,2:3] = Categorical BIC.U,BIC.W
  colnames(allBIC) <- c("S","U","W")
  tmpS  <-  whichClust(NULL,which(.allRoles=="S"),"Excl",
                       .fullData,.g1,.g2,.varTypes,.mM,.which.pri,.kp,.tol,.max.iter,.mc)      #?

  if(is.na(tmpS$clustFlag)){ ##Valid clustering, proceed with full BIC calculation
   allBIC[1,"S"] <- tmpS$bic ##ClusterBIC
   uw <- makeRoles(.allRoles,.varTypes,baseStruct); 
   if("U" %in% .allRoles){   ##If variables in U, grab BIC from Normal U vars, BIC from Categorical U Vars
    tmpU <- uLik(.allRoles,.SRRoles,.varTypes,.fullData,.N,uw$uStruct)
    for(x in c(na.omit(unlist(uw$uStruct$N))))
      allBIC[which(na.omit(unlist(uw$uStruct$N))==x),which(colnames(allBIC)=="U")] <-
        eval(parse(text=paste("tmpU$uRoleN$",x,"$bic",sep="")))
    for(x in c(na.omit(unlist(uw$uStruct$C))))
      allBIC[3+which(na.omit(unlist(uw$uStruct$C))==x),which(colnames(allBIC)=="U")] <-
        eval(parse(text=paste("tmpU$uRoleC$",x,"$bic",sep="")))
   }
   if("W" %in% .allRoles){  #If variables in U, grab BIC from Normal U vars, BIC from Categorical U Vars     ?  W
    tmpW <- wLik(.allRoles,.varTypes,.fullData,.N,uw$wStruct)
    for(y in na.omit(unlist(uw$wStruct$N)))   
           allBIC[which(na.omit(unlist(uw$wStruct$N))==y),which(colnames(allBIC)=="W")] <-
             eval(parse(text=paste("tmpW$wRoleN$",y,"$bic",sep="")))   
    for(y in na.omit(unlist(uw$wStruct$C)))   
           allBIC[3+which(na.omit(unlist(uw$wStruct$C))==y),which(colnames(allBIC)=="W")] <-
             eval(parse(text=paste("tmpW$wRoleC$",y,"$bic",sep="")))  
   }
   rownames(allBIC) <- na.omit(c(unlist(baseStruct$N),unlist(baseStruct$C)))
  
   if("U" %in% .allRoles){
     xN <- uw$uStruct$N; xC <- uw$uStruct$C
     if(c(!all(is.na(xN)) & all(is.na(xC))) | (c(all(is.na(xN)) & !all(is.na(xC))))){ ##Unmixed U: All Norm or All Cat
     rowBIC <- na.omit(c(unlist(uw$uStruct)))
     tmpUBIC <- matrix(NA,length(rowBIC),1)
     dimnames(tmpUBIC) <- list(rowBIC,"Nul W")
     for(x in which(!is.na(allBIC[,"U"])))   tmpUBIC[rownames(allBIC)[x],1] <- sum(allBIC[1,"S"],allBIC[x,"U"])
    }else{                                   ##Mixed U: Norm and Cat, Must add components from each "structure" (EEX = EXX+XEX, etc)
     rowBIC <- paste(rep(substr(xN,1,1),each=length(xC)),substr(xC,2,2),sep="")
     tmpUBIC <- matrix(NA,length(rowBIC),1)
     dimnames(tmpUBIC) <- list(rowBIC,"Nul W") 
     for(.xN in xN){
       for(.xC in xC) tmpUBIC[paste(substr(.xN,1,1),substr(.xC,2,2),sep=""),1] <- sum(allBIC[1,"S"],allBIC[.xN,"U"],allBIC[.xC,"U"]) 
    }}}
    if("W" %in% .allRoles){
    yN <- uw$wStruct$N; yC <- uw$wStruct$C 
    if(c(!all(is.na(yN)) & all(is.na(yC))) | (c(all(is.na(yN)) & !all(is.na(yC))))){ ##Unmixed W: All Norm or All Cat
     colBIC <- na.omit(c(unlist(uw$wStruct)))
     tmpWBIC <- matrix(NA,1,length(colBIC))
     dimnames(tmpWBIC) <- list("Nul U",colBIC)
     for(y in which(!is.na(allBIC[,"W"])))   tmpWBIC[1,rownames(allBIC)[y]] <- sum(allBIC[1,"S"],allBIC[y,"W"])
    }else{                                   ##Mixed W: Norm and Cat, Must add components from each "structure" (EEX = EXX+XEX, etc)
     colBIC <- paste(rep(substr(yN,1,1),each=length(yC)),substr(yC,2,2),sep="")  
     tmpWBIC <- matrix(NA,1,length(colBIC))
     dimnames(tmpWBIC) <- list("Nul U",colBIC) 
     for(.yN in yN){
       for(.yC in yC) tmpWBIC[1,paste(substr(.yN,1,1),substr(.yC,2,2),sep="")] <- sum(allBIC[1,"S"],allBIC[.yN,"W"],allBIC[.yC,"W"])
    }}} 
   if(("U" %in% .allRoles) & !("W" %in% .allRoles)) resBIC <- tmpUBIC ##Only U, No W
   if(("W" %in% .allRoles) & !("U" %in% .allRoles)) resBIC <- tmpWBIC ##Only W, No U
   if(("U" %in% .allRoles) & ("W" %in% .allRoles)){                   ##Both U, W
     resBIC <- matrix(NA,nrow(tmpUBIC),ncol(tmpWBIC))
     dimnames(resBIC) <- list(rownames(tmpUBIC),colnames(tmpWBIC))
     for(x in 1:nrow(resBIC)){
      for(y in 1:ncol(resBIC)) resBIC[x,y] <- tmpUBIC[x,1] + tmpWBIC[1,y]
   }}
   if(!("U" %in% .allRoles) & !("W" %in% .allRoles)){                 ##No U, W
     resBIC <- matrix(NA,1,1)
     dimnames(resBIC) <- list("No U","No W")
     resBIC["No U","No W"] <- allBIC[1,"S"]
   }  
  return(list(clust=tmpS,resBIC=resBIC,allBIC=allBIC))
 }else{ ##Clustering Error
  cat(paste("Cannot compute clustering on S:",.mM,", ",.g1,.g2," \n",sep="") )
  return(list(clust=list(bic=NA),resBIC=matrix(NA,3,3),allBIC=matrix(NA,1,1)))
}}
fullBIC <- cmpfun(.fullBIC)
##---------------------------------------------------------------
##-------------------------------
.stepExcl <- function(.clustVar,.g1,.g2,.modelName,.nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc){
  uV  <- rep("Sc",.nV);
  uV[.clustVar] <- "S"
  ## Cluster on current .clustVar variables
  tS <- whichClust(NULL,.clustVar,"Excl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
  if(!is.null(tS$model)){
    res <- matrix(NA,sum(uV=="S"),9)
    for(i in 1:(sum(uV=="S"))){
    var = tc = trg = NULL
    var <- makeVar.S(i,uV,"Excl")                                   ## Propose variable from exclusion, remnVar = .clustVar-propVar
    tc  <- whichClust(var$propVar,var$remnVar,"Excl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
    if(!is.null(tc$model)){
      if(.varTypes[var$propVar]=="C" & substr(tS$model,6,6)=="J"){       ##If PROPOSED cluster tc was calculated on catData$Repr, use catData$Repr for reg
        useC <- .clustVar[.varTypes[.clustVar]=="C"]
        useN <- .clustVar[.varTypes[.clustVar]=="N"]
        newVar <- length(.fullData)+1
        mD <- makeData(var$propVar,var$remnVar,"Incl",.varTypes,.fullData)
        tmpData <- cbind(.fullData,mD$catData$Repr)
        trg1 <- whichReg(newVar,useN,c(.varTypes,"C"),tmpData,.N)
        if(substr(tc$model,6,6)=="J"){                                   ##If PREVIOUS cluster tS was calculated on catData$Repr, use catData$Repr for reg
          newVar <- length(.fullData)+1
          mD2 <- makeData(NULL,var$remnVar,"Excl",.varTypes,.fullData)
          tmpData2 <- cbind(.fullData,mD2$catData$Repr)
          #trg2 = step.regFunc(newVar,useN,c(.varTypes,"C"),tmpData2,.N,tc$model)$bic
          trg2 <- whichReg(newVar,useN,c(.varTypes,"C"),tmpData2,.N)$bic
        }else{                                                          ##Else PREVIOUS cluster tS was calculated on catData$Orig, use catData$Orig for reg
          useC <- var$remnVar[var$remnVar %in% which(.varTypes=="C")] #.varTypes[var$remnVar].clustVar]=="C"]
          tmptrg2 <- rep(NA,length(useC))
          for(uC in useC) tmptrg2[uC] <- whichReg(uC,useN,.varTypes,.fullData,.N)$bic#step.regFunc(uC,useN,.varTypes,.fullData,.N,tc$model)$bic
          trg2 <- sum(tmptrg2,na.rm=TRUE)  
        }
        trg = list(bic=trg1$bic-trg2)
        }else{ trg  <- whichReg(var$propVar,var$remnVar,.varTypes,.fullData,.N) }    #trg  = step.regFunc(var$propVar,var$remnVar,.varTypes,.fullData,.N,tS$model) }     
       res[i,] <- c(paste(c(var$remnVar),collapse=","),round(tS$bic),round(tc$bic),round(trg$bic),round(tS$bic - {tc$bic+trg$bic}),tS$model,tc$model,"Excl","-")
    }else{
      res[i,] <- c(paste(c(var$remnVar),collapse=","),round(tS$bic),round(tc$bic),NA,NA,tS$model,NA,"Excl","-")
    }}}else{
        res <- matrix(NA,1,9)
        res[1,] <- c(paste(.clustVar,collapse=","),NA,NA,NA,NA,NA,NA,"Excl","-")
    }
  colnames(res)<- c("ClustVar","ClustTot","ClustLess","RegFull","Tot-(Less+Reg)","ClustTot.Model","ClustLess.Model","Step","Decision")
  tstBIC       <- as.numeric(res[,"Tot-(Less+Reg)"])
  whichExcl    <- ifelse(sum(tstBIC<excl.Crit,na.rm=TRUE)>0,which(tstBIC <excl.Crit & tstBIC ==  min(tstBIC[tstBIC<excl.Crit],na.rm=TRUE))[1],NA)
  if(!is.na(whichExcl)){
    whichRemnVar <- as.numeric(strsplit(res[whichExcl,"ClustVar"],",")$ClustVar )
    acc          <- TRUE
  }else{
    whichRemnVar <- .clustVar
    acc          <- FALSE
  } 
  ##If exclusion results in all normal values being excluded, reject step
  if(sum(.varTypes[whichRemnVar]=="N")==0){
    whichRemnVar <- .clustVar
    acc          <- FALSE
  }
  res[whichExcl,"Decision"] <- ifelse(acc,"Accept","Reject")
  return(list(res=res,rowExcl=whichExcl,remnVar=whichRemnVar,accept=acc))
}
stepExcl <- cmpfun(.stepExcl)
##-------------------------------
.stepIncl <- function(.clustVar,.g1,.g2,.modelName,.nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc){
  ##Include, .prevStep=NULL or .prevStep = bic
  uV  <- rep("Sc",.nV)
  uV[.clustVar] <- "S"
  tc  <- whichClust(NULL,.clustVar,"Excl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
  if(!is.null(tc$model)){
    res <- matrix(NA,sum(uV=="Sc"),9)
    for(i in 1:(sum(uV=="Sc"))){
    var = tS = trg = NULL
    var <- makeVar.S(i,uV,"Incl")
    tS  <- whichClust(var$propVar,var$remnVar,"Incl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
    if(!is.null(tS$model)){
      if(.varTypes[var$propVar]=="C" & substr(tS$model,6,6)=="J"){       ##If PROPOSED cluster tS was calculated on catData$Repr, use catData$Repr for reg
        useC <- .clustVar[.varTypes[.clustVar]=="C"]
        useN <- .clustVar[.varTypes[.clustVar]=="N"]
        newVar  <- length(.fullData)+1
        mD      <- makeData(var$propVar,var$remnVar,"Incl",.varTypes,.fullData)
        tmpData <- cbind(.fullData,mD$catData$Repr)
        trg1 <- whichReg(newVar,useN,c(.varTypes,"C"),tmpData,.N)
        if(substr(tc$model,6,6)=="J"){                                   ##If PREVIOUS cluster tc was calculated on catData$Repr, use catData$Repr for reg
          newVar <- length(.fullData)+1
          mD2 <- makeData(NULL,var$remnVar,"Excl",.varTypes,.fullData)
          tmpData2 <- cbind(.fullData,mD2$catData$Repr)
          trg2 <- whichReg(newVar,useN,c(.varTypes,"C"),tmpData2,.N)$bic
        }else{                                                           ##Else PREVIOUS cluster tc was calculated on catData$Repr, use catData$Repr for reg
          tmptrg2 = rep(NA,length(.varTypes))
          for(uC in useC) tmptrg2[uC] <- whichReg(uC,useN,.varTypes,.fullData,.N)$bic
          trg2 <- sum(tmptrg2,na.rm=TRUE)  
        }
       trg = list(bic=trg1$bic-trg2)
      }else{ trg  <- whichReg(var$propVar,var$remnVar,.varTypes,.fullData,.N)}
      if(substr(tc$model,5,5)=="I" & substr(tS$model,5,5)=="R"){
        tmp.tc <- tc$BIC[,substr(colnames(tc$BIC),5,5)=="R"] 
        tmp.ind <- which(max(tmp.tc,na.rm=TRUE),arr.ind=TRUE)
        res[i,] <- c(paste(c(var$remnVar,var$propVar),collapse=","),round(tS$bic),round(tmp.tc[tmp.ind[1],tmp.ind[2]]),
                    round(trg$bic),round(tS$bic - {tmp.tc+trg$bic}),tS$model,colnames(tmp.tc)[tmp.ind[2]],"Incl","-") 
      }else{  res[i,] <- c(paste(c(var$remnVar,var$propVar),collapse=","),round(tS$bic),round(tc$bic),round(trg$bic),round(tS$bic-{tc$bic+trg$bic}),tS$model,tc$model,"Incl","-") }
     }else{
      res[i,] <- c(paste(c(var$remnVar,var$propVar),collapse=","),round(tS$bic),round(tc$bic),NA,NA,NA,tc$model,"Incl","-")
    }}}else{
       res <- matrix(NA,1,9)
       res[1,] <- c(paste(.clustVar,collapse=","),NA,NA,NA,NA,NA,NA,"Incl","-")
    }
  colnames(res)<- c("ClustVar","ClustTot","ClustLess","RegFull","Tot-(Less+Reg)","ClustTot.Model","ClustLess.Model","Step","Decision")
  tstBIC       <- as.numeric(res[,"Tot-(Less+Reg)"])
  whichIncl    <- ifelse(sum(tstBIC>incl.Crit,na.rm=TRUE)>0,which(tstBIC >incl.Crit & tstBIC ==  max(tstBIC[tstBIC>incl.Crit],na.rm=TRUE))[1],NA)
  if(!is.na(whichIncl)){
    whichRemnVar <- as.numeric(strsplit(res[whichIncl,"ClustVar"],",")$ClustVar )
    acc <- TRUE
  }else{
    whichRemnVar <- .clustVar
    acc <- FALSE
  }
  res[whichIncl,"Decision"] <- ifelse(acc,"Accept","Reject")
  return(list(res=res,rowIncl=whichIncl,remnVar=whichRemnVar,accept=acc))
}
stepIncl <- cmpfun(.stepIncl)
##---------------------------------------------------------------
## Stepwise implement
##---------------------------------------------------------------
##---------------------------------------------------------------
#.startVar=1:7; .g1=.g2=2;.modelName="EII";
.stepSelect.Maug <- function(.startVar,.g1,.g2,.modelName,.fullData,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc,.intermed){
  nV = dim(.fullData)[2]; .N = dim(.fullData)[1]
  ##-----
  ##Determine S,Sc
  iE = iI = vector("list",nV*2)
  i  = 1
  iE[[i]] <- stepExcl(.startVar,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc);
  iI[[i]] <- list(remnVar = c(iE[[1]]$remnVar),accept=FALSE)
  iter =  ifelse(iE[[1]]$accept,TRUE,FALSE)
  s.flag1 <- s.flag2 <- ifelse(iter,TRUE,FALSE)
  while(iter & i < nV*2){
    i <- i+1
    iE[[i]] <- stepExcl(iI[[i-1]]$remnVar,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc);
    if(iE[[i]]$accept | iI[[i-1]]$accept){ iI[[i]] <- stepIncl(iE[[i]]$remnVar,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc);
    }else{                                 iI[[i]] <- list(remnVar=iE[[i]]$remnVar,accept=FALSE)}
    if(iE[[i]]$accept){  
	   if(iI[[i]]$accept==FALSE & length(iE[[i]]$remnVar)==1)   s.flag1 = FALSE #Reject Include, Only 1 var in S
	   if(length(iI[[i]]$remnVar)==1)                   	      s.flag2 = FALSE #Only 1 var in S
	   if(length(iI[[i]]$remnVar)==nV)                   	      s.flag2 = FALSE #Back to all in S
    }else{##Nothing to Exclude	
      if(iI[[i-1]]$accept==FALSE | iI[[i]]$accept==FALSE)     s.flag1 = FALSE #Reject consec Excl/Incl | Incl/Excl
    }
    iter <- ifelse(s.flag1==FALSE | s.flag2 == FALSE,FALSE,TRUE)
   }##End "S" loop
  iE <- iE[1:i]; iI <- iI[1:i]
  varRoles = rep("Sc",nV)
  varRoles[iE[[i]]$remnVar] = "S"
  cat("End S, Sc \n")
  if(!iter){
    tmpSW <- tmpSRUW <- list(useRoles=rep("S",nV),useR=rep(NA,nV))
  }else{
  ##-----
  ##Determine SRUW from S,Sc
  tmpSW <- makeSW(varRoles,.fullData,.N,.varTypes)
  tmpSRUW <- makeSRUW(varRoles,.fullData,.N,.varTypes)
  }
  cat("End S,R,U,W \n")
  ##-----
  ##Calculate Overall BIC
  calcBICSW <- fullBIC(tmpSW$useRoles,tmpSW$useR,.g1,.g2,.modelName,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)
  calcBICSRUW <- fullBIC(tmpSRUW$useRoles,tmpSRUW$useR,.g1,.g2,.modelName,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)
  cat("End Calc BIC \n")
  
  tmpRoleSW   <- tmpSW$useRoles
  tmpRoleSRUW <- tmpSRUW$useRoles
  tmpRoleSRUW[which(tmpSRUW$useR=="R")] <- "R"
  if(!.intermed){
	 return(list(SW=list(fBIC=calcBICSW,allRoles=tmpRoleSW),SRUW=list(fBIC=calcBICSRUW,allRoles=tmpRoleSRUW)))
  }else{ 
    allSteps = NULL
    for(s in 1:max(length(iE),length(iI))){
      if(!is.null(iI[[s]])) allSteps = rbind(allSteps,iE[[s]]$res)
      if(!is.null(iE[[s]])) allSteps = rbind(allSteps,iI[[s]]$res)
    }      
    return(list(SW=list(fBIC=calcBICSW,allRoles=tmpRoleSW),SRUW=list(fBIC=calcBICSRUW,allRoles=tmpRoleSRUW),
					                    steps=allSteps))
  }
}

stepSelect.Maug <- cmpfun(.stepSelect.Maug)
##--------------
.uniStep <- function(.startVar,.g1,.g2,.modelName,.nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc){
  noClust <- uniClust <- numeric()
  useV = length(.startVar)
  for(i in 1:useV){
	 noClust[i] <- whichClust(NULL,.startVar[i],"Excl",.fullData,1,1,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)$bic
	 if(is.na(noClust[i])){noClust[i] = 0} else{noClust[i]} 
	 }       # ? solved I add
  for(i in 1:useV){
    uniClust[i] <- whichClust(NULL,.startVar[i],"Excl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)$bic   
  if(is.na(uniClust[i])){uniClust[i] = 0} else{uniClust[i]} 
    }      # ? solved I add
  #res <- cbind(paste(.startVar),round(uniClust),round(noClust),rep(NA,useV),round(uniClust-noClust),rep("-",useV),rep("-",useV),"Uni-Incl","-")
  # if uniClust = NA, then uniClust-noClust = NA, then tstBIC = NA
  res <- cbind(paste(.startVar),round(uniClust),round(noClust),rep(NA,useV),round(uniClust-noClust),rep("-",useV),rep("-",useV),"Uni-Incl","-")
  colnames(res)= c("ClustVar","ClustTot","ClustLess","RegFull","Tot-(Less+Reg)","ClustTot.Model","ClustLess.Model","Step","Decision") 
  tstBIC <- as.numeric(res[,"Tot-(Less+Reg)"])
  whichIncl <- which(tstBIC==max(tstBIC,na.rm=TRUE))[1]
  res[whichIncl,"Decision"] <- "Accept"
  return(list(res=res,remnVar=whichIncl,accept=TRUE))
}
uniStep <- cmpfun(.uniStep)
##-------------------------------
.biStep <- function(.clustVar,.g1,.g2,.modelName,.nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc){
  uV  = rep("Sc",.nV)
  uV[.clustVar] = "S"
  res = matrix(NA,sum(uV=="Sc"),9)
  tc <- whichClust(NULL,.clustVar,"Excl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
  for(i in 1:(sum(uV=="Sc"))){
    var = tS = tr1 = tr2 = NULL
    var = makeVar.S(i,uV,"Incl")
    tS  = whichClust(var$propVar,var$remnVar,"Incl",.fullData,.g1,.g2,.varTypes,.modelName,.which.pri,.kp,.tol,.max.iter,.mc)
    if(!is.null(tS$model)){
      #trg  = step.regFunc(var$propVar,var$remnVar,.varTypes,.fullData,.N,tS$model)
      trg  = whichReg(var$propVar,var$remnVar,.varTypes,.fullData,.N)
      res[i,] = c(paste(c(var$remnVar,var$propVar),collapse=","),round(tS$bic),round(tc$bic),round(trg$bic),round(tS$bic - {tc$bic+trg$bic}),tc$model,tS$model,"Bi-Incl","-")
    }else{
      if(is.null( tc$model)){ res[i,] = c(paste(c(var$remnVar,var$propVar),collapse=","),NA,round(tc$bic),NA,NA,NA,NA,"Bi-Incl","-")}    #begin to change
      else{
      res[i,] = c(paste(c(var$remnVar,var$propVar),collapse=","),NA,round(tc$bic),NA,NA,tc$model,NA,"Bi-Incl","-")}                       #end
      #orginal:  res[i,] = c(paste(c(var$remnVar,var$propVar),collapse=","),NA,round(tc$bic),NA,NA,tc$model,NA,"Bi-Incl","-")
    }
  }
  colnames(res)= c("ClustVar","ClustTot","ClustLess","RegFull","Tot-(Less+Reg)","ClustTot.Model","ClustLess.Model","Step","Decision")  
  tstBIC       = as.numeric(res[,"Tot-(Less+Reg)"])
  #tstBIC    = as.numeric(res[which(.varTypes[-.clustVar]=="N"),"Tot-(Less+Reg)"])
  whichIncl    = which(tstBIC ==  max(tstBIC,na.rm=TRUE))[1]
  
  if(!is.na(whichIncl)){
      whichRemnVar = as.numeric(strsplit(res[whichIncl,"ClustVar"],",")$ClustVar )
      acc = TRUE
    }else{
      whichRemnVar = .clustVar
      acc = FALSE
    }
  res[whichIncl,"Decision"] <- ifelse(acc,"Accept","-")
  return(list(res=res,rowIncl=whichIncl,remnVar=whichRemnVar,accept=acc))
}
biStep <- cmpfun(.biStep)
##--------------
.switchStep <- function(.iI,.g1,.g2,.modelName,.nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc){
  res <- matrix(NA,.nV,9)
  colnames(res)= c("ClustVar","ClustTot","ClustLess","RegFull","Tot-(Less+Reg)","ClustTot.Model","ClustLess.Model","Step","Decision")  
  res[,"Step"] <- "Switch"
  res[,"Decision"] <- "-"
  keep <- .iI[[2]]$remnVar[-which(.iI[[2]]$remnVar==.iI[[1]]$remnVar)]
  for(m in which(.varTypes=="C")){
    cl <- whichClust(NULL,c(keep,m),"Excl",.fullData,.g1,.g2,.varTypes,.modelName,
                     .which.pri,.kp,.tol,.max.iter,.mc)
    r1 <- whichReg(.iI[[1]]$remnVar,c(keep,m),.varTypes,.fullData,.N)$bic
    res[m,"ClustVar"] <- paste(c(keep,m),collapse=",")
    res[m,"ClustTot"] <- round(cl$bic + r1)
    res[m,"ClustLess"]<- as.numeric(.iI[[2]]$res[.iI[[2]]$rowIncl,"ClustTot"])
    res[m,"RegFull"]  <- round(whichReg(m,c(.iI[[1]]$remnVar,keep),.varTypes,.fullData,.N)$bic)
    res[m,"Tot-(Less+Reg)"] <- as.numeric(res[m,"ClustTot"]) - {
      as.numeric(res[m,"ClustLess"])+ as.numeric(res[m,"RegFull"])}
    res[m,"ClustTot.Model"] <- cl$model
    res[m,"ClustLess.Model"] <- .iI[[2]]$res[.iI[[2]]$rowIncl,"ClustTot.Model"]
  }
  switchVar <- ifelse(sum(as.numeric(res[,"Tot-(Less+Reg)"])>0,na.rm=TRUE)>0,
                      which.max(as.numeric(res[,"Tot-(Less+Reg)"])),NA)
  if(!(is.na(switchVar))){
    res[switchVar,"Decision"] <- "Accept"
    tmp <- list(res=na.omit(res),rowSwitch=switchVar,
                remnVar=as.numeric(unlist(strsplit(res[switchVar,"ClustVar"],","))),accept=TRUE)
  }else{
    tmp <- list(res=na.omit(res),rowSwitch=switchVar,remnVar=.iI[[2]]$remnVar,accept=FALSE)       
  }
  return(tmp)
}
switchStep <- cmpfun(.switchStep)
##--------------
.stepSelect.Dean <- function(.startVar,.g1,.g2,.modelName,.fullData,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc,.intermed){
  nV =length(.startVar)
  .N = dim(.fullData)[1]
  ##-----
  ##Determine S,S
  iE = iI = vector("list",nV*2)
  i  = 1
  iI[[i]] <- uniStep(.startVar,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)     #? solved
  iE[[i]] <- list(remnVar = c(iI[[1]]$remnVar),accept=FALSE)
  if(!is.na(iI[[i]]$remnVar)){
   i <- i+1
   iI[[i]] <- biStep(iE[[i-1]]$remnVar,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)  #? solved
   iE[[i]] <- list(remnVar = c(iI[[i]]$remnVar),accept=FALSE)
   iter <-   TRUE
   s.flag1 <- s.flag2 <- ifelse(iter,TRUE,FALSE)
  }else{iter <- FALSE}
  ##Switch
  if(iter & all(.varTypes[iI[[i]]$remnVar]=="N") & any(.varTypes =="C")){
   i <- i+1
   iI[[i]] <- switchStep(iI,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)  
   iE[[i]] <- list(remnVar = c(iI[[i]]$remnVar),accept=FALSE) 
  }
   while(iter){
    i <- i+1
    iI[[i]] <- stepIncl(iE[[i-1]]$remnVar,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)
    if(length(iI[[i]]$remnVar)>1 & (iI[[i]]$accept | iE[[i-1]]$accept)){
	   iE[[i]] <- stepExcl(iI[[i]]$remnVar,.g1,.g2,.modelName,nV,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)
    }else{ iE[[i]] <- list(remnVar=iI[[i]]$remnVar,accept=FALSE) }
    if(iI[[i]]$accept){  
	   if(iE[[i]]$accept==FALSE & length(iE[[i]]$remnVar)==nV) s.flag1 = FALSE #Reject Exclude, All var in S
     if(length(iE[[i]]$remnVar)==nV)                   	     s.flag2 = FALSE #All var in S   
     if(length(iI[[i]]$remnVar)==1)                          s.flag2 = FALSE #Back to 1 var in S
    }else{##Nothing to Exclude	
     if(iE[[i-1]]$accept==FALSE | iE[[i]]$accept==FALSE)     s.flag1 = FALSE #Reject consec Excl/Incl | Incl/Excl
     if(length(iI[[i-1]]$remnVar)==1)                        s.flag2 = FALSE
    }
    iter <- ifelse(s.flag1==FALSE | s.flag2 == FALSE,FALSE,TRUE)
  }##End "S" loop
  iI <- iI[1:i]; iE <- iE[1:i]
  varRoles = rep("Sc",nV)
  varRoles[iE[[i]]$remnVar] = "S"
  cat("End S, Sc \n")
 
  ##-----
  ##Determine SRUW from S,Sc
  tmpSW <- makeSW(varRoles,.fullData,.N,.varTypes)
  tmpSRUW <- makeSRUW(varRoles,.fullData,.N,.varTypes)
  cat("End S,R,U,W \n")
  ##-----
  ##Calculate Overall BIC
  calcBICSW <- fullBIC(tmpSW$useRoles,tmpSW$useR,.g1,.g2,.modelName,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)
  calcBICSRUW <- fullBIC(tmpSRUW$useRoles,tmpSRUW$useR,.g1,.g2,.modelName,.fullData,.N,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc)
  cat("End Calc BIC \n")     #?
  
  tmpRoleSW   <- tmpSW$useRoles
  tmpRoleSRUW <- tmpSRUW$useRoles
  tmpRoleSRUW[which(tmpSRUW$useR=="R")] <- "R"
  if(!.intermed){
 	 return(list(SW=list(fBIC=calcBICSW,allRoles=tmpRoleSW),SRUW=list(fBIC=calcBICSRUW,allRoles=tmpRoleSRUW)))
  }else{ 
    allSteps = NULL
    for(i in 1:max(length(iE),length(iI))){
      if(!is.null(iI[[i]])) allSteps = rbind(allSteps,iI[[i]]$res)
      if(!is.null(iE[[i]])) allSteps = rbind(allSteps,iE[[i]]$res)
    }      
    return(list(SW=list(fBIC=calcBICSW,allRoles=tmpRoleSW),SRUW=list(fBIC=calcBICSRUW,allRoles=tmpRoleSRUW),
					steps=allSteps))}#list(iE=iE,iI=iI))) }
}
stepSelect.Dean <- cmpfun(.stepSelect.Dean)
##--------------------
.stepSelect <- function(.startVar,.g1,.g2,.mM,.fullData,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc,.intermed, .meth){
 if(.meth=="Dean"){
	s.Select = stepSelect.Dean(.startVar,.g1,.g2,.mM,.fullData,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc,.intermed)
 }else{ #.meth="Maug"
	s.Select = stepSelect.Maug(.startVar,.g1,.g2,.mM,.fullData,.varTypes,.which.pri,.kp,.tol,.max.iter,.mc,.intermed)
 }
 return(s.Select)
}
stepSelect <- cmpfun(.stepSelect)
