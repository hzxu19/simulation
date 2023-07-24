command=TRUE
if(command){
  if(getwd()!="/u/kevans/Research Files/newVS"){
    if(getwd()=="C:/Documents and Settings/kevans/My Documents"){
      setwd("C://Documents and Settings//kevans//My Documents//My Dropbox//tmp//newVS//")
    }else{
      setwd("C://Users//K80//Documents//My Dropbox//tmp//newVS//")
    }}
  
  source("vsLoad.R")
  source("http://rtm.wustl.edu/code/sendEmail.R")
  source("vsFunction.R")
  
  Args <- commandArgs(TRUE)
  #Args <- c("21","1","1","1","1","1","5")
  ##-------------------------------------------------
  ##--Arguments set by "Args" and/or sourced from vsLoad.R
  runPre <- ifelse(as.numeric(Args[2])==1,TRUE,FALSE)
  runDV  <- ifelse(as.numeric(Args[3])==1,TRUE,FALSE)
  runMV  <- FALSE #ifelse(Args[3]==as.numeric(1),TRUE,FALSE)
  runDF  <- FALSE #ifelse(Args[4]==as.numeric(1),TRUE,FALSE)
  runMF  <- ifelse(as.numeric(Args[4])==1,TRUE,FALSE)
  whichRun  <- ifelse(as.numeric(Args[5])==0,"Null","Step") 
  if(whichRun=="Null"){ whichReg <- regFunc }else{ whichReg <- stepReg }
  source("vsFunction.R")
  R         <- substr(whichRun,1,1)
  
  d.i       <- as.numeric(Args[1])
  tmpD      <- read.csv(file=paste("d",d.i,".csv",sep=""),header=FALSE)        ##data
  d         <- data.frame(tmpD)
  N         <- sum(d[,8]==1)
  resN      <- paste("res",d.i,sep="") #.",Args[3],Args[4],Args[5],sep="")         ##Results file name
  u         <- d[,7]
  G         <- TruthVec[TruthVec$Ind==d.i,"G"]               ##True Number of Clusters (for prelim check)
  useS      <- paste(unlist(TrueS[TrueS$Ind==d.i,2+(1:6)]))
  saveRes   <- FALSE
  ##-------------------------------------------------
  ##--PRESETS----------------------------------------
  .startVar <- 1:6                               ##Starting Point
  if(d.i < 100){
   d[,5]     <- factor(d[,5]); 
   d[,6]     <- factor(d[,6])
   varTypes  <- c("N","N","N","N","C","C")        ##Variable Types
   nN        <- 4
   nC        <- 2
  tryMods    <- tryMods1 <- paste(rep(.mclustModel,each=4),c(".II",".IJ",".RI",".RJ"),sep="") ##Range of Models to consider
  }
  if(d.i > 100){
  varTypes   <- c("N","N","N","N","N","N")        ##Variable Types
   nN        <- 6
   nC        <- 0
   tryMods   <- tryMods1 <- .mclustModel							   ##Range of Models to consider
  }
  which.pri <- "MOD"                             ##Prior: NULL,"MOD","DEF"
  mc        <- TRUE
  g1        <- 2;
  g2        <- 4;                                ##Range of G values to consider
  nSim1     <- as.numeric(Args[6])
  nSim2     <- as.numeric(Args[7])
  resCol    <- c("Meth","Sim",rep(paste(1:(nN+nC),varTypes,sep=":"),2),"G","Mod","bic","Rand","Time (min)")
}

Flag <-0
if(typeof(d)!="list"){
 Flag <- 1;   cat("Data not a data.frame \n")
}
if(is.null(resN)){
 Flag <- 1;   cat("Results file not specified \n")
}
if(is.null(G)){
 Flag <- 1;   cat("G not specified \n")
}
if(is.null(.startVar)){
 Flag <- 1;   cat(".startVar not specified \n")
}
if(is.null(varTypes)){
 Flag <- 1;   cat("varTypes not specified \n")
}
if(is.null(g1) | is.null(g2)){
 Flag <- 1;   cat("g1 and/or g2 not specified \n")
}
if(is.null(tryMods)){
 Flag <- 1;   cat("tryMods not specified \n")
}
if(Flag==0){
##----------------------------------------------------------
##Prelim
if(runPre){

 pre <- numeric()
 for(nS in nSim1:nSim2){
  use.d <- d[d[,8]==nS,.startVar]
  u     <- d[d[,8]==nS,7]
 # mD <- makeData(NULL,1:(nN+nC),"Excl",varTypes,use.d)
  a1=Sys.time()
  all <- whichClust(NULL,.startVar,"Excl",use.d,g1,g2,varTypes,tryMods,which.pri,kp,tol,max.iter,mc); b1=Sys.time(); a2=Sys.time()
  nor <- whichClust(NULL,1:nN,"Excl",use.d,g1,g2,varTypes,tryMods,which.pri,kp,tol,max.iter,mc); b2=Sys.time(); a3=Sys.time()
  tru <- whichClust(NULL,which(useS=="S"),"Excl",use.d,g1,g2,varTypes,tryMods,which.pri,kp,tol,max.iter,mc); b3=Sys.time(); a4=Sys.time()
  truG <- whichClust(NULL,which(useS=="S"),"Excl",use.d,G,G,varTypes,tryMods,which.pri,kp,tol,max.iter,mc); b4=Sys.time()
  pre <- rbind(pre,
   c("All",nS,rep("S",nN+nC),rep(NA,nN+nC),all$G,all$model,round(all$bic,2),round(adjustedRandIndex(u,all$class),2),round(difftime(b1,a1,units='mins'),2)),
   c("Nor",nS,rep("S",nN),NA,NA,rep(NA,nN+nC),nor$G,nor$model,round(nor$bic,2),round(adjustedRandIndex(u,nor$class),2),round(difftime(b2,a2,units='mins'),2)),
   c("Tru",nS,useS,rep(NA,nN+nC),tru$G,tru$model,round(tru$bic,2),round(adjustedRandIndex(u,tru$class),2),round(difftime(b3,a3,units='mins'),2)),
   c("TruG",nS,useS,rep(NA,nN+nC),truG$G,truG$model,round(truG$bic,2),round(adjustedRandIndex(u,truG$class),2),round(difftime(b4,a4,units='mins'),2)))
 }
write.table(pre,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
sendEmail(subject="Prelim",text=paste(d.i,": Done!",sep=""),address="keighto@gmail.com")
}
##---------------------------------------------------------------------------------
##-----------------------------------------DV--------------------------------------
##D(Forward), M(Backward), V=Model, G varied at each step
run.DV = run.MV = list()
if(runDV){
  ##Dean
  for(nS in nSim1:nSim2){
    res.DV = matrix(NA,1,length(resCol)) 
    .DV1 <- .DVran <- numeric(); a = Sys.time()
    use.d <- d[d[,8]==nS,.startVar]
    u     <- d[d[,8]==nS,7]
    run.DV[[nS]] <-  stepSelect(.startVar,g1,g2,tryMods,use.d,varTypes,which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=TRUE,"Dean"); b=Sys.time()
    .DV1 <-  c(run.DV[[nS]]$SW$fBIC$clust$G,run.DV[[nS]]$SW$fBIC$clust$model,round(run.DV[[nS]]$SW$fBIC$clust$bic,2))
    .DVran <- round(adjustedRandIndex(u,run.DV[[nS]]$SW$fBIC$clust$class),2)
    res.DV[1,] <- c(paste("DV.",R,sep=""),nS,run.DV[[nS]]$SW$allRoles,run.DV[[nS]]$SRUW$allRoles,.DV1,.DVran,round(difftime(b,a,units='mins'),2))
foo <- rbind(foo,res.DV)
#    if(saveRes) save(run.DV,file=paste("DV.",R,d.i,"_",N,".rdata",sep=""))
#    write.table(res.DV,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
  }
  sendEmail(subject="DV",text=paste(d.i,":",nS,":",R),address="keighto@gmail.com")
}
##----------------------------------------MV----------------------------------------
if(runMV){
  ##Maug
  for(nS in nSim1:nSim2){
    res.MV = matrix(NA,1,length(resCol)) 
    .MV1 <- .MVran <- numeric(); a = Sys.time()
    use.d = d[d[,8]==nS,.startVar]
    u     = d[d[,8]==nS,7]
    run.MV[[nS]] <-  stepSelect(.startVar,g1,g2,tryMods,use.d,varTypes,
                                which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=TRUE,"Maug"); b=Sys.time()
    .MV1 <-  c(run.MV[[nS]]$SW$fBIC$clust$G,run.MV[[nS]]$SW$fBIC$clust$model,round(run.MV[[nS]]$SW$fBIC$clust$bic,2))
    .MVran <- round(adjustedRandIndex(table(u,run.MV[[nS]]$SW$fBIC$clust$class)),2)
    res.MV[1,] <- c(paste("MV.",R,sep=""),nS,run.MV[[nS]]$SW$allRoles,run.MV[[nS]]$SRUW$allRoles,.MV1,.MVran,round(difftime(b,a,units='mins'),2))
    if(saveRes) save(run.MV,file=paste("MV.",R,d.i,"_",N,".rdata",sep=""))
    write.table(res.MV,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
  }
  sendEmail(subject="MV",text=paste(d.i,":",nS,":",R),address="keighto@gmail.com")
}
##---------------------------------------------------------------------------------
##------------------------------------------DF-------------------------------------
tryMods   <- paste(rep(unique(substr(tryMods1,1,3)),each=2),c(".I",".R"),sep="") ##Range of Models to consider
##Dean
run.DF.S <- run.DF.SW <- run.DF.SRUW <- list()
if(runDF){
  for(nS in nSim1:nSim2){
  res.DF.S <- res.DF.SW <- res.DF.SRUW <- matrix(NA,1,length(resCol))
  tmp.DF <- list()
  modBic <- array(NA,dim=c(g2-g1+1,length(tryMods),3))
  dimnames(modBic) <- list(paste(g1:g2),tryMods,c("S","SW","SRUW")); 
  use.d <- d[d[,8]==nS,.startVar]
  u     <- d[d[,8]==nS,7]
  aT=Sys.time()
   for(fg in g1:g2){
    tmp.DF[[fg]] <- list()
     for(fm in tryMods){
      tmp.DF[[fg]][[fm]] <- list()
      tmp.DF[[fg]][[fm]] <-  stepSelect(.startVar,fg,fg,paste(fm,c("I","J"),sep=""),use.d,varTypes,which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=TRUE,"Dean")
      if(length(tmp.DF[[fg]][[fm]])>0){
       modBic[paste(fg),fm,"S"]    <- tmp.DF[[fg]][[fm]]$SW$fBIC$clust$bic
       modBic[paste(fg),fm,"SW"]   <- max(tmp.DF[[fg]][[fm]]$SW$fBIC$resBIC,na.rm=TRUE)
       modBic[paste(fg),fm,"SRUW"] <- max(tmp.DF[[fg]][[fm]]$SRUW$fBIC$resBIC,na.rm=TRUE)
  }}}; bT=Sys.time()
  .bst <- matrix(NA,3,2)
  .bst[1,] <- which(modBic[,,"S"]==max(modBic[,,"S"],na.rm=TRUE),arr.ind=TRUE)[1,]
  .bst[2,] <- which(modBic[,,"SW"]==max(modBic[,,"SW"],na.rm=TRUE),arr.ind=TRUE)[1,]
  .bst[3,] <- which(modBic[,,"SRUW"]==max(modBic[,,"SRUW"],na.rm=TRUE),arr.ind=TRUE)[1,]
  tmp <- NULL
  for(b in 1:3){
   for(i in 1:2) tmp <- c(tmp,dimnames(modBic)[[i]][.bst[b,i]])
  }
  run.DF.S[[nS]]    <- tmp.DF[[as.numeric(tmp[1])]][[tmp[2]]]
  run.DF.SW[[nS]]   <- tmp.DF[[as.numeric(tmp[3])]][[tmp[4]]]
  run.DF.SRUW[[nS]] <- tmp.DF[[as.numeric(tmp[5])]][[tmp[6]]]
  .DF1.S      <-  c(run.DF.S[[nS]]$SW$fBIC$clust$G,run.DF.S[[nS]]$SW$fBIC$clust$model,round(run.DF.S[[nS]]$SW$fBIC$clust$bic,2))
  .DF1.SW     <-  c(run.DF.SW[[nS]]$SW$fBIC$clust$G,run.DF.SW[[nS]]$SW$fBIC$clust$model,round(run.DF.SW[[nS]]$SW$fBIC$clust$bic,2))
  .DF1.SRUW   <-  c(run.DF.SRUW[[nS]]$SW$fBIC$clust$G,run.DF.SRUW[[nS]]$SW$fBIC$clust$model,round(run.DF.SRUW[[nS]]$SW$fBIC$clust$bic,2))
  .DFran.S    <- round(adjustedRandIndex(u,run.DF.S[[nS]]$SW$fBIC$clust$class),2)
  .DFran.SW   <- round(adjustedRandIndex(u,run.DF.SW[[nS]]$SW$fBIC$clust$class),2)
  .DFran.SRUW <- round(adjustedRandIndex(u,run.DF.SRUW[[nS]]$SW$fBIC$clust$class),2)
  res.DF.S[1,]    <- c(paste("DF.S.",R,sep=""),nS,run.DF.S[[nS]]$SW$allRoles,run.DF.S[[nS]]$SRUW$allRoles,.DF1.S,.DFran.S,round(difftime(bT,aT,units='mins'),2))
  res.DF.SW[1,]   <- c(paste("DF.SW.",R,sep=""),nS,run.DF.SW[[nS]]$SW$allRoles,run.DF.SW[[nS]]$SRUW$allRoles,.DF1.SW,.DFran.SW,round(difftime(bT,aT,units='mins'),2))
  res.DF.SRUW[1,] <- c(paste("DF.SRUW.",R,sep=""),nS,run.DF.SRUW[[nS]]$SW$allRoles,run.DF.SRUW[[nS]]$SRUW$allRoles,.DF1.SRUW,.DFran.SRUW,round(difftime(bT,aT,units='mins'),2))
  run.DF <- list(DF.S=run.DF.S,DF.SW=run.DF.SW,DF.SRUW=run.DF.SRUW)
 if(saveRes) save(run.DF,file=paste("DF.",R,d.i,"_",N,".rdata",sep=""))
write.table(res.DF.S,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
write.table(res.DF.SW,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
write.table(res.DF.SRUW,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
}
  sendEmail(subject="DF",text=paste(d.i,":",nS,":",R),address="keighto@gmail.com")
}
##------------------------------------------------------------------------------------------------
##--------------------------------------------MF--------------------------------------------------
tryMods   <- paste(rep(unique(substr(tryMods1,1,3)),each=2),c(".I",".R"),sep="") ##Range of Models to consider
if(runMF){
 # paste(tryMods,d.i,whichRun,nSim1,nSim2,resN,sep=",")
 run.MF.S <- run.MF.SW <- run.MF.SRUW <- list()
 for(nS in nSim1:nSim2){
  res.MF.S <- res.MF.SW <- res.MF.SRUW <- matrix(NA,1,length(resCol))
  tmp.MF <- list()
  modBic <- array(NA,dim=c(g2-g1+1,length(tryMods),3))
  dimnames(modBic) <- list(paste(g1:g2),tryMods,c("S","SW","SRUW"));
  use.d <- d[d[,8]==nS,.startVar]
  u     <- d[d[,8]==nS,7]
  aT=Sys.time()
  for(fg in g1:g2){
   tmp.MF[[fg]] <- list()
   for(fm in tryMods){
    tmp.MF[[fg]][[fm]] <- list()
    tmp.MF[[fg]][[fm]] <-  stepSelect(.startVar,fg,fg,paste(fm,c("I","J"),sep=""),use.d,varTypes,
		which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=TRUE,"Maug")
    if(length(tmp.MF[[fg]][[fm]])>0){
     modBic[paste(fg),fm,"S"]    <- tmp.MF[[fg]][[fm]]$SW$fBIC$clust$bic
     modBic[paste(fg),fm,"SW"]   <- max(tmp.MF[[fg]][[fm]]$SW$fBIC$resBIC,na.rm=TRUE)
     modBic[paste(fg),fm,"SRUW"] <- max(tmp.MF[[fg]][[fm]]$SRUW$fBIC$resBIC,na.rm=TRUE)
   }}}; bT=Sys.time()
  .bst <- matrix(NA,3,2)
  .bst[1,] <- which(modBic[,,"S"]==max(modBic[,,"S"],na.rm=TRUE),arr.ind=TRUE)[1,]
  .bst[2,] <- which(modBic[,,"SW"]==max(modBic[,,"SW"],na.rm=TRUE),arr.ind=TRUE)[1,]
  .bst[3,] <- which(modBic[,,"SRUW"]==max(modBic[,,"SRUW"],na.rm=TRUE),arr.ind=TRUE)[1,]
  tmp <- NULL
  for(b in 1:3){
   for(i in 1:2) tmp <- c(tmp,dimnames(modBic)[[i]][.bst[b,i]])
  }
  run.MF.S[[nS]]    <- tmp.MF[[as.numeric(tmp[1])]][[tmp[2]]]
  run.MF.SW[[nS]]   <- tmp.MF[[as.numeric(tmp[3])]][[tmp[4]]]
  run.MF.SRUW[[nS]] <- tmp.MF[[as.numeric(tmp[5])]][[tmp[6]]]
  .MF1.S      <-  c(run.MF.S[[nS]]$SW$fBIC$clust$G,run.MF.S[[nS]]$SW$fBIC$clust$model,round(run.MF.S[[nS]]$SW$fBIC$clust$bic,2))
  .MF1.SW     <-  c(run.MF.SW[[nS]]$SW$fBIC$clust$G,run.MF.SW[[nS]]$SW$fBIC$clust$model,round(run.MF.SW[[nS]]$SW$fBIC$clust$bic,2))
  .MF1.SRUW   <-  c(run.MF.SRUW[[nS]]$SW$fBIC$clust$G,run.MF.SRUW[[nS]]$SW$fBIC$clust$model,round(run.MF.SRUW[[nS]]$SW$fBIC$clust$bic,2))
  .MFran.S    <- round(adjustedRandIndex(u,run.MF.S[[nS]]$SW$fBIC$clust$class),2)
  .MFran.SW   <- round(adjustedRandIndex(u,run.MF.SW[[nS]]$SW$fBIC$clust$class),2)
  .MFran.SRUW <- round(adjustedRandIndex(u,run.MF.SRUW[[nS]]$SW$fBIC$clust$class),2)
  res.MF.S[1,]    <- c(paste("MF.S.",R,sep=""),nS,run.MF.S[[nS]]$SW$allRoles,run.MF.S[[nS]]$SRUW$allRoles,.MF1.S,.MFran.S,round(difftime(bT,aT,units='mins'),2))
  res.MF.SW[1,]   <- c(paste("MF.SW.",R,sep=""),nS,run.MF.SW[[nS]]$SW$allRoles,run.MF.SW[[nS]]$SRUW$allRoles,.MF1.SW,.MFran.SW,round(difftime(bT,aT,units='mins'),2))
  res.MF.SRUW[1,] <- c(paste("MF.SRUW.",R,sep=""),nS,run.MF.SRUW[[nS]]$SW$allRoles,run.MF.SRUW[[nS]]$SRUW$allRoles,.MF1.SRUW,.MFran.SRUW,round(difftime(bT,aT,units='mins'),2))
  run.MF <- list(MF.S=run.MF.S,MF.SW=run.MF.SW,MF.SRUW=run.MF.SRUW)
# if(saveRes) save(run.MF,file=paste("MF.",R,d.i,"_",N,".rdata",sep=""))
write.table(res.MF.S,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
write.table(res.MF.SW,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
write.table(res.MF.SRUW,file=paste(resN,".txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
}
#sendEmail(subject="MF",text=paste(d.i,":",nS,":",R),address="keighto@gmail.com")
}
}
#q()
