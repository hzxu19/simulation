Args <- commandArgs(TRUE)

if(getwd()!="/u/kevans/Research Files/newVS"){
if(getwd()=="C:/Documents and Settings/kevans/My Documents"){
 setwd("C://Documents and Settings//kevans//My Documents//My Dropbox//tmp//newVS//")
}else{
 setwd("C://Users//K80//Documents//My Dropbox//tmp//newVS//")
}}

source("vsLoad.R")
source("vsFunctionStepReg.R")
whichReg <- stepReg
source("vsFunction.R")
##REVISIT 24, 41, 44, 51, 52 ~53, 54, 61,62,63
##Plan: 1. Look at DV.S results (from school) -- How "close" was the decision -
##For atleast 1, 2 -- Sample 25

sav.plot  = FALSE
#-------------------------------------------------
tryMods1 = paste(rep(.mclustModel,each=4),c(".II",".RI",".IJ",".RJ"),sep="") 
.startVar=1:6; g1=2;g2=4;varTypes=c(rep("N",4),"C","C"); which.pri="MOD"
nSim=20
nSim1=21; nSim2=40
#-------------------------------------------------
names1 = cbind("Ind","Name","G","N","Roles")
names2 = cbind("Ind","Name","N1","N2","N3","N4","C1","C2")
write.table(names1,file="loadData.txt",row.names=FALSE,col.names=FALSE,sep=",")
write.table(names2,file="loadRole.txt",row.names=FALSE,col.names=FALSE,sep=",")
##-------------------------------------------------
##-----------------1: EII.II (1 Cat)------------------
##mixClust works ok with N=250
##VS working
 ##Set + Write out Truth to external files --> Used in vsLoad, vsRun, vsSimRes
 for(d.i in c(11:14)){
  if(d.i==11) d.ind="EII.II"; 
  if(d.i==12) d.ind="EEI.II"; 
  if(d.i==13) d.ind="EEE.II"; 
  if(d.i==14) d.ind="VII.II";
  if(d.i==15) d.ind="EII.II"; 

  G=3; N=250; 
  if(d.i < 15){  SRUW= c("S=(N1,N2,C1),U=(N3,N4,C2)"); useS = cbind(d.i,d.ind,"S","S","U","U","S","U") 
  }else{ SRUW= c("S=(N1,N2,C1),W=(N3,N4,C2)"); useS = cbind(d.i,d.ind,"S","S","W","W","S","W")  }    

  write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  ##---------
  sN = which(varTypes=="N" & useS[3:8]=="S"); l.sN = length(sN)
  nN = 4; nC = 2
  prob=list()
  prob[[1]] = rbind(c(.9,.05,.05),c(.05,.9,.05),c(.05,.05,.9))
  prob[[2]] = rbind(c(.95,.05),c(.05,.95),c(.95,.05),c(.05,.95))
  mu = cbind(c(3,3),c(1.5,1.5),c(-1.5,-1.5))
  .var= c(1,1)
  G.var = array(NA,dim=c(l.sN,l.sN,G))
  if(d.i==11 | d.i==15){ 
    if(d.i==15) prob[[2]] = rbind(rep(1/4,4))  
    mu = cbind(c(3,3),c(1.5,1.5),c(-1.5,-1.5))
    .var= c(1,1)
    G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)
  }
  if(d.i==12){ 
    mu = cbind(c(3,3),c(1.5,1.5),c(-1.5,-1.5))
    .var= c(1,1)
    G.var[,,1:G] <- c(1.5,.5)*diag(l.sN)
  }
  if(d.i==13){
     mu = cbind(c(2.5,2.5),c(1.5,1.5),c(-1.5,-1.5))
     .var= c(1,1)
     G.var[,,1:G] <- rbind(c(1.5,-.5),c(-.5,.75))  
  }
  if(d.i==14){ 
    mu = cbind(c(3,3),c(1.5,1.5),c(-1.5,-1.5))
    .var= c(1,1)
    G.var[,,1:2] = .5*diag(l.sN); G.var[,,3] = 2*diag(l.sN) 
  }
  ##--------
  d <- matrix(NA,N*nSim,nN+nC+2) 
   for(nS in nSim1:nSim2){
    u = sample(1:G,N,replace=TRUE)
    useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL
    chk2 <- NULL
    for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u[i]],.var*G.var[,,u[i]])
    for(i in 1:N) useCat1[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][u[i],],replace=TRUE)
    if(d.i<15){ ##Create U .variables
      for(i in 1:N) useNorm[i,3:4] <- rmvnorm(1,useNorm[i,1:2]*sample(seq(.7,1.3,by=.05),2),
                                              matrix(c(1,.9,.9,1),2,2))
      for(i in 1:N) chk2[i]    <- ifelse(useNorm[i,1]<=0,1,2)
      for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[2]]),1,prob=prob[[2]][chk2[i],],replace=TRUE)     
    }else{      ##Create W .variables
      useNorm[1:N,3:4] <- rmvnorm(N,c(2,2),diag(2)*2)
      for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[2]]),1,prob=prob[[2]],replace=TRUE)    
    }
    d[(1:N)+(nS-(nSim1-1)-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
  }
#write.table(d,file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}
##--------
  if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
  for(nS in 1:nSim){
    if(!sav.plot) dev.new();
    my.pairs(d[d[,8]==nS,1:(nN+nC)],d[d[,8]==nS,7],NA);   mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,": ",SRUW,": N=",N,sep=""))  
    }; dev.off(); 

if(chk){
a = vs = list(); dd = data.frame(d); dd[,5] <- factor(d[,5]); dd[,6] <- factor(d[,6])
#a0=whichClust(NULL,c(1,2),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
truth <- FALSE; nS <- nSim1-1; i <- 0
need <- 2
while(!truth & i<need ){
  nS <- nS+1
  tmp <- whichClust(NULL,c(1,2,5),"Excl",subset(dd[,.startVar],d[,8]==nS),g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
  truth <- tmp$G==G & tmp$model==d.ind
  i <- i+ifelse(truth,1,0)
}
write.table(subset(d,d[,8]<=nS),file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
for(nS in nSim1:nSim2) vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}
##-------------------------------------------------
##-----------------2: EII.II (2 Cat)------------------
##mixClust works ok with N=250
##VS working
  ##Set + Write out Truth to external files --> Used in vsLoad, vsRun, vsSimRes
for(d.i in 21:24){
  if(d.i==21) d.ind="EII.II"; 
  if(d.i==22) d.ind="EEI.II"; 
  if(d.i==23) d.ind="EEE.II"; 
  if(d.i==24) d.ind="VII.II";
  if(d.i==25) d.ind="EII.II"; 
  G=4; N=250; 
  if(d.i<25){     SRUW= c("S=(N1,N2,C1,C2),U=(N3,N4)"); useS = cbind(d.i,d.ind,"S","S","U","U","S","S")   
  }else{          SRUW= c("S=(N1,N2,C1,C2),W=(N3,N4)"); useS = cbind(d.i,d.ind,"S","S","W","W","S","S") }
  #write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  #write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  ##---------
  sN = which(varTypes=="N" & useS[3:8]=="S"); l.sN = length(sN)
  nN = 4; nC = 2
  prob=list()
  prob[[1]] = rbind(c(.9,.05,.05),c(.05,.9,.05),c(.05,.05,.9))
  prob[[2]] = rbind(c(.9,.05,.05),c(.05,.9,.05),c(.05,.9,.05))
  .var= c(1,1)
  G.var = array(NA,dim=c(l.sN,l.sN,G))
  if(d.i==21 | d.i==25){ 
    mu = cbind(c(-3,-4),c(4,2),c(0,-4),c(4,0))
    .var= c(1,1)
    G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)
  }
  if(d.i==22){ 
    mu = cbind(c(-2,-4),c(4,2),c(0,-4),c(4,0))
    .var= c(1,1)
    G.var[,,1:G] <- c(1.5,.75)*diag(l.sN)
  }
  if(d.i==23){
    mu = cbind(c(-2,-4),c(4,2),c(0,-4),c(4,0))
    .var= c(1,1)
   G.var[,,1:G] <- rbind(c(.5,-.75),c(-.75,2))  
  }
  if(d.i==24){ 
    mu = cbind(c(-1.5,-4),c(4,3),c(0,-4),c(4,0))
    .var= c(1,1)
    G.var[,,c(1,3)] = .5*diag(l.sN); G.var[,,c(2,4)] = 1.5*diag(l.sN) 
  }
##--------
  d <- matrix(NA,N*nSim,nN+nC+2) 
   for(nS in nSim1:nSim2){
     useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL; chk1 <- chk2 <- NULL
    u = sample(1:G,N,replace=TRUE)
    for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u[i]],.var*G.var[,,u[i]])
    for(i in 1:N) chk1[i] <- ifelse(u[i]==1 | u[i]==3,1,ifelse(u[i]==2,2,3))
    for(i in 1:N) chk2[i] <- ifelse(u[i]==2 | u[i]==4,1,ifelse(u[i]==3,2,3))
    for(i in 1:N) useCat1[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][chk1[i],],replace=TRUE)
    for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][chk2[i],],replace=TRUE)
    if(d.i<25){ ##U .variables
      for(i in 1:N) useNorm[i,3:4] <- rmvnorm(1,useNorm[i,1:2]*sample(seq(.7,1.3,by=.05),2),
                                              matrix(c(1,.9,.9,1),2,2))
    }else{ useNorm[1:N,3:4] <- rmvnorm(N,c(2,2),diag(2)*2) } ##W .variables  
    d[(1:N)+(nS-(nSim1-1)-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
  }
#write.table(d,file=paste("d",d.i,".csv",sep=""),col.names=FALSE,row.names=FALSE,sep=",")
}
##--------
if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
 for(nS in 1:nSim){
  if(!sav.plot) dev.new();
   my.pairs(d[d[,8]==nS,1:(nN+nC)],d[d[,8]==nS,7],NA);   mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,": ",SRUW,": N=",N,sep=""))
}; dev.off(); 

if(chk){
a = vs = list(); dd = data.frame(d); dd[,5] <- factor(d[,5]); dd[,6] <- factor(d[,6])
#a0=whichClust(NULL,c(1,2),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
truth <- FALSE; nS <- nSim1-1; i <- 0
#need <- 2
while(i!=need & nS<nSim2 ){
  nS <- nS+1
  tmp <- whichClust(NULL,c(1,2,5,6),"Excl",subset(dd[,.startVar],d[,8]==nS),g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
  truth <- tmp$G==G & tmp$model==d.ind
  i <- i+ifelse(truth,1,0)
}
if(i==need & nS <=nSim2){
write.table(subset(d,d[,8]<=nS ),file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
}

a = vs = list(); d = data.frame(d); d[,5] <- factor(d[,5]); d[,6] <- factor(d[,6])
a0=whichClust(NULL,c(1,2),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
for(nS in 1:nSim) a[[nS]] = whichClust(NULL,c(1,2,5,6),"Excl",d[d[,8]==nS,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
for(nS in 1:nSim) vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}
##-----------------   3  --------------------------
##-----------------EII.IJ (2 Cat)------------------
##mixClust, VS works N=250
  ##Set + Write out Truth to external files --> Used in vsLoad, vsRun, vsSimRes
  G=3; N=250; 
for(d.i in 31:34){
  if(d.i==31) d.ind="EII.IJ"; 
  if(d.i==32) d.ind="EEI.IJ"; 
  if(d.i==33) d.ind="EEE.IJ"; 
  if(d.i==34) d.ind="VII.IJ";
  #d.i=35; d.ind="EII.IJ"; 

  if(d.i<35){     SRUW= c("S=(N1,N2,C*=(C1,C2)),U=(N3,N4)"); useS = cbind(d.i,d.ind,"S","S","U","U","S","S")  
  }else{          SRUW= c("S=(N1,N2,C*=(C1,C2)),W=(N3,N4)"); useS = cbind(d.i,d.ind,"S","S","W","W","S","S") }
  # write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  # write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  ##---------
  sN = which(varTypes=="N" & useS[3:8]=="S"); l.sN = length(sN)
  nN = 4; nC = 2
  prob=list()
  prob[[1]] <- rbind(rep(1/2,2),c(.95,.05),c(.05,.95))
  prob[[2]] <- rbind(c(.95,.05),c(.05,.95))
  .var= c(1,1)
  G.var = array(NA,dim=c(l.sN,l.sN,G))
  if(d.i==31 | d.i==35){ 
    mu = cbind(c(0,0),c(-1.5,-1.5),c(1.75,0))
    .var= c(1,1)*.75
    G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)
  }
  if(d.i==32){ 
  #  mu = cbind(c(0,0),c(-1.5,-1.5),c(1.75,0))
  #  .var= c(.5,2)
  #  G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)  
      mu = cbind(c(0,0),c(-1.5,-1.5),c(1.75,-.25))
      .var= c(.5,2)
      G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)    
    
  }
  if(d.i==33){
    mu = cbind(c(0,0),c(-2,-2),c(1.75,0))
    .var= c(1,1)
    G.var[1:l.sN,1:l.sN,1:G] = matrix(c(.5,-.6,-.6,2),2,2) 

   }
  if(d.i==34){ 
    mu = cbind(c(1,1),c(-1.5,-1.5),c(3,1))
    .var= c(1,1)*.75
    G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)
    G.var[,,1] <- 2*diag(l.sN)
   }
##--------
  d <- matrix(NA,N*nSim,nN+nC+2) 
  for(nS in nSim1:nSim2){
    useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL; chk2 <- NULL
    u = sample(1:G,N,replace=TRUE)
    for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u[i]],.var*G.var[,,u[i]])
    for(i in 1:N) useCat1[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][u[i],],replace=TRUE)
    for(i in 1:N) chk2[i] <- ifelse(u[i]==1,ifelse(useCat1[i]==1,1,2),ifelse(useCat1[i]==1,2,1))
    for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[2]]),1,prob=prob[[2]][chk2[i],],replace=TRUE)
    #if(d.i<35){ ##U .variables
     for(i in 1:N) useNorm[i,3:4] <- rmvnorm(1,useNorm[i,1:2]*sample(seq(.7,1.3,by=.05),2),
                                              matrix(c(1,.9,.9,1),2,2))
   # }else{ useNorm[1:N,3:4] <- rmvnorm(N,c(2,2),diag(2)*2) } ##W .variables
    d[(1:N)+(nS-(nSim1-1)-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
  }
#write.table(d,file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}
##--------
#  if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
#  for(nS in 1:nSim){
#    if(!sav.plot) dev.new();
  #   my.pairs(d[d[,8]==nS,1:(nN+nC)],d[d[,8]==nS,7],NA);    mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,": ",SRUW,": N=",N,sep=""))
  # }; dev.off(); 

if(rerun){
a = vs = list(); dd = data.frame(d); dd[,5] <- factor(d[,5]); dd[,6] <- factor(d[,6])
truth <- FALSE; nS <- nSim1-1; i <- 0
#need <- 4
while(i!=need & nS<nSim2 ){
  nS <- nS+1
  tmp <- whichClust(NULL,c(1,2,5,6),"Excl",subset(dd[,.startVar],d[,8]==nS),g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
  truth <- tmp$G==G & tmp$model==d.ind
  i <- i+ifelse(truth,1,0)
}
if(i==need & nS <=nSim2){
write.table(subset(d,d[,8]<=nS),file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
}
}
if(chk){
  used=d
  a = vs = list();  d = data.frame(d); d[,5] <- factor(d[,5]); d[,6] <- factor(d[,6])
  #a1=whichClust(NULL,c(1,2),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
  #for(nS in 1:nSim) 
  a[[nS]] = whichClust(NULL,c(1,2,5,6),"Excl",d[d[,8]==nS,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
  #for(nS in 1:nSim) 
    vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}
##-----------------   4  --------------------------
##-----------------Exx.RI (1 Cat)------------------
##mixClust , VS work with N=250
##Set + Write out Truth to external files --> Used in vsLoad, vsRun, vsSimRes
for(d.i in c(41,44)){
G=2;
   if(d.i==41){ d.ind="EXX.RI"; N=200 }
   if(d.i==44){ d.ind="VXX.RI"; N=250 }
   #if(d.i==45) d.ind="EXX.RI"; 

  if(d.i<45){
   SRUW= c("S=(N1,C1),U=(N3,N4,C2),W=(N2)"); useS = cbind(d.i,d.ind,"S","W","U","U","S","U")  
  }else{   SRUW= c("S=(N1,C1),W=(N2,N3,N4,C2)"); useS = cbind(d.i,d.ind,"S","W","W","W","S","W")  }
 
 #write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
 #write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
##---------
sN = which(varTypes=="N" & useS[3:8]=="S"); 
if(d.i==41 | d.i==44 | d.i==45){ l.sN = length(sN)+1 }else{ l.sN = length(sN) }
nN = 4; nC = 2
prob=list()
prob[[1]] = rbind(c(.55,.45),c(.5,.5))
prob[[2]] = rbind(c(.9,.05,.05),c(.05,.4,.45))
G.var = array(NA,dim=c(l.sN,l.sN,G))
if(d.i==41 | d.i==45){ 
  if(d.i==45) prob[[2]] = rep(1/G,G) 
  mu = cbind(c(-2,0),c(-1,0),c(1,0),c(2,0))
  .var= c(1,2)*.5
  G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)
}
if(d.i==44){ 
  #mu = cbind(c(-1.75,0),c(-1,0),c(1.5,0),c(3.5,0))
  mu = cbind(c(-1.75,0),c(-1,0),c(2.25,0),c(3.5,0))
  .var= c(1,1)
  G.var[,,1] = diag(l.sN)*c(.25,1);  
  G.var[,,2] = diag(l.sN)*c(1,1); 
}
##--------
##--------
d <- matrix(NA,N*(nSim2-nSim1+1),nN+nC+2) 
for(nS in nSim1:nSim2){
  useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL; chk2 <- NULL
  u = sample(1:G,N,replace=TRUE)
  for(i in 1:N) useCat1[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][u[i],],replace=TRUE)
  u2 <- 2*u-(useCat1%%2)
  for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u2[i]],.var*G.var[,,u[i]])
  if(d.i<45){ ##Create U .variables
    for(i in 1:N) useNorm[i,3:4] <- rmvnorm(1,useNorm[i,1:2]*sample(seq(.7,1.3,by=.05),2),
                                            matrix(c(1,.9,.9,1),2,2))
    for(i in 1:N) chk2[i] <- ifelse(useNorm[i,1]>11,1,2)
    for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[2]]),1,prob=prob[[2]][chk2[i],],replace=TRUE)
  }else{      ##Create W .variables
    useNorm[1:N,3:4] <- rmvnorm(N,c(2,2),diag(2)*2)
    for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[2]]),1,prob=prob[[2]],replace=TRUE)    
  }
  d[(1:N)+(nS-(nSim1-1)-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
}
#write.table(d,file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}
if(rerun){
a = vs = list(); dd = data.frame(d); dd[,5] <- factor(d[,5]); dd[,6] <- factor(d[,6])
truth <- FALSE; nS <- nSim1-1; i <- 0
#need <- 4
while(i!=need & nS<nSim2 ){
  nS <- nS+1
  tmp <- whichClust(NULL,c(1,5,6),"Excl",subset(dd[,.startVar],d[,8]==nS),g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
  truth <- tmp$G==G & tmp$model==d.ind
  i <- i+ifelse(truth,1,0)
}
if(i==need & nS <=nSim2){
write.table(subset(d,d[,8]<=nS),file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
}
}
##--------
if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
 for(nS in 1:nSim){
 if(!sav.plot) dev.new();
 used <- d[d[,8]==nS,];    u2 <- 2*used[,7]-(used[,5]%%2)
 my.pairs(used[,1:6],used[,7],NA);   mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,":",SRUW,"N=",N,"(True Class)",sep=""))
 my.pairs(used[,1:6],u2,NA);    mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,":",SRUW,"N=",N,"(Condit on C1)",sep=""))
}; dev.off(); 



if(chk){
  a = vs = list();   d = data.frame(d); d[,5] <- factor(d[,5]); d[,6] <- factor(d[,6])
  a0=whichClust(NULL,c(1),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
  for(nS in 1:5) a[[nS]] = whichClust(NULL,c(1,5),"Excl",d[d[,8]==nS,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
  for(nS in 1:2) vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}
##-----------------   5  --------------------------
##-----------------Exx.RI (1 Cat)------------------
for(d.i in 51:54){
G=3; N=500; 

  if(d.i==51) d.ind="EII.RI"; 
  if(d.i==52) d.ind="EEI.RI"; 
  if(d.i==53) d.ind="EEE.RI"; 
  if(d.i==54) d.ind="VII.RI";
  if(d.i==55) d.ind="EII.RI"; 

 if(d.i<55){
  SRUW= c("S=(N1,N2,C1,C2),U=(N3,N4)"); useS = cbind(d.i,d.ind,"S","S","U","U","S","S")  
 }else{
  SRUW= c("S=(N1,N2,C1,C2),W=(N3,N4)"); useS = cbind(d.i,d.ind,"S","S","W","W","S","S")  
 }
 #write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
 #write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
##---------
sN = which(varTypes=="N" & useS[3:8]=="S"); l.sN = length(sN)
nN = 4; nC = 2
prob=list()
prob[[1]] = rbind(c(.50,.5),c(.5,.5),c(.5,.5))
prob[[2]] = rbind(c(.5,.5),c(.5,.5),c(.5,.5))
G.var = array(NA,dim=c(l.sN,l.sN,G))
if(d.i==51 | d.i==55){ 
  mu11 = c(-1.5,-3.25)
  uc1 = c(-1.25,.75); uc2=c(.9,1.15)
  mu1 = cbind(mu11,mu11+uc1,mu11+uc2,mu11+uc1+uc2); 
  mu = cbind(mu1,mu1+c(0,-4),mu1+c(4,-1))
  .var= c(1,1)*.5
  G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)
}
if(d.i==52){ 
  mu11 = c(-1.5,-3.25)
  uc1 = c(-2.5,.75); uc2=c(1.25,1)
  mu1 = cbind(mu11,mu11+uc1,mu11+uc2,mu11+uc1+uc2); 
  mu = cbind(mu1,mu1+c(0,-4),mu1+c(8,-3))
  .var= c(4,.5)*.5
  G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN)
}
if(d.i==53){    
  mu11 = c(-1.5,-3.25)
  uc1 = c(-2.5,.75); uc2=c(1.25,1)
  mu1 = cbind(mu11,mu11+uc1,mu11+uc2,mu11+uc1+uc2);   
  mu = cbind(mu1,mu1+c(0,-5),mu1+c(7,0))
  .var= c(1,1)
  G.var[1:l.sN,1:l.sN,1:G] = rbind(c(1.5,-.75),c(-.75,.75))
}
if(d.i==54){ 
  mu11 = c(-1.5,-3.25)
  uc1 = c(-1.25,.75); uc2=c(.9,1.15)
  mu1 = cbind(mu11,mu11+uc1,mu11+uc2,mu11+uc1+uc2); 
  mu = cbind(2*mu1,mu1+c(0,-8),1.5*mu1+c(5,-1))
  .var= c(1,1)
  G.var[,,1] = 1.5*diag(l.sN); G.var[,,2] = .5*diag(l.sN); G.var[,,3] = .75*diag(l.sN)     
}
##--------
d <- matrix(NA,N*nSim,nN+nC+2) 
for(nS in nSim1:nSim2){
  u = sample(1:G,N,replace=TRUE)
  useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL; u2 <- NULL
  for(i in 1:N) useCat1[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][u[i],],replace=TRUE)
  for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[2]]),1,prob=prob[[2]][u[i],],replace=TRUE)
  for(i in 1:N) u2[i] <- 4*(u[i]-1)+1*ifelse(useCat1[i]==1,useCat2[i],useCat2[i]+2)
  for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u2[i]],.var*G.var[,,u[i]])
  if(d.i<55){ ##Create U .variables
    for(i in 1:N) useNorm[i,3:4] <- rmvnorm(1,useNorm[i,1:2]*sample(seq(.75,1.25,by=.05),2),
                                            matrix(c(1,.9,.9,1),2,2))
  }else{ useNorm[1:N,3:4] <- rmvnorm(N,c(2,2),diag(2)*2) } ##Create W .variables
  d[(1:N)+(nS-(nSim1-1)-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
}
#write.table(d,file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}
if(rerun){
a = vs = list(); dd = data.frame(d); dd[,5] <- factor(d[,5]); dd[,6] <- factor(d[,6])
truth <- FALSE; nS <- nSim1-1; i <- 0
#need <- 4
while(i!=need & nS<nSim2 ){
  nS <- nS+1
  tmp <- whichClust(NULL,c(1,2,5,6),"Excl",subset(dd[,.startVar],d[,8]==nS),g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
  truth <- tmp$G==G & tmp$model==d.ind
  i <- i+ifelse(truth,1,0)
}
if(i==need & nS <=nSim2){
write.table(subset(d,d[,8]<=nS),file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
}
}
##----------  
if(sav.plot) pdf(file=paste("DatPlot",d.ind,".pdf",sep=""))
for(nS in 1:nSim){
  if(!sav.plot) dev.new();
  used <- d[d[,8]==nS,]; u2 <- NULL
  for(i in 1:N) u2[i] <- 4*(used[i,7]-1)+1*ifelse(used[i,5]==1,used[i,6],used[i,6]+2)
  my.pairs(used[,1:6],used[,7],NA);    mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind," :R,S=(N1,N2,C1,C2),W=(N3,N4),N=",N,"(True Class)",sep=""))
  my.pairs(used[,1:3],u2,NA);   mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind," :R,S=(N1,N2,C1,C2),W=(N3,N4),N=",N,"(Condit on 5:C)",sep=""))
}; dev.off(); 

if(chk){
  d = data.frame(d); d[,5] <- factor(d[,5]); d[,6]<- factor(d[,6]);   a = vs = list()
  a0=whichClust(NULL,c(1,2),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
  for(nS in 1:5) a[[nS]] = whichClust(NULL,c(1,2,3,5,6),"Excl",d[d[,8]==nS,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
  for(nS in 1:5) vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}

##---------------------   6  --------------------
##--------------------EII.RJ---------------------
##mixClust, vs work N=500
##mixClust, vs work N=250
 ##Set + Write out Truth to external files --> Used in vsLoad, vsRun, vsSimRes
 for(d.i in c(61,64)){
  G=2; N=250;
 if(d.i==61) d.ind="EXX.RJ"; 
  #if(d.i==62) d.ind="EEI.RJ"; 
  #if(d.i==63) d.ind="EEE.RJ"; 
  if(d.i==64) d.ind="VXX.RJ";
  #if(d.i==65) d.ind="EII.RJ"; 

  if(d.i<65){
    SRUW= c("S=(N1,C*=(C1,C2)),U=(N3,N4),W=(N2)"); useS = cbind(d.i,d.ind,"S","W","U","U","S","S")
  }else{   
    SRUW= c("S=(N1,C*=(C1,C2)),W=(N2,N3,N4)"); useS = cbind(d.i,d.ind,"S","W","W","W","S","S")
  }
  # write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  #write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  ##---------
  sN = which(varTypes=="N" & useS[3:8]=="S"); l.sN = length(sN) +1
  nN = 4; nC = 2
  prob=list()
  prob[[1]] <- rbind(rep(1/2,2),rep(1/2,2))
  prob[[2]] <- rbind(c(.95,.05),c(.05,.95))
  G.var = array(NA,dim=c(l.sN,l.sN,G))
  if(d.i==61 | d.i==65){ 
   # mu = cbind(c(-2,0),c(-3,0),c(-1,0),c(3,0),c(1,0),c(2,0))
    mu = cbind(c(-2,0),c(-3,0),c(-1,0),c(2.75,0),c(1,0),c(2,0))
    .var= c(1,2)*.5
    G.var[1:l.sN,1:l.sN,1:G] = diag(l.sN) 
  }
 # if(d.i==62){   }
 # if(d.i==63){   }
  if(d.i==64){ 
    mu = cbind(c(-2.25,0),c(-3,0),c(-1.75,0),c(3,0),c(1,0),c(2,0))
    .var= c(1,1)
    G.var[,,1] = diag(l.sN)*c(.25,1);  
    G.var[,,2] = diag(l.sN)*c(1,1); 
    
  } 
##--------
d <- matrix(NA,N*(nSim2-nSim1+1),nN+nC+2) 
 for(nS in nSim1:nSim2){
    useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL; chk2 <- chk1 <- u2 <- NULL
    u = sample(1:G,N,replace=TRUE)
    for(i in 1:N) useCat1[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][u[i],],replace=TRUE)
    for(i in 1:N) useCat2[i] <- sample(1:ncol(prob[[2]]),1,prob=prob[[1]][u[i],],replace=TRUE)
    for(i in 1:N) chk2[i] <- ifelse(useCat1[i]==useCat2[i],1,ifelse(useCat1[i]<useCat2[i],2,3))
    u2 <- 3*(u-1)+chk2
    for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u2[i]],.var*G.var[,,u[i]])
    if(d.i<45){ ##Create U .variables
      for(i in 1:N) useNorm[i,3:4] <- rmvnorm(1,useNorm[i,1:2]*sample(seq(.7,1.3,by=.05),2),
                                              matrix(c(1,.9,.9,1),2,2))
    }else{      ##Create W .variables
      useNorm[1:N,3:4] <- rmvnorm(N,c(2,2),diag(2)*2)
    }
    d[(1:N)+(nS-(nSim1-1)-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
  }
#write.table(d,file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}

if(rerun){
a = vs = list(); dd = data.frame(d); dd[,5] <- factor(d[,5]); dd[,6] <- factor(d[,6])
truth <- FALSE; nS <- nSim1-1; i <- 0
#need <- 4
while(i!=need & nS<nSim2 ){
  nS <- nS+1
  tmp <- whichClust(NULL,c(1,5,6),"Excl",subset(dd[,.startVar],d[,8]==nS),g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
  truth <- tmp$G==G & tmp$model==d.ind
  i <- i+ifelse(truth,1,0)
}
if(i==need & nS <=nSim2){
write.table(subset(d,d[,8]<=nS),file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
}
}
##--------
if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
 for(nS in 1:nSim){
 if(!sav.plot) dev.new();
  used <- d[d[,8]==nS,];    
  for(i in 1:N) chk2[i] <- ifelse(used[i,5]==used[i,6],1,ifelse(used[i,5]<used[i,6],2,3))
  u2 <- 3*(used[,7]-1)+chk2
  my.pairs(used[,1:6],used[,7],NA);   mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,":",SRUW,"N=",N,"(True Class)",sep=""))
  my.pairs(used[,1:6],u2,NA);    mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,":",SRUW,"N=",N,"(Condit on C1)",sep=""))
}; dev.off();

if(chk){
  a = vs = list(); d = data.frame(d); d[,5] <- factor(d[,5]); d[,6] <- factor(d[,6])
  a0=whichClust(NULL,c(1),"Excl",d[d[,8]==1,1:6],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
  for(nS in 1:5) a[[nS]] = whichClust(NULL,c(1,5,6),"Excl",d[d[,8]==nS,1:6],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
  for(nS in 1:5) vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,1:6],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}

##------------------------------------------------------
##-----------------Exx.CI (1 Cat): UNEVEN GROUP PROPORTIONS
##mixClust , VS work with N=250
if(d.ind=="EXX.RI"){
  ##Set + Write out Truth to external files --> Used in vsLoad, vsRun, vsSimRes
  d.i= d.ind="EXX.RI"; G=2; N=250; SRUW= c("S=(N1,C1)),W=(N2,N3,N4,C2)")
  TrueS = cbind(d.i,d.ind,"S","W","W","W","S","W")  
  write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(TrueS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  ##---------
  prob=list()
  prob[[1]] = rbind(c(.55,.45),c(.5,.5))
  prob[[2]] = rep(1/G,G)#rbind(c(.9,.05,.05),c(.05,.9,.05),c(.05,.9,.05))
  mu = cbind(c(-2,0),c(-1,0),c(1,0),c(2,0))
  .var= c(1,2)*.5
  nN = 2+2; nC = 2
  d <- matrix(NA,N*nSim,nN+nC+2); 
  for(nS in 1:nSim){
    useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL
    u = sample(1:G,N,replace=TRUE,prob=c(.25,.75))
    for(i in 1:N) useCat1[i] <- sample(1:ncol(prob[[1]]),1,prob=prob[[1]][u[i],],replace=TRUE)
    u2 <- 2*u-(useCat1%%2)
    for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u2[i]],.var*diag(2))
    useNorm[1:N,3:4] <- rmvnorm(N,c(2,2),diag(2)*2)
    for(i in 1:N) useCat2[i] <- sample(1:length(prob[[2]]),1,prob=prob[[2]],replace=TRUE)
    d[(1:N)+(nS-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
  }
if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
  for(nS in 1:nSim){
    if(!sav.plot) dev.new();
    used <- d[d[,8]==nS,]
    u2 <- 2*used[,7]-(used[,5]%%2)
    my.pairs(used[,1:6],used[,7],NA);   mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,":",SRUW,"N=",N,"(True Class)",sep=""))
    my.pairs(used[,1:6],u2,NA);    mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,":",SRUW,"N=",N,"(Condit on C1)",sep=""))
  }; dev.off(); 
  write.csv(d,file=paste("d",d.i,".csv",sep=""))
}
if(chk){
  a = vs = list()
  d = data.frame(d); d[,5] <- factor(d[,5]); d[,6] <- factor(d[,6])
  a0=whichClust(NULL,c(1),"Excl",d[d[,8]==1,.startVar],1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
  for(nS in 1:5) a[[nS]] = whichClust(NULL,c(1,5),"Excl",d[d[,8]==nS,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
  for(nS in 1:5) vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}
my.pairs(used[,1:6],vs[[nS]]$SW$fBIC$clust$class,NA);    mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,":",SRUW,"N=",N,"(Condit on C1)",sep=""))


##------------------------------------------------------
##-----------------Seychelles Sim
  if(d.i==71) d.ind="EII.II"; 
  if(d.i==12) d.ind="EEI.II"; 
  if(d.i==13) d.ind="EEE.II"; 
  if(d.i==14) d.ind="VII.II";
  if(d.i==15) d.ind="EII.II"; 

  G=4; N=250; 
  #if(d.i < 15){  SRUW= c("S=(N1,N2,C1),U=(N3,N4,C2)"); useS = cbind(d.i,d.ind,"S","S","U","U","S","U") 
  #}else{ SRUW= c("S=(N1,N2,C1),W=(N3,N4,C2)"); useS = cbind(d.i,d.ind,"S","S","W","W","S","W")  }    

  #write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  #write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  ##---------
  #sN = which(varTypes=="N" & useS[3:8]=="S"); l.sN = length(sN)
  sN <- c(1,2,3); l.sN <- length(sN)
  nN = 4; nC = 2
  #prob=list()
  #prob[[1]] = rbind(c(.9,.05,.05),c(.05,.9,.05),c(.05,.05,.9))
  #prob[[2]] = rbind(c(.95,.05),c(.05,.95),c(.95,.05),c(.05,.95))
  uProb <- c(.1,.2,.35,.35)
  mu = cbind(c(1.5,2,2.25),c(1,1.5,1.5),c(1,2,3),c(1,3,3))
  .var= c(.5,.5,.375)
  G.var = array(NA,dim=c(l.sN,l.sN,G))
  G.var[,,1:G] <- diag(l.sN)
   ##--------
  d <- matrix(NA,N*nSim,nN+nC+2) 
   for(nS in 1:nSim){
    u = sample(1:G,N,replace=TRUE,prob=uProb)
    useNorm <- matrix(NA,N,nN); useCat1 <- useCat2 <- NULL
    chk2 <- NULL
    for(i in 1:N) useNorm[i,1:3] <- rmvnorm(1,mu[,u[i]],.var*G.var[,,u[i]])
    useNorm[which(u==1),1] <- 0
    useNorm[which(u==1 | u==3),2] <- 0
    useNorm[,1:2] <- abs(useNorm[,1:2])
    useNorm[which(u!=1),1] <- useNorm[which(u!=1),1]+.125
    for(i in 1:N) useCat1[i] <- ifelse(useNorm[i,1]==0,1,2)
    for(i in 1:N) useCat2[i] <- ifelse(useNorm[i,2]==0,1,2)
    useNorm[,4] <- rnorm(N,useNorm[,3]*sample(seq(.6,1.4,by=.05),N,replace=TRUE),
                          sd(useNorm[,3])*1.25)
    d[(1:N)+(nS-1)*N,] <- cbind(useNorm,factor(useCat1),factor(useCat2),u,rep(nS,N))
  }
#write.table(d,file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
#}
##--------
#  if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
#  for(nS in 1:nSim){
#    if(!sav.plot) dev.new();
my.pairs(d[d[,8]==nS,1:(nN+nC)],d[d[,8]==nS,7],NA);  # mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,": ",SRUW,": N=",N,sep=""))  
#    }; dev.off(); 
#
 for(i in 1:5) print(mean(d[d[,8]==i,2]==0))
#if(chk){
dd <- d
  a = vs = list(); d = data.frame(d); d[,5] <- factor(d[,5]); d[,6] <- factor(d[,6])
#ay=whichClust(NULL,c(1,2),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
ax=whichClust(NULL,c(1,2,3),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
#a0=whichClust(NULL,c(2,3,5,6),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
#a1=whichClust(NULL,c(1,2,3,5),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  
a2=whichClust(NULL,c(2,3,5),"Excl",d[d[,8]==1,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)  

table(ay$class,d[d[,8]==1,7])
table(ax$class,d[d[,8]==1,7])
table(a0$class,d[d[,8]==1,7])
table(a1$class,d[d[,8]==1,7])
table(a2$class,d[d[,8]==1,7])
ax$bic; a0$bic; a1$bic; a2$bic
my.pairs(dd[dd[,8]==1,1:(nN+nC)],dd[dd[,8]==1,7],NA);
my.pairs(dd[dd[,8]==1,1:(nN+nC)],a2$class,NA);


for(nS in 1:nSim) a[[nS]] = whichClust(NULL,c(1,2,5),"Excl",d[d[,8]==nS,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
for(nS in 1:nSim) vs[[nS]]<- stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
# }
