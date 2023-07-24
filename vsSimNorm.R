Args <- commandArgs(TRUE)

if(getwd()!="/u/kevans/Research Files/newVS"){
if(getwd()=="C:/Documents and Settings/kevans/My Documents"){
 setwd("C://Documents and Settings//kevans//My Documents//My Dropbox//tmp//newVS//")
}else{
 setwd("C://Users//K80//Documents//My Dropbox//tmp//newVS//")
}}

source("vsLoad.R")
source("vsFunctionRegStep.R")
whichReg <- stepReg
source("vsFunction.R")

##REVISIT 24, 41, 44, 51, 52 ~53, 54, 61,62,63
##Plan: 1. Look at DV.S results (from school) -- How "close" was the decision -
##For atleast 1, 2 -- Sample 25

sav.plot  = FALSE
nSim1=21; nSim2=40; nSim=20
#-------------------------------------------------
tryMods1 = .mclustModel
.startVar=1:6; g1=2;g2=4;varTypes=rep("N",6); which.pri="MOD"
#-------------------------------------------------
#names1 = cbind("Ind","Name","G","N","Roles")
#names2 = cbind("Ind","Name","N1","N2","N3","N4","N5","N6")
#write.table(names1,file="loadData.txt",row.names=FALSE,col.names=FALSE,sep=",")
#write.table(names2,file="loadRole.txt",row.names=FALSE,col.names=FALSE,sep=",")
##-------------------------------------------------
##-----------------101 EII------------------
##mixClust works ok with N=250
##VS working
 ##Set + Write out Truth to external files --> Used in vsLoad, vsRun, vsSimRes
for(d.i in 101:105){
 nSim=20
 if(d.i==101) d.ind="EII"; 
  if(d.i==102) d.ind="EEI"; 
  if(d.i==103) d.ind="EEE"; 
  if(d.i==104) d.ind="VII";
  if(d.i==105) d.ind="VVV"; 

  G=3; N=250; 
  SRUW= c("S=(N1,N2),U=(N3,N4,N5,N6)"); useS = cbind(d.i,d.ind,"S","S","U","U","W","W") 
#  write.table(cbind(d.i,d.ind,G,N,SRUW),file="loadData.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
#  write.table(useS,file="loadRole.txt",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
  ##---------
  sN = which(varTypes=="N" & useS[3:8]=="S"); l.sN = length(sN)
  nN = 6;
  mu = cbind(c(3,3),c(0,0),c(-2,-2))
  var= c(1,1)
  Gvar = array(NA,dim=c(l.sN,l.sN,G))
  if(d.i==101){ 
    mu = cbind(c(3,3),c(1.5,1.5),c(-1.5,-1.5))
    var= c(1,1)
    Gvar[1:l.sN,1:l.sN,1:G] = diag(l.sN)
  }
  if(d.i==102){ 
    mu = cbind(c(3,3),c(1.5,1.5),c(-1.5,-1.5))
    var= c(1,1)
    Gvar[,,1:G] <- c(1.5,.5)*diag(l.sN)
  }
  if(d.i==103){
     mu = cbind(c(3,3),c(.5,.5),c(-1.5,-1.5))
     var= c(1,1)
     Gvar[,,1:G] <- rbind(c(1.5,-.5),c(-.5,.75))  
  }
  if(d.i==104){ 
    mu = cbind(c(3,3),c(.5,.5),c(-2,-2))
    var= c(1,1)
    Gvar[,,1:2] = .5*diag(l.sN); Gvar[,,3] = 2*diag(l.sN) 
  }
  if(d.i==105){ 
    mu = cbind(c(-1,-1),c(1,1),c(3,-3))
    var= c(1,1)
    Gvar[,,1] = rbind(c(.5,0,0,.5)); 
    Gvar[,,2] = rbind(c(2,-.5,-.5,.5)); 
    Gvar[,,3] = rbind(c(.5,.75,.75,3)) 
  }

  ##--------
  d <- matrix(NA,N*nSim,nN+2) 
   for(nS in nSim1:nSim2){
    u = sample(1:G,N,replace=TRUE)
    useNorm <- matrix(NA,N,nN); 
    for(i in 1:N) useNorm[i,1:2] <- rmvnorm(1,mu[,u[i]],var*Gvar[,,u[i]])
    for(i in 1:N) useNorm[i,3:4] <- rmvnorm(1,useNorm[i,1:2]*sample(seq(.7,1.3,by=.05),2),
                                              matrix(c(1,.9,.9,1),2,2))
    useNorm[,5:6] <- rmvnorm(N,c(0,0),c(2)*diag(2))
    d[(1:N)+(nS-(nSim1-1)-1)*N,] <- cbind(useNorm,u,rep(nS,N))
  }
#write.table(d,file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}
##--------
#  if(sav.plot) pdf(file=paste("DatPlot",d.i,".pdf",sep=""))
#  for(nS in 1:nSim){
#    if(!sav.plot) dev.new();
    my.pairs(d[d[,8]==nS,1:(nN+nC)],d[d[,8]==nS,7],NA);   mtext(side=3,outer=TRUE,cex=1.25,paste(d.ind,": ",SRUW,": N=",N,sep=""))  
#    }; dev.off(); 
if(chk){
a = vs = list(); d = data.frame(d);

truth <- FALSE; nS <- nSim1-1; i <- 0
#need <- 2
while(i!=need & nS<nSim2 ){
  nS <- nS+1
  tmp <- whichClust(NULL,c(1,2),"Excl",subset(d[,.startVar],d[,8]==nS),g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)
  truth <- tmp$G==G & tmp$model==d.ind
  i <- i+ifelse(truth,1,0)
}
if(i==need & nS <=nSim2){
write.table(subset(d,d[,8]<=nS ),file=paste("d",d.i,".csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
}


a = vs = list(); d = data.frame(d);
for(nS in 1:nSim) a[[nS]] = whichClust(NULL,c(1,2),"Excl",d[d[,8]==nS,.startVar],g1,g2,varTypes,tryMods1,"MOD",kp,tol,max.iter,TRUE)                                      
for(nS in 1:nSim) vs[[nS]]= stepSelect(.startVar,2,g2,tryMods1,d[d[,8]==nS,.startVar],varTypes,which.pri,.kp,.tol,.max.iter,TRUE,TRUE,"Dean")
}
