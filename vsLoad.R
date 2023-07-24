
#options(repos=structure(c(CRAN="http://cran.mtu.edu/")))
#install.packages(c("mvtnorm","e1071","epitools","flexclust"))
library(mvtnorm)
library(MASS)
library(mclust)
library(e1071)
library(epitools)
library(nnet)
#library(xtable)
library(flexclust)
library(compiler)

.kp        <<- kp        <<- .01
.tol       <<- tol       <<- 1e-5
.max.iter  <<- max.iter  <<- 50
.which.pri <<- which.pri <<- "MOD"
incl.Crit  <<- 0
excl.Crit  <<- 0
tr <- function(x) sum(diag(x))

.mclustModel <<- c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV")
.mclustUni   <<- c("EXX","VXX")
.sphrModel   <<- .mclustModel[1:2]
.diagModel   <<- .mclustModel[3:6]
.ellpModel   <<- .mclustModel[7:10]
.mixdModel   <<- c("II","RI","IJ","RJ")
.mixmodModel <<- c("Gaussian_pk_L_I","Gaussian_pk_Lk_I","Gaussian_pk_L_B","Gaussian_pk_Lk_B","Gaussian_pk_L_Bk","Gaussian_pk_Lk_Bk",
                 "Gaussian_pk_L_C","Gaussian_pk_L_Dk_A_Dk","Gaussian_pk_Lk_Dk_A_Dk","Gaussian_pk_Lk_Ck")

baseStruct <<- list()
baseStruct$N <- c(NA,"EXX","DXX","VXX")
baseStruct$C <- c(NA,"XEX","XVX")
baseStruct$M <- c(NA,"EEX","DEX","VEX","EVX","DVX","VVX")

source("vsModDiagN3.R")
source("vsModMixN3.R")

TruthVec <- read.table(file="loadData.txt",header=TRUE,sep=",")
TrueS    <- read.table(file="loadRole.txt",header=TRUE,sep=",")    

##-------------------------------------------------------------
##Universal Functions
my.pairs <- function(data,class,out){
  K = dim(data)[2]
  G = length(table(class))
  par(mfrow=c(K,K),mar=c(4,4,0,0),oma=c(0,0,2,0))
  for(k in 1:K){
    for(.k in 1:K){
      if(k ==.k){
        if(length(unique(data[,k]))>.3*nrow(data)){ ##Continuous Diag
          plot(density(data[class==1,k],bw=.5), type="n",
               ylim=c(0,max(density(data[class==1,k])$y*1.5)),
               xlim=range(data[,k]),
               main="",xlab=paste(k),ylab="")
          for(g in 1:G) lines(density(data[class==g,k]),col=g+2)
          legend("topleft",bty="n",paste(c(1:G)),col=3:(G+2),lty=1)
        }
        else{                                      ##Factor Diag
          barplot(table(data[,k],class))
        }}
      else{                                      ##NonDiag
        if((length(unique(data[,k]))>.3*nrow(data)) | (length(unique(data[,.k]))>.3*nrow(data)) ){
          
          plot(data[,c(.k,k)],xlab=paste(.k),ylab=paste(k))
          for(g in 1:G) points(data[class==g,c(.k,k)],col=g+2,pch=g)
         # if(!is.na(out)[1])  points(data[out,c(.k,k),drop=FALSE],pch=4,col=2)
        #  plot(data[,c(.k,k)],xlab=paste(.k),ylab=paste(k))
         #  for(g in 1:G) points(data[class==g,c(.k,k)],col=g+2,pch=g)
        }
        else{
          tt = table(data[,.k],data[,k],class)
          ttP = NULL
          for(g in 1:G) ttP <- rbind(ttP,tt[,,g])
          colnames(ttP) <- sort(unique(data[,k]))
          rownames(ttP) <- paste(rep(paste("G",1:G,sep=""),each=nrow(tt)),":",sort(unique(data[,.k])),sep="")
          image(1:nrow(ttP),1:ncol(ttP),xlab=paste("G:",.k),ylab=paste(k),
                ttP/sum(ttP),col=grey(c(length(ttP):0)/length(ttP)),axes=FALSE)
          axis(1,at=1:nrow(ttP),lab=rownames(ttP))
          axis(2,at=1:ncol(ttP),lab=colnames(ttP))
        }}
    }}
}

#my.pairs2(used[,1:6],used[,7],3,NA,TRUE);  mtext(side=3,outer=TRUE,cex=1.25,paste(info$Name,":",info$Roles,"N=",info$N,"(True Class)",sep=""))
#data = used[,1:6]; class = used[,7]; G=3; col=TRUE
##u2 <- NULL; 
#for(i in 1:info$N) u2[i] <- 4*(used[i,7]-1)+1*ifelse(used[i,5]==1,ifelse(used[i,6]==1,1,2),ifelse(used[i,6]==1,3,4))   
#class=u2

.pt <<- c(3,2,1,4,9,8,6,7)
my.pairs2 <- function(data,class,G,varType,vNames,out,col){
  K <- dim(data)[2]
  .pG <- length(table(class))
  .R <- .pG/G
  pg <- cbind(rep(1:G,each=.R),rep(1:.R,G),1:.pG)
  .u <- ifelse(.R==1,1,2)
  if(is.null(names)) names <- paste(1:K)
    if(col){ 
      u.col <- rainbow(.pG) 
      u.lwd <- rep(1,.R)
    }else{ 
      u.col <- grey(seq(.1,.75,length=.R)) 
      u.lwd <- seq(1,2,length=G)
    }
  par(mfrow=c(K,K),mar=c(2,2,0,0),oma=c(0,2,4,0))
  for(k in 1:K){
    for(.k in 1:K){
      if( k ==.k ){
        if(varType[k]=="N"){ ##Continuous Diag
          plot(density(data[class==1,k],bw=.5),type="n",
               ylim=c(0,max(density(data[class==1,k])$y*2)),
               xlim=range(data[,k]),
               main="",xlab="",ylab="")
          if(k==1 )  mtext(side=2,vNames[k],line=2,cex=1.25)
          if(.k==1 ) mtext(side=3,vNames[k],line=0,cex=1.25,las=1) 
          for(g in 1:.pG) lines(density(data[class==g,k]),col=u.col[pg[g,.u]],lwd=u.lwd[pg[g,1]])
          legend("topleft",bty="n",paste(c(1:G)),col=u.col,lty=1,cex=.5)
        }else{                                      ##Factor Diag
          barplot(table(data[,k],class))
          if( k==1 ) mtext(side=1,paste(k),line=0,cex=.75)
        }}else{                                      ##NonDiag
          x1 <- varType[.k] =="N"
          y1 <- varType[k]  =="N"
          if(x1 | y1){
           xrange <- range(data[,.k])
           yrange <- range(data[,k])
          if( !x1 )  xrange <- xrange+c(-.25,.25)
          if( !y1 )  yrange <- yrange+c(-.25,.25)           
            plot(data[,c(k,.k)],type="n",xlab="",ylab="",xlim=xrange,ylim=yrange)
            if(k==1)  mtext(side=3,vNames[.k],line=0,cex=1.25)
            if(.k==1 ) mtext(side=2,vNames[k],line=2,cex=1.25,las=0)             
            if(k < .k) for(g in 1:.pG) points(data[class==g,c(.k,k)],col=u.col[pg[g,.u]],pch=.pt[pg[g,1]])
            if(k > .k) for(g in .pG:1) points(data[class==g,c(.k,k)],col=u.col[pg[g,.u]],pch=.pt[pg[g,1]])
            if( !is.na(out)[1] ) points(data[out,c(k,.k),drop=FALSE],pch=.pt[pg[class[out],1]],col=u.col[pg[class[out],.u]],cex=2)
          }else{
            tt = table(data[,.k],data[,k],class)
            ttP = NULL
            for(g in 1:.pG) ttP <- rbind(ttP,tt[,,g])
            colnames(ttP) <- sort(unique(data[,k]))
            rownames(ttP) <- NULL
            #rownames(ttP) <- paste(rep(paste("G",1:G,sep=""),each=nrow(tt)),":",sort(unique(data[,.k])),sep="")
            image(1:nrow(ttP),1:ncol(ttP),xlab="",ylab="",
                  ttP/sum(ttP),col=grey(c(length(ttP):0)/length(ttP)),axes=FALSE)
            #axis(1,at=1:nrow(ttP),lab=rownames(ttP))
            axis(2,at=1:ncol(ttP),lab=colnames(ttP))
            if(k==1)  mtext(side=3,vNames[.k],line=0,cex=1.25)
            if(.k==1 ) mtext(side=2,vNames[k],line=2,cex=1.25,las=0)                        
          }}
    }}
}

#source("shaded.R")
.pt <<- c(3,2,1,4,9,8,6,7)
my.pairs3 <- function(data,class,G,varType,vNames,out,col){
  K <- dim(data)[2]
  .pG <- length(table(class))
  .R <- .pG/G
  pg <- cbind(rep(1:G,each=.R),rep(1:.R,G),1:.pG)
  .u <- ifelse(.R==1,1,2)
  if(is.null(names)) names <- paste(1:K)
  if(col){ 
    u.col <- rainbow(G) 
    u.lwd <- rep(3,length=G)
  }else{ 
    u.col <- grey(seq(.1,.75,length=.R)) 
    u.lwd <- seq(1,2,length=G)
  }
  par(mfrow=c(K,K),mar=c(2,2,0,0),oma=c(0,2,4,0))
  for(k in 1:K){
    for(.k in 1:K){
      if( k ==.k ){
        if(varType[k]=="N"){ ##Continuous Diag
          plot(density(data[class==1,k],bw=.5),type="n",
               ylim=c(0,max(density(data[class==1,k])$y*2)),
               xlim=range(data[,k]),
               main="",xlab="",ylab="")
          if(k==1 )  mtext(side=2,vNames[k],line=2,cex=1.25)
          if(.k==1 ) mtext(side=3,vNames[k],line=0,cex=1.25,las=1) 
          if(.u==1){
            for(g in 1:.pG) lines(density(data[class==g,k]),col=u.col[pg[g,1]],lwd=u.lwd[pg[g,1]])
          }
          if(.u==2){
            for(g in 1:.pG) lines(density(data[class==g,k]),col=u.col[pg[g,1]],lwd=u.lwd[pg[g,1]],lty=pg[g,.u])
          }
          legend("topleft",bty="n",paste(c(1:G)),col=u.col,lty=1,cex=.5)
        }else{                                      ##Factor Diag
          if(.u==1) barplot(t(table(data[,k],class)),col=u.col)
          if(.u==2) barplot(t(table(data[,k],class)),col=rep(u.col,each=.R),
                            density=rep(seq(-1,20,length=.pG/G)))
          if( k==1 ) mtext(side=1,paste(k),line=0,cex=.75)
        }}else{                                      ##NonDiag

          x1 <- varType[.k] =="N"
          y1 <- varType[k]  =="N"
          xrange <- range(data[,.k])
          yrange <- range(data[,k])
          if(sum(c(x1,y1))==2){
            plot(data[,c(k,.k)],type="n",xlab="",ylab="",xlim=xrange,ylim=yrange)
            if(k==1)  mtext(side=3,vNames[.k],line=0,cex=1.25)
            if(.k==1 ) mtext(side=2,vNames[k],line=2,cex=1.25,las=0)             
            if(k < .k) for(g in 1:.pG) points(data[class==g,c(.k,k)],col=u.col[pg[g,.u]],pch=.pt[pg[g,1]])
            if(k > .k) for(g in .pG:1) points(data[class==g,c(.k,k)],col=u.col[pg[g,.u]],pch=.pt[pg[g,1]])
            if( !is.na(out)[1] ) points(data[out,c(k,.k),drop=FALSE],pch=.pt[pg[class[out],1]],col=u.col[pg[class[out],.u]],cex=2)
          }
          if(sum(c(x1,y1))==1 & .u==1){
            if(k < .k){
            if( !x1 )  xrange <- xrange+c(-.25,.25)
            if( !y1 )  yrange <- yrange+c(-.25,.25)
             cats <- c(sapply(1:max(data[,.k]), function(w)
                    sapply(1:G, function(g) paste("l=",w,", g=",g,sep=""))))
             tBox1 <- paste("l=",data[,.k],", g=",class,sep="")
             tBox2 <- sapply(1:nrow(data), function(i) which(cats==tBox1[i]))
             useBox = factor(x = tBox2, 1:length(cats), labels = cats)

             boxplot(data[,k]~useBox,col=u.col)
            if(k==1)  mtext(side=3,vNames[.k],line=0,cex=1.25)
            if(.k==1 ) mtext(side=2,vNames[k],line=2,cex=1.25,las=0)             
          }else{ plot(1,1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n") }}
          if(sum(c(x1,y1))==1 & .u==2){
              if(k < .k){
                if( !x1 )  xrange <- xrange+c(-.25,.25)
                if( !y1 )  yrange <- yrange+c(-.25,.25)
                cats <- c(sapply(1:max(data[,.k]), function(w)
                  sapply(1:G, function(g) paste("l=",w,", g=",g,sep=""))))
                tBox1 <- paste("l=",data[,.k],", g=",pg[class,1],sep="")
                tBox2 <- sapply(1:nrow(data), function(i) which(cats==tBox1[i]))
                useBox = factor(x = tBox2, 1:length(cats), labels = cats)
                
                boxplot(data[,k]~useBox,col=rep(u.col,times=.u)) 
                if(k==1)  mtext(side=3,vNames[.k],line=0,cex=1.25)
                if(.k==1 ) mtext(side=2,vNames[k],line=2,cex=1.25,las=0)             
              }else{
                plot(1,1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
                 if(.k==1 ) mtext(side=2,vNames[k],line=2,cex=1.25)
          }}
          if(sum(c(x1,y1))==0){
            if(k < .k){
              tt = table(data[,.k],data[,k],class)
              ttP <- NULL
              for(g in 1:(.pG-1)) ttP = cbind(ttP,tt[,,g],rep(0,dim(tt)[1]))
              ttP <- t(cbind(ttP,tt[,,.pG]))
              image(1:nrow(ttP),1:ncol(ttP),xlab="",ylab="",
                    ttP/sum(ttP),col=grey(c(length(ttP):0)/length(ttP)),axes=FALSE) 
              abline(v=which(apply(ttP,1,sum)==0))
              ylab <- paste("w=",1:max(data[,.k]),sep="")
              xlab <- NULL
              xlab <- c(sapply(1:G, function(g)
                sapply(1:c(max(data[,k])+1), function(w)
                  ifelse(w<=max(data[,k]),paste("l=",w,", g=",g,sep="")," "))))
              xlab <- xlab[-length(xlab)]
            }else{
              tt = table(data[,k],data[,.k],class)
              ttP <- NULL
              for(g in 1:(.pG-1)) ttP = cbind(ttP,tt[,,g],rep(0,dim(tt)[1]))
              ttP <- t(cbind(ttP,tt[,,.pG]))
              image(1:nrow(ttP),1:ncol(ttP),xlab="",ylab="",
                    ttP/sum(ttP),col=grey(c(length(ttP):0)/length(ttP)),axes=FALSE)
             abline(v=which(apply(ttP,1,sum)==0))
              ylab <- paste("l=",1:max(data[,k]),sep="")
              xlab <- NULL
              xlab <- c(sapply(1:G, function(g)
                sapply(1:c(max(data[,.k])+1), function(w)
                  ifelse(w<=max(data[,.k]),paste("l=",w,", g=",g,sep="")," "))))
              xlab <- xlab[-length(xlab)]
            }
            axis(1,at=1:nrow(ttP),lab=xlab[1:nrow(ttP)],cex.axis=.8)
            axis(2,at=1:ncol(ttP),lab=ylab,cex.axis=.8)
            if(k==1)  mtext(side=3,vNames[.k],line=0,cex=1.25)
            if(.k==1 ) mtext(side=2,vNames[k],line=2,cex=1.25,las=0)                        
          }}
    }}
}


.lineN <<- c(c(1,4),c(1:4),c(1:4))
.colN  <<- c(rep(grey(.8),2),rep(grey(.6),4),rep(grey(.1),4))
plot.bic <- function(clust.obj,data.name,yrange,xrange,ylab,leg){  
  useg <- as.numeric(rownames(clust.obj$BIC))
  plot(clust.obj$BIC,type="n",bty="l",main="",ylim=yrange,xlim=c(min(useg),max(useg)+1),ylab=ylab,xaxt="n",xlab="G") 
  mtext(data.name,side=3,outer=FALSE,cex=1,line=1)
  mtext(paste("(G=",length(table(clust.obj$class)),",Model=",clust.obj$model,")"),side=3,outer=FALSE,
        cex=.9,line=0)
  DB <- dim(clust.obj$BIC)[2]
  for(s in 1:DB) lines(useg,clust.obj$BIC[,s],lty=.lineN[s],col=.colN[s],lwd=1.5)
  axis(1, at=c(xrange[1]:xrange[2]),labels=paste(xrange[1]:xrange[2]))
  if(leg==TRUE) legend("topright",colnames(clust.obj$BIC),col=.colN[1:DB],lty=.lineN[1:DB],
                       cex=.75,bty="n",lwd=1.5)
  
}

modName <- list()
modName[[1]] <- .mclustModel
modName[[2]] <- .mixdModel
modName[[3]] <- c("Sphr","Diag","Ellip")
modName[[4]] <- paste(rep(modName[[3]],each=4),"-",.mixdModel,sep="")
.lineS <<- rep(c(1:3),each=4)
.pchS  <<- rep(c(1:3),each=4)
.colS  <<- rep(grey(c(.05,.3,.6,.9)),3)
.lineA <<- rep(c(1:4),10)
.pchA  <<- rep(c(1:4),10)
.colA  <<- rep(rainbow(10),each=4)

plot.bic2 <- function(clust.obj,data.name,yrange,ylab,leg,sum,all){  
  useg <- as.numeric(rownames(clust.obj$BIC))
  summ <- mod <- tmpBIC <- list(); DB <- NULL
  mod <- matrix(NA,4,ncol(clust.obj$BIC))
  nor <- substr(dimnames(clust.obj$BIC)[[2]],1,3)
  mix <- substr(dimnames(clust.obj$BIC)[[2]],5,6)
  if(length(grep("X",nor))>0){ 
	.useMclust <- .mclustUni; nMod <- 2
	}else{	     .useMclust <- .mclustModel; nMod <- 4 }
  for(i in 1:length(.mclustModel)) mod[1,which(nor==.useMclust[i])] <- i
  for(i in 1:length(.mixdModel))   mod[2,which(mix==.mixdModel[i])]   <- i
  mod[3,which(nor %in% .sphrModel)] <- 1
  mod[3,which(nor %in% .diagModel)] <- 2
  mod[3,which(nor %in% .ellpModel)] <- 3
  mod[4,] <- mod[2,]+4*(mod[3,]-1)
  for(i in 1:nMod) DB[i] <- max(mod[i,],na.rm=TRUE)
  for(i in 1:nMod){
    tmpBIC[[i]] <- matrix(NA,length(useg),DB[i])
    for(j in 1:DB[i]) tmpBIC[[i]][,j] <- rowMeans(clust.obj$BIC[,mod[i,]==j,drop=FALSE],na.rm=TRUE)
  } 
  xrange <- range(useg)
  .useg  <- paste(useg)
  if(leg){ xrange[2] <- xrange[2]+1; .useg <- c(.useg,"") } 
  if(sum){
  plot(clust.obj$BIC,type="n",main="",ylim=range(tmpBIC[[nMod]],na.rm=TRUE),
       xlim=xrange,ylab=ylab,xlab="# Groups",xaxt="n",bty="L") 
  for(s in 1:ncol(tmpBIC[[4]])){
    dt <- data.frame(A=useg,B=tmpBIC[[4]][,s])
    xy <- list(x=dt$A,y=dt$B)
    if(sum(!is.na(dt$B))>1)  xy <- approx(x=dt$A,y=dt$B,xout=seq(xrange[1],xrange[2],by=.1))
    points(xy$x,xy$y,pch=.pchS[s],col=.colS[s],cex=.75)
  }
  if(leg==TRUE) legend(xrange[2]-.45,yrange[2]-10,modName[[4]],pch=.pchS,col=.colS,bty="n",cex=.75)  
  }
  if(all){
    plot(clust.obj$BIC,type="n",main="",ylim=range(tmpBIC[[nMod]],na.rm=TRUE),
         xlim=xrange,ylab=ylab,xlab="# Groups",xaxt="n",bty="L") 
    for(s in 1:ncol(clust.obj$BIC)){
     dt <- data.frame(A=useg,B=clust.obj$BIC[,s])
     xy <- list(x=dt$A,y=dt$B)
     if(sum(!is.na(dt$B))>1)  xy <- approx(x=dt$A,y=dt$B,xout=seq(xrange[1],xrange[2],by=.1))
     points(xy$x,xy$y,pch=.pchA[s],col=.colA[s],cex=.75)
    }
    if(leg==TRUE){
      legend(xrange[2]-.95,yrange[2]-10,colnames(clust.obj$BIC)[1:20],pch=.pchA[1:20],col=.colA[1:20],bty="n",cex=.75)  
      legend(xrange[2]-.45,yrange[2]-10,colnames(clust.obj$BIC)[21:40],pch=.pchA[21:40],col=.colA[21:40],bty="n",cex=.75)  
  }  }
  mtext(data.name,side=3,outer=FALSE,cex=1,line=1)
  mtext(paste("(G=",length(table(clust.obj$class)),",Model=",clust.obj$model,")"),side=3,outer=FALSE,
        cex=.9,line=0)
  axis(1, at=c(xrange[1]:xrange[2]),labels=.useg)
  }
#plot.bic2(clust.obj,"VII.RJ",range(clust.obj$BIC,na.rm=TRUE),"BIC",T,T,F)
#par(mfrow=c(1,2))
#plot.bic2(clust.obj,"VII.RJ - Summarized over Normal Structures",range(clust.obj$BIC,na.rm=TRUE),"BIC",T,T,F)
#plot.bic2(clust.obj,"VII.RJ - All Models",range(clust.obj$BIC,na.rm=TRUE),"BIC",T,F,T)

