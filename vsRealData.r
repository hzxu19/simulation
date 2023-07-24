if(getwd()!="/u/kevans/Research Files/tmpVS2"){
if(getwd()=="C:/Documents and Settings/kevans/My Documents"){
 setwd("C://Documents and Settings//kevans//My Documents//My Dropbox//tmp//tmpVS2//")
}else{
 setwd("C://Users//K80//Documents//My Dropbox//tmp//newVS//")
}}

#library(Rmixmod)
source("vsLoad.R")
source("vsFunction.R")
whichReg <- stepReg
source("vsFunction.R")

##---------------------------------------------
##IrisCheck
##----------------------
##Recreate Raftery & Dean
.fullData = iris[,1:4]
.varTypes = rep("N",4)
.startVar=1:4
dean.iris <- stepSelect(.startVar,2,9,.mclustModel,.fullData,.varTypes,
		which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=TRUE,"Dean")
## dean.iris <- Correctly matches variables 2:4
maug.iris <- stepSelect(.startVar,2,9,.mclustModel,.fullData,.varTypes,
		which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=FALSE,"Maug")
dean.iris$SW$allRoles; dean.iris$SRUW$allRoles
maug.iris$SW$allRoles; maug.iris$SRUW$allRoles
## maug.iris matches with dean.iris when allowed to cover 2-9, all models

.cat <- factor(c(rep(1,100),rep(2,50)))
.fullData = data.frame(cbind(iris[,1:4],.cat)); dim(.fullData)
.varTypes= c(rep("N",4),"C")
.startVar=1:5
dean.irisC <- stepSelect(.startVar,2,9,.mclustModel,.fullData,.varTypes,
		which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=TRUE,"Dean")
## matches with dean.iris, puts "C" into "W" or "U"
maug.irisC <- stepSelect(.startVar,2,9,.mclustModel,.fullData,.varTypes,
		which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=FALSE,"Maug")
dean.irisC$SW$allRoles; dean.irisC$SRUW$allRoles
maug.irisC$SW$allRoles; maug.irisC$SRUW$allRoles
## 
##--------------------
##CrabCheck
##------
.cat <- factor(c(rep(1,100),rep(2,100)))
.fullData = data.frame(cbind(crabs[,4:8],.cat)); dim(.fullData)
.fullData = crabs[,4:8]
.varTypes= c(rep("N",5))#,"C")
.startVar=1:5#6
dean.crabC52 <- stepSelect(.startVar,2,5,.mclustModel,.fullData,.varTypes,
		NULL,.kp,.tol,.max.iter,.mc=TRUE,.intermed=TRUE,"Dean")
##
maug.crabC <- stepSelect(.startVar,2,9,.mclustModel,.fullData,.varTypes,
		.which.pri,.kp,.tol,.max.iter,.mc=TRUE,.intermed=FALSE,"Dean")
##
dean.crabC$SW$allRoles; dean.crabC$SRUW$allRoles
maug.crabC$SW$allRoles; maug.crabC$SRUW$allRoles

source("./clusterfunctions/mod_clvarselnosampgr.r")
chk1 <- my.clustvarsel(d[[2]][,1:4], 2,5, emModels1 = c("E", "V"), emModels2 = c("EII", 
                                                           "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV"), 
          samp = FALSE, sampsize = 2000, allow.EEE = FALSE, forcetwo = TRUE, 
          search = "greedy", upper = 0, lower = -10, itermax = 100, prior = NULL,
          prior2=NULL) 
source("./clusterfunctions/clustvarsel.r")

chk2 <- clustvarsel(crabs[,4:8],2:5,emModels1=c("E","V"),emModels2=c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),
       
#Prostate BELOW
##----------------------------------------
##ProstateCheck
##------------------------------------------
pp = read.table(file="./data/prostate.csv",header=TRUE,sep=",")
##NOTE: Might consider re-adding in observations which are missing some variables
##      if variables turn out to not be important for clustering
##      Might add something to the paper about how missing data could also be handled
##      Why do all of the analyses mention how SG could be treated as continuous vs discrete
##       but say nothing about DBP which only takes on one more level?
##      Proposed reparametrization for a number of variables
##      .O suffix generally indicates original levels of variables - But anything with <=12 levels is treated as discrete
##      .M suffix indicates redefined levels for discrete variables (i.e. anything with <=12 levels is candidate for redefinition
comp.pp = na.omit(pp); dim(comp.pp)
comp.pp = comp.pp[-which(comp.pp$ekg==""),]
comp.pp[["stage2"]] <- stage2 <- as.numeric(comp.pp[["stage"]])
comp.pp[["status1"]] <- status1 <- sapply(1:nrow(comp.pp),function(i)
	ifelse(as.character(comp.pp[["status"]][i])=="alive",1,
	ifelse(as.character(comp.pp[["status"]][i])=="dead - prostatic ca",2,3)))
newOutcome1 <- c("Alive","Dead - Pros","Dead - All Other")
comp.pp[["status2"]] <- status2 <- sapply(1:nrow(comp.pp),function(i)
	ifelse(as.character(comp.pp[["status"]][i])=="alive",1,
	ifelse(as.character(comp.pp[["status"]][i])=="dead - prostatic ca",2,
	ifelse(as.character(comp.pp[["status"]][i])=="dead - heart or vascular"
		 | as.character(comp.pp[["status"]][i])=="dead - pulmonary embolus"
             | as.character(comp.pp[["status"]][i])=="dead - cerebrovascular",3,4))))
newOutcome2 <- c("Alive","Dead - Pros","Dead - Card","Dead - Oth")
comp.pp[["status3"]] <- status3 <- sapply(1:nrow(comp.pp),function(i)
	ifelse(as.character(comp.pp[["status"]][i])=="alive",1,2))
newOutcome3 <- c("Alive","Dead")
comp.pp[["status4"]] <- status4 <- sapply(1:nrow(comp.pp),function(i)
  ifelse(as.character(comp.pp[["status"]][i])=="dead - prostatic ca",1,2))
newOutcome4 <- c("Dead - Prostate","Other")
comp.pp[["rx2"]] <- rx2 <- sapply(1:nrow(comp.pp), function(i)
			 ifelse(as.character(comp.pp[["rx"]][i])=="placebo",1,2))
comp.pp[["rx4"]] <- rx4 <- as.numeric(comp.pp[["rx"]])
comp.pp[["rx3"]] <- rx3 <- sapply(1:nrow(comp.pp), function(i)
                         ifelse(rx4[i]<=2,1,2))

trueClass <- cbind(stage2,status1,status2,status3,status4,rx2,rx3,rx4)

useCol = c("age","wt","pf","hx","sbp","dbp","ekg","hg","sz","sg","ap","bm")
useP = comp.pp[,useCol]
datP = data.frame(useP)
frx2 <- factor(rx2)
frx3 <- factor(rx3)
frx4 <- factor(rx4)
fstatus1 <- factor(status1)
fstatus2 <- factor(status2)
fstatus3 <- factor(status3)
fstatus4 <- factor(status4)

sq.sz <- sqrt(datP$sz)
lg.ap <- log(datP$ap)
pfNew <- rep(0,nrow(datP))
pfNew[which(datP$pf=="normal activity")] <- 1
ekgNew <- datP$ekg
ekgNew <- gsub("recent MI","MI",ekgNew)
ekgNew <- gsub("old MI","MI",ekgNew)
dbpNew <- rep(0,nrow(datP))
for(i in 1:nrow(datP)) dbpNew[i] <- ifelse(datP[i,"dbp"]<=6,1,
                                     ifelse(datP[i,"dbp"]<=9,datP[i,"dbp"]-5,
                                      ifelse(datP[i,"dbp"]>9,5,999)))
sgNew <-  rep(0,nrow(datP))
for(i in 1:nrow(datP)) sgNew[i] <- ifelse(datP[i,"sg"]<=7,1,
                                     ifelse(datP[i,"sg"]<=13,datP[i,"sg"]-6,
                                       ifelse(datP[i,"sg"]>13,8,999)))
table(datP[,"dbp"],dbpNew)
table(datP[,"sg"],  sgNew)
pfNew.N <- as.numeric(pfNew)
pf.N    <- as.numeric(datP$pf)
hx.N    <- as.numeric(datP$hx)
ekg.N   <- as.numeric(datP$ekg)
ekgNew.N<- as.numeric(factor(as.character(ekgNew)))
bm.N    <- as.numeric(datP$bm)
dbpNew.N<- as.numeric(dbpNew)
dbp.N   <- as.numeric(datP$dbp)
sgNew.N <- as.numeric(sgNew)
sg.N    <- as.numeric(datP$sg)

pfNew.F <- factor(pfNew)
pf.F    <- factor(datP$pf)
hx.F    <- factor(datP$hx)
ekg.F   <- factor(datP$ekg)
ekgNew.F<- factor(ekgNew)
bm.F    <- factor(datP$bm)
dbpNew.F<- factor(dbpNew)
dbp.F   <- factor(datP$dbp)
sgNew.F <- factor(sgNew)
sg.F    <- factor(datP$sg)

datPlot <- cbind(useP[,!(names(useP) %in% c("sz","ap","pf","ekg","hx","bm","dbp","sg"))],
			sq.sz,lg.ap,pf.N,pfNew.N,hx.N,ekg.N,ekgNew.N,bm.N,dbp.N,dbpNew.N,sg.N,sgNew.N)
datClst <- cbind(useP[,!(names(useP) %in% c("sz","ap","pf","ekg","hx","bm","dbp","sg"))],
			sq.sz,lg.ap,pf.F,pfNew.N,pfNew.F,hx.F,ekg.F,
                 ekgNew.F,bm.F,dbp.N,dbpNew.F,sg.N,sgNew.F)
nVarsO <- c("age","wt","sbp","dbp.N","hg","sq.sz","lg.ap","sg.N")
nVarsM <- c("age","wt","sbp","hg","sq.sz","lg.ap")
cVarsPO <- c("pf.N","hx.N","ekg.N","bm.N")
cVarsCO <- c("pf.F","hx.F","ekg.F","bm.F")
cVarsPM <- c("pfNew.N","hx.N","ekgNew.N","bm.N","dbpNew.N","sgNew.N")
cVarsCM <- c("pfNew.F","hx.F","ekgNew.F","bm.F","dbpNew.F","sgNew.F")

datPlot.O <- datPlot[,c(nVarsO,cVarsPO)]
datPlot.M <- datPlot[,c(nVarsM,cVarsPM)]
datClst.O <- datClst[,c(nVarsO,cVarsCO)]
datClst.M <- datClst[,c(nVarsM,cVarsCM)]
vT.O <- sapply(1:ncol(datClst.O),function(i)ifelse(colnames(datClst.O)[i] %in% nVarsO,"N","C"))
vT.M <- sapply(1:ncol(datClst.M),function(i)ifelse(colnames(datClst.M)[i] %in% nVarsM,"N","C"))
vT.Op <- sapply(1:ncol(datPlot.O),function(i)ifelse(colnames(datPlot.O)[i] %in% nVarsO,"N","C"))
vT.Mp <- sapply(1:ncol(datPlot.M),function(i)ifelse(colnames(datPlot.M)[i] %in% nVarsM,"N","C"))

##-----------------------
locIndep.O <- whichClust(NULL,1:ncol(datClst.O),"Excl",datClst.O,1,4,vT.O,
                       c("VVI.II"),"MOD",kp,tol,max.iter,TRUE)
locIndep.M <- whichClust(NULL,1:ncol(datClst.M),"Excl",datClst.M,1,4,vT.M,
                       c("VVI.II"),"MOD",kp,tol,max.iter,TRUE)
table(locIndep.O$class,locIndep.M$class)
##-----------------------
tryMods   <-  paste(rep(.mclustModel,each=4),
                  c(".II",".IJ",".CI",".CJ"),sep="") ##Range of Models to consider
bstModel.O <- whichClust(NULL,1:ncol(datClst.O),"Excl",datClst.O,2,4,vT.O,
				tryMods,"MOD",kp,tol,max.iter,TRUE)
bstModel.M <- whichClust(NULL,1:ncol(datClst.M),"Excl",datClst.M,2,4,vT.M,
				tryMods,"MOD",kp,tol,max.iter,TRUE)

table(bstModel.O$class,bstModel.M$class)
table(locIndep.O$class,bstModel.O$class)
table(locIndep.M$class,bstModel.M$class)
##-----------------------
tryMods  <- paste(rep(.mclustModel,each=4),
                  c(".II",".CI"),sep="") ##Range of Models to consider
deanVS.O <- stepSelect(c(1:ncol(datClst.O)),2,4,tryMods,datClst.O,vT.O,
			"MOD",kp,tol,max.iter,TRUE,TRUE,"Dean")
deanVS.M <- stepSelect(c(1:ncol(datClst.M)),2,4,tryMods,datClst.M,vT.M,
			"MOD",kp,tol,max.iter,TRUE,TRUE,"Dean")

##-----------------------
tryMods   <- paste(rep(.mclustModel,each=4),
                  c(".II",".IJ",".CI",".CJ"),sep="") ##Range of Models to consider
deanVS2.O <- stepSelect(c(1:ncol(datClst.O)),2,4,tryMods,datClst.O,vT.O,
			"MOD",kp,tol,max.iter,TRUE,TRUE,"Dean")
deanVS2.M <- stepSelect(c(1:ncol(datClst.M)),2,4,tryMods,datClst.M,vT.M,
			"MOD",kp,tol,max.iter,TRUE,TRUE,"Dean")
save <- FALSE

deanVS2VVV.M <- whichClust(NULL,1:2,
                           "Excl",datClst.M[,c(1,6)],3,3,c("N","N"),c("VVV"),"MOD",kp,tol,max.iter,TRUE)

colnames(datClst.O)[deanVS.O$SW$allRoles=="S"]
colnames(datClst.M)[deanVS.M$SW$allRoles=="S"]
colnames(datClst.O)[deanVS2.O$SW$allRoles=="S"]
colnames(datClst.M)[deanVS2.M$SW$allRoles=="S"]
deanVS.O$SW$fBIC$clust$G
deanVS.O$SW$fBIC$clust$par$mean
my.pairs2(datClst.O[,colnames(datClst.O)[deanVS.O$SW$allRoles=="S"]],
          deanVS.O$SW$fBIC$clust$class,
          deanVS.O$SW$fBIC$clust$G,
          vT.O[deanVS.O$SW$allRoles=="S"],
          colnames(datClst.O)[deanVS.O$SW$allRoles=="S"],NA,TRUE)
table(deanVS.O$SW$fBIC$clust$class)

if(save){
  saveList <- list(locIndep.O,locIndep.M,bstModel.O,bstModel.M,deanVS.O,
                   deanVS.M,deanVS2.O,deanVS2.M,deanVS2VVV.M)
save(saveList,file="saveProstate.rdata")
}
load <- TRUE
if(load){
  load("saveProstate.rdata")
  locIndep.O <- saveList[[1]]
  locIndep.M <- saveList[[2]]
  bstModel.O <- saveList[[3]]
  bstModel.M <- saveList[[4]]
  deanVS.O   <- saveList[[5]]
  deanVS.M   <- saveList[[6]]
  deanVS2.O  <- saveList[[7]]
  deanVS2.M  <- saveList[[8]]
  deanVS2VVV.M<- saveList[[9]]
}

for(i in 1:4) print(paste(saveList[[i]]$G,saveList[[i]]$model),sep=":")
for(i in c(5,7)) print(paste(saveList[[i]]$SW$fBIC$clust$G,saveList[[i]]$SW$fBIC$clust$model,
                          paste(colnames(datClst.O)[which(saveList[[i]]$SW$allRoles=="S")],collapse=", "),sep=":"))
for(i in c(6,8)) print(paste(saveList[[i]]$SW$fBIC$clust$G,saveList[[i]]$SW$fBIC$clust$model,
                          paste(colnames(datClst.M)[which(saveList[[i]]$SW$allRoles=="S")],collapse=", "),sep=":"))
##----------------------------------------
##Plotting
graphics.off()

windows(height=3,width=5.5)#,family="Times")
#pdf("C://Users//K80//Documents//My Dropbox//tmp//text//paperProstatePlot.pdf",height=3,width=5.5)
pdf("C://Users//K80//Documents//My Dropbox//tmp//text//thesisProstatePlot.pdf",height=3,width=5.5)
                    x <- "lg.ap"; xN <- "Log(AP)"
y0 <- c("age","dbp.N"); y0N <- c("Age","DBP (Orig)")
y1 <- c("age","dbpNew.N"); y1N <- c("Age","DBP (Mod)")
par(mfrow=c(1,1),mar=c(2,2,.5,4),oma=c(1.5,1.5,.5,1.5))#,pin=c(4,2))
                    myPCH <- c(21,22,24); 
                    myBG <- c(grey(.7),NA); 
                    myCol <- rep(1,3) #grey(c(.5,.2,.7)
useDat <- datPlot.M
useClust <- deanVS.M$SW$fBIC$clust
bgCol <- cbind(useClust$class,comp.pp$status4)
yy <- 1
  plot(useDat[,c(x,y0[yy])],type="n",bty="l",xlab="",ylab="",xlim=c(-2.5,8))
  mtext(side=1,paste(xN),cex=1,line=2); mtext(side=2,paste(y0N[yy]),cex=1,line=2)
  for(i in 1:nrow(useDat)) points(useDat[i,c(x,y0[yy])],
                                pch=myPCH[bgCol[i,1]],#col=myCol[bgCol[i,1]]),
                                bg=myBG[bgCol[i,2]])
abline(v=0,lty=3,lwd=2)
                    
par(xpd=NA)
 legend(7,77,bty="n",
         c("Filled characters","represent prostate","related death"),pch=rep(NA,3),cex=.8)
par(xpd=NA)
legend(7,91,bty="n",c("Clusters","Older, Stage 3","Older, Stage 4","Younger"),
       pch=c(NA,myPCH),col=c(NA,myCol),cex=.8,
       y.intersp=.9)

par(xpd=NA) 
legend(1,95,bty="n",c("Log(AP)=0 used to clinically","diagnose Stage 3, 4"),
       lty=c(3,NA),lwd=c(2,NA),cex=.8, y.intersp=.8)

dev.off(); dev.off(); graphics.off()

source("../outFiles/outlierID2.R")
source("../outFiles/newMetric.R")
                    
eig <- outlierID2(datClst.M[,c("lg.ap","age")],3,3,TRUE)
cls <- eig$nullClust$class
use <- order(cls)
mmm <- new.metric2(eig$min.eigs,cls,1,5)
ind <- out.ind(eig$min.eigs,cls,mmm$met)
my.pairs2(datClst.M[,c("lg.ap","age")],cls,max(cls),c("N","N"),c("Log(AP)","Age"),NA,TRUE)
my.pairs2(datClst.M[,c("lg.ap","age")],cls,max(cls),c("N","N"),c("Log(AP)","Age"),which(ind[,2]==1),TRUE)
look <- cbind(datClst.M[ind[,2]==1,c("lg.ap","age")],eig$min.eigs[ind[,2]==1],cls[ind[,2]==1])
look[order(look[,4]),]
round(eig$nullClust$par$mean,2)
cbind(datClst.M[ind[,2]==1,][order(look[,4]),],look[order(look[,4]),4])


trxStatus <- cbind(frx4,fstatus2,eig$nullClust$class)
table(trxStatus[,1],trxStatus[,2],trxStatus[,3])
trxStatus[ind[,2]==1,][order(look[,3]),]

for(g in 1:3) print(summary(datClst.M[eig$nullClust$class==g,]))
                    
pdf("prostateEigsPresentation.pdf")
plot(eig$min.eigs[use],pch=c(3,2,1)[cls[use]],col=rainbow(3)[cls[use]],bty="l",ylab="",xlab="")
for(m in 2){   
  for(g in 1:max(cls)) if(!is.na(mmm$met[m,g])) lines(which(cls[use]==g),
                                                      rep(mmm$met[m,g],sum(cls==g)),lwd=2,
                                                      col=rainbow(3)[g],lty=c(1,2,4,3)[m])
}
mtext(side=1,line=2,"Observation Index",cex=1.5)
mtext(side=2,line=2,"Eig",cex=1.5)
dev.off()

x <- "lg.ap"; xN <- "Log(AP)"
y0 <- c("age","dbp.N"); y0N <- c("Age","DBP (Orig)")
y1 <- c("age","dbpNew.N"); y1N <- c("Age","DBP (Mod)")
for(i in 1:2){
pdf(paste("prostatePresentation",i,".pdf",sep=""),height=4,width=6)
par(mfrow=c(1,1),mar=c(2,2,.5,2),oma=c(1.5,1.5,.5,1.5))#,pin=c(4,2))
useDat <- datPlot.M
useClust <- deanVS.M$SW$fBIC$clust
yy <- 1
plot(useDat[,c(x,y0[yy])],type="n",bty="l",xlab="",ylab="",ylim=c(45,95),xlim=c(-2.5,8))
mtext(side=1,paste(xN),cex=1,line=2); mtext(side=2,paste(y0N[yy]),cex=1,line=2)
for(g in 1:useClust$G) points(useDat[useClust$class==g,c(x,y0[yy])],
                              pch=c(3,2,1)[g],col=rainbow(3)[g])
if(i==2){
  for(g in 1:useClust$G) points(useDat[useClust$class==g & ind[,2]==1,c(x,y0[yy]),drop=FALSE],
                              pch=c(3,2,1)[g],col=rainbow(3)[g],cex=3)
}
abline(v=0,lty=3,lwd=2)
par(xpd=NA)
legend(7,91,bty="n",c("Clusters","Older, Stage 3","Older, Stage 4","Younger"),
       pch=c(NA,2,3,1),col=c(NA,rainbow(3)[c(2,1,3)]),cex=.8,
       y.intersp=.75)

par(xpd=NA) 
legend(1,97,bty="n",c("Log(AP)=0 used to clinically diagnose Stage 3, 4",""),
       lty=c(3,NA),lwd=c(2,NA),cex=.8, y.intersp=.75)

dev.off(); #dev.off()
}
pdf("thesis_Prostate.pdf",height=4,width=6)
                    x <- "lg.ap"; xN <- "Log(AP)"
                    y0 <- c("age","dbp.N"); y0N <- c("Age","DBP (Orig)")
                    y1 <- c("age","dbpNew.N"); y1N <- c("Age","DBP (Mod)")
                    par(mfrow=c(1,2),mar=c(4,4,1,0),oma=c(0,0,1,1))
                       
                    useDat <- datPlot.M
                    useClust <- locIndep.M
                    for(yy in 1){
                      plot(useDat[,c(x,y1[yy])],type="n",bty="l",xlab="",ylab="")
                      mtext(side=1,paste(xN),cex=1.5,line=2.5); 
                      mtext(side=2,paste(y1N[yy]),cex=1.5,line=2)
                      mtext(side=3,"All Variables",cex=1.5)
                      for(g in 1:useClust$G) points(useDat[useClust$class==g,c(x,y1[yy])],
                                                    pch=c(3,2)[g],col=rainbow(useClust$G)[g])
                    }
                    
                    
                    useDat <- datPlot.M
                    useClust <- deanVS.M$SW$fBIC$clust
                    yy <- 1
                    plot(useDat[,c(x,y0[yy])],type="n",bty="l",xlab="",ylab="",ylim=c(45,95),xlim=c(-2.5,8))
                    mtext(side=1,paste(xN),cex=1.5,line=2.5); 
                    mtext(side=2,paste(y1N[yy]),cex=1.5,line=2)
                    mtext(side=3,"Variable Selection",cex=1.5)
                    for(g in 1:useClust$G) points(useDat[useClust$class==g,c(x,y0[yy])],
                                                  pch=c(3,2,1)[g],col=rainbow(3)[g])
dev.off()                    
##--------------------------------------------------------------

methodNames <- c("Original Variables, No Grouping","Modified Variables, No Grouping",
                 "Local Independence Clustering - Original Variables","Local Independence Clustering - Modified Variables",
                 "Best Model - All Original Variables","Best Model - All Modified Variables",
                 "VS - Restricted Models (no J) - All Original Variables","VS - Restricted Models (no J) - All Modified Variables",
                 "VS - All Models - All Original Variables","VS - All Models - All Modified Variables")


pdf("csProstatePlot1.pdf",height=4,width=6)
x <- "lg.ap"; xN <- "Log(AP)"
y0 <- c("age","dbp.N"); y0N <- c("Age","DBP (Orig)")
y1 <- c("age","dbpNew.N"); y1N <- c("Age","DBP (Mod)")
par(mfrow=c(1,2),mar=c(4,4,1,0),oma=c(0,0,3,1))
useDat <- datPlot.O
useClust <- locIndep.O
for(yy in 2:1){
  plot(useDat[,c(x,y0[yy])],type="n",bty="l",xlab="",ylab="")
  mtext(side=1,paste(xN),cex=2,line=2.5); mtext(side=2,paste(y0N[yy]),cex=2,line=2)
  for(g in 1:useClust$G) points(useDat[useClust$class==g,c(x,y0[yy])],col=rainbow(useClust$G)[g])
}
mtext(side=3,"Original Data",line=2,outer=TRUE,cex=1.5)
mtext(side=3,"All variables (Local Independence)",line=1,cex=1.5,outer=TRUE)
mtext(side=3,"G=2, M=VVI.II",line=0,outer=TRUE,cex=1.5)

useDat <- datPlot.M
useClust <- locIndep.M
for(yy in 2:1){
  plot(useDat[,c(x,y1[yy])],type="n",bty="l",xlab="",ylab="")
  mtext(side=1,paste(xN),cex=2,line=2.5); mtext(side=2,paste(y1N[yy]),cex=2,line=2)
  for(g in 1:useClust$G) points(useDat[useClust$class==g,c(x,y1[yy])],col=rainbow(useClust$G)[g])
}
mtext(side=3,"Modified Data",line=2,outer=TRUE,cex=1.5)
mtext(side=3,"All variables (Local Independence)",line=1,cex=1.5,outer=TRUE)
mtext(side=3,"G=2, M=VVI.II",line=0,cex=1.5,outer=TRUE)
dev.off()


pdf("csProstatePlot2.pdf",height=4,width=6)
x <- "lg.ap"; xN <- "Log(AP)"
y0 <- c("age","dbp.N"); y0N <- c("Age","DBP (Orig)")
y1 <- c("age","dbpNew.N"); y1N <- c("Age","DBP (Mod)")
par(mfrow=c(1,2),mar=c(4,4,.5,0),oma=c(0,0,3,1))

useDat <- datPlot.O
useClust <- deanVS2.O$SW$fBIC$clust
for(yy in 2:1){
  plot(useDat[,c(x,y0[yy])],type="n",bty="l",xlab="",ylab="")
  mtext(side=1,paste(xN),cex=2,line=2.5); mtext(side=2,paste(y0N[yy]),cex=2,line=2)
  for(g in 1:useClust$G) points(useDat[useClust$class==g,c(x,y0[yy])],col=rainbow(useClust$G)[g])
}
mtext(side=3,paste("Original Data"),line=2,outer=TRUE,cex=1.5)
mtext(side=3,"Variable Selection: DBP and Log(AP)",line=1,cex=1.5,outer=TRUE)
mtext(side=3,paste("G=",useClust$G,", M=",useClust$model,sep=""),line=0,cex=1.5,outer=TRUE)

useDat <- datPlot.M
useClust <- deanVS2.M$SW$fBIC$clust
for(yy in 2:1){
  plot(useDat[,c(x,y1[yy])],type="n",bty="l",xlab="",ylab="")
  mtext(side=1,paste(xN),cex=2,line=2.5); mtext(side=2,paste(y1N[yy]),cex=2,line=2)
  for(g in 1:useClust$G) points(useDat[useClust$class==g,c(x,y1[yy])],col=rainbow(useClust$G)[g])
}
mtext(side=3,paste("Modified Data"),line=2,outer=TRUE,cex=1.5)
mtext(side=3,"Variable Selection: Age and Log(AP)",line=1,cex=1.5,outer=TRUE)
mtext(side=3,paste("G=",useClust$G,", M=",useClust$model,sep=""),line=0,cex=1.5,outer=TRUE)
dev.off()


pdf("useProstatePlots.pdf")
my.pairs2(datPlot.O[,which(deanVS2.O$SW$allRoles=="S")],deanVS.O$SW$fBIC$clust$class,deanVS.O$SW$fBIC$clust$G,
          vT.Op,colnames(datPlot.O)[which(deanVS.O$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Variable Selection - Original Data",cex=1)
my.pairs2(datPlot.M[,which(deanVS2.M$SW$allRoles=="S")],deanVS.M$SW$fBIC$clust$class,deanVS.M$SW$fBIC$clust$G,
          vT.Mp,colnames(datPlot.M)[which(deanVS.M$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Variable Selection - Modified Data",cex=1)
my.pairs2(datPlot.M[,which(deanVS2.M$SW$allRoles=="S")],deanVS2VVV.M$class,
          deanVS2VVV.M$G, vT.Mp,colnames(datPlot.M)[which(deanVS.M$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Variable Selection - Modified Data (VVV structure)",cex=1)

par(mfrow=c(1,2),mar=c(4,4,2,0),oma=c(0,0,2,1))
.grey <- grey(seq(.75,.1,length=deanVS2.M$SW$fBIC$clust$G)) 
plot(datPlot.M[,c("age","lg.ap")],type="n",xlab="Age",ylab="Log(AP)",bty="l")
mtext(side=3,line=1,"Variable Selection")
mtext(side=3,line=0,paste("Modified Variables: G=",deanVS2.M$SW$fBIC$clust$G,
                         "M=",deanVS2.M$SW$fBIC$clust$model),cex=.75)
for(g in 1:deanVS2.M$SW$fBIC$clust$G) points(
  datPlot.O[deanVS2.M$SW$fBIC$clust$class==g,c("age","lg.ap")],pch=.pt[g],col=.grey[g])

.grey <- grey(seq(.75,.1,length=deanVS.O$SW$fBIC$clust$G)) 
plot(datPlot.O[,c("dbp.N","lg.ap")],type="n",ylab="",xlab="DBP",bty="l")
mtext(side=3,line=1,"Variable Selection")
mtext(side=3,line=0,paste("Original Variables G=",deanVS2.M$SW$fBIC$clust$G,
      "M=",deanVS2.M$SW$fBIC$clust$model),cex=.75)
 for(g in 1:deanVS.O$SW$fBIC$clust$G) points(
   datPlot.O[deanVS.O$SW$fBIC$clust$class==g,c("dbp.N","lg.ap")],pch=.pt[c(3,1,2,4)][g],col=.grey[g])
par(mfrow=c(1,1))
plot.bic(saveList[[i]]$SW$fBIC$clust,methodNames[i+2],r,c(2,4),"BIC",TRUE)
dev.off()


pdf("prostatePlots.pdf")
my.pairs2(datPlot.O,rep(1,nrow(datPlot.O)),1,vT.Op,colnames(datPlot.O),NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Original Variables, No Grouping",cex=2)
my.pairs2(datPlot.M,rep(1,nrow(datPlot.M)),1,vT.Mp,colnames(datPlot.M),NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Modified Variables, No Grouping",cex=2)
my.pairs2(datPlot.O,locIndep.O$class,locIndep.O$G,vT.Op,colnames(datPlot.O),NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Local Independence Clustering - Original Variables",cex=2)
my.pairs2(datPlot.M,locIndep.M$class,locIndep.M$G,vT.Mp,colnames(datPlot.M),NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Local Independence Clustering - Modified Variables",cex=2)
my.pairs2(datPlot.O,bstModel.O$class,bstModel.O$G,vT.Op,colnames(datPlot.O),NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Best Model - All Original Variables",cex=2)
my.pairs2(datPlot.M,bstModel.M$class,bstModel.M$G,vT.Mp,colnames(datPlot.M),NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"Best Model - All Modified Variables",cex=2)
my.pairs2(datPlot.O[,which(deanVS.O$SW$allRoles=="S")],deanVS.O$SW$fBIC$clust$class,deanVS.O$SW$fBIC$clust$G,
          vT.Op,colnames(datPlot.O)[which(deanVS.O$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"VS - Restricted Models (no J) - All Original Variables",cex=1)
my.pairs2(datPlot.M[,which(deanVS.M$SW$allRoles=="S")],deanVS.M$SW$fBIC$clust$class,deanVS.M$SW$fBIC$clust$G,
          vT.Mp,colnames(datPlot.M)[which(deanVS.M$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"VS - Restricted Models (no J) - All Modified Variables",cex=1)
my.pairs2(datPlot.O[,which(deanVS2.O$SW$allRoles=="S")],deanVS2.O$SW$fBIC$clust$class,deanVS2.O$SW$fBIC$clust$G,
          vT.Op,colnames(datPlot.O)[which(deanVS2.O$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"VS - All Models - All Original Variables",cex=2)
my.pairs2(datPlot.M[,which(deanVS2.M$SW$allRoles=="S")],deanVS2.M$SW$fBIC$clust$class,deanVS2.M$SW$fBIC$clust$G,
          vT.Mp,colnames(datPlot.M)[which(deanVS2.M$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=1.75,outer=TRUE,"VS - All Models - With Modified Variables",cex=2)


plot.bic(saveList[[i]]$SW$fBIC$clust,methodNames[i],r,c(2,4),"BIC",TRUE)

newClass <- NULL
for(i in 1:nrow(datPlot.M)){
 m <- max(deanVS2.M$SW$fBIC$clust$z[i,])
 o <- order(deanVS2.M$SW$fBIC$clust$z[i,],decreasing=TRUE)
 newClass[i] <- ifelse(m>.55,o[1],o[2])
}
my.pairs2(datPlot.M[,which(deanVS2.M$SW$allRoles=="S")],
          (2*newClass-1)+(stage2-3),deanVS2.M$SW$fBIC$clust$G*2,
          vT.Mp,colnames(datPlot.M)[which(deanVS2.M$SW$allRoles=="S")],NA,TRUE)
my.pairs2(datPlot.M[,which(deanVS2.M$SW$allRoles=="S")],
          stage2-2,2,
          vT.Mp,colnames(datPlot.M)[which(deanVS2.M$SW$allRoles=="S")],NA,TRUE)
mtext(side=3,line=2.5,outer=TRUE,"VS - All Models, Modified Variables",cex=1.5)
mtext(side=3,line=1.25,outer=TRUE,"3 groups separated by stage3, stage4",cex=.75)
dev.off()
##----------------------------------------
frx2 <- factor(rx2)
frx3 <- factor(rx3)
frx4 <- factor(rx4)
fstatus1 <- factor(status1)
fstatus2 <- factor(status2)
fstatus3 <- factor(status3)
fstatus4 <- factor(status4)

resTable <- rep(list(),8)
for(r in 1:8) resTable[[r]] <- rep(list(),8)
tabFunc <- function(trx,stat,statName,stage,resList){
  l <- length(unique(trx))
  tmp <- vector("list",length(resList))
 for(i in c(1:4,9)){
 for(g in 1:resList[[i]]$G){
 for(r in 3:4){
  tt <- table(trx[which(stage==r & resList[[i]]$class==g)],
              stat[which(stage==r & resList[[i]]$class==g)])
  tmp[[i]] <- rbind(tmp[[i]],cbind(rep(g,l),rep(r,l),c(1:l),tt/apply(tt,1,sum),apply(tt,1,sum)))
}}
colnames(tmp[[i]]) <- c("G","Stage","Tx",statName,"Row N")
rownames(tmp[[i]]) <- NULL
}
for(i in 5:8){
for(g in 1:resList[[i]]$SW$fBIC$clust$G){
for(r in 3:4){
  tt <-table(trx[which(stage==r & resList[[i]]$SW$fBIC$clust$class==g)],
             stat[which(stage==r & resList[[i]]$SW$fBIC$clust$class==g)])
  tmp[[i]] <- rbind(tmp[[i]],cbind(rep(g,l),rep(r,l),c(1:l),tt/apply(tt,1,sum),apply(tt,1,sum)))
}}
colnames(tmp[[i]]) <- c("G","Stage","Tx",statName,"Row N")
rownames(tmp[[i]]) <- NULL
}
  return(tmp)
}

t1a <- tabFunc(frx4,fstatus2,newOutcome2,stage2,saveList)
t1a2 <- tabFunc(frx4,fstatus2,newOutcome2,rep(3,nrow(datClst.M)),saveList)
t1b <- tabFunc(frx3,fstatus2,newOutcome2,stage2,saveList)
t1c <- tabFunc(frx2,fstatus2,newOutcome2,stage2,saveList)
t2a <- tabFunc(frx4,fstatus3,newOutcome3,stage2,saveList)
t2b <- tabFunc(frx3,fstatus3,newOutcome3,stage2,saveList)
t2c <- tabFunc(frx2,fstatus3,newOutcome3,stage2,saveList)
t3a <- tabFunc(frx4,fstatus4,newOutcome4,stage2,saveList)
t3b <- tabFunc(frx3,fstatus4,newOutcome4,stage2,saveList)
t3c <- tabFunc(frx2,fstatus4,newOutcome4,stage2,saveList)
useT<- tabFunc(frx4,fstatus1,newOutcome1,stage2,saveList)

useT2 <-list()
for(i in 1:length(t1a2)){ 
  useT2[[i]] <- t1a2[[i]][!apply(
  t1a2[[i]],1,function(x) sum(is.na(x)))>0,-2]
}

table(saveList[[6]]$SW$fBIC$clust$class[age>65],stage2[age>65])
table(saveList[[2]]$class[age>65],stage2[age>65])
table(comp.pp$rx4[saveList[[6]]$SW$fBIC$clust$class==3],
      comp.pp$status2[saveList[[6]]$SW$fBIC$clust$class==3],
      comp.pp$stage[saveList[[6]]$SW$fBIC$clust$class==3])
                    
sink("prostateTables.tex")
cat("\\section{Suggested tables: 4 trx, 3 out}")
for(i in c(1,3,5,2,6,8)){
  dx <- xtable(useT2[[i]],caption=paste(methodNames[i+2],": Compare to stage, treatment (4), outcome classes"))
  digits(dx) <- c(rep(0,3),rep(2,4),0)
  print(dx,include.colnames=TRUE,include.rownames=FALSE,caption.placement="top")
}

cat("\\section{Suggested tables: 4 trx, 3 out}")
for(i in c(1,3,5,2,6)){
  dx <- xtable(useT[[i]],caption=paste(methodNames[i+2],": Compare to stage, treatment (4), outcome classes"))
  digits(dx) <- c(rep(0,4),rep(2,3),0)
  print(dx,include.colnames=TRUE,include.rownames=FALSE,caption.placement="top")
}

cat("\\clearpage \\\\")
cat("\\section{Original versions of tables: 2 trx, 4 out}")
for(i in c(1,3,5,2,6)){
  dx <- xtable(t1b[[i]],caption=paste(methodNames[i+2],": Compare to stage, treatment (4), outcome classes"))
  digits(dx) <- c(rep(0,4),rep(2,4),0)
  print(dx,include.colnames=TRUE,include.rownames=FALSE,caption.placement="top")
}
sink()
cat("\\clearpage \\\\")
cat("\\section{Original versions of tables: 4 trx, 4 out }")
for(i in c(1,3,5,2,6)){
  dx <- xtable(resTable2[[i]],caption=paste(methodNames[i+2],": Compare to stage, treatment (2A), outcome classes"))
  digits(dx) <- c(rep(0,4),rep(2,4),0)
  print(dx,include.colnames=TRUE,include.rownames=TRUE,caption.placement="top")
}
cat("\\clearpage \\\\")
for(i in 1:8){
  dx <- xtable(resTable3[[i]],caption=paste(methodNames[i+2],": Compare to stage, treatment (my2), outcome classes"))
  digits(dx) <- c(rep(0,4),rep(2,4),0)
  print(dx,include.colnames=TRUE,include.rownames=TRUE,caption.placement="top")
}
sink()

##C1,S3

##---------------------------------------------
pdf("prostateBIC.pdf")
x11(width=6.5,height=5)
for(i in c(1,3)){
par(mfrow=c(1,2))
r <- range(saveList[[i]]$BIC,saveList[[i+1]]$BIC,na.rm=TRUE)

plot.bic(saveList[[i]]$SW$fBIC$clust,methodNames[i],r,c(2,4),"BIC",TRUE)
plot.bic2(saveList[[i]],methodNames[i],r,"",TRUE,TRUE,FALSE)
plot.bic2(saveList[[i+1]],methodNames[i+1],r,"BIC",TRUE,FALSE,TRUE)
plot.bic2(saveList[[i+1]],methodNames[i+1],r,"",TRUE,TRUE,FALSE)
}
for(i in c(5,7)){
par(mfrow=c(1,2))
r <- range(saveList[[i]]$SW$fBIC$clust$BIC,saveList[[i+1]]$SW$fBIC$clust$BIC,na.rm=TRUE)
plot.bic2(saveList[[i]]$SW$fBIC$clust,methodNames[i],r,"BIC",TRUE,FALSE,TRUE)
plot.bic2(saveList[[i]]$SW$fBIC$clust,methodNames[i],r,"",TRUE,TRUE,FALSE)
plot.bic2(saveList[[i+1]]$SW$fBIC$clust,methodNames[i+1],r,"BIC",TRUE,FALSE,TRUE)
plot.bic2(saveList[[i+1]]$SW$fBIC$clust,methodNames[i+1],r,"",TRUE,TRUE,FALSE)
}
dev.off()
