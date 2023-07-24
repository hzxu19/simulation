##Checks - A is DxD with entries only on diagonal
##         L is Gx1 for VEI, VVI
##         z is NxG, data is NxD

##----------------------------------------
##             Hyperpriors
##----------------------------------------
tol = 1*10^(-5)
.mu.hyp <- function(.data) sum(.data,na.rm=TRUE)/length(!is.na(.data)) #mean(colMeans(.data))
mu.hyp <- cmpfun(.mu.hyp)
##---------------------
kp.hyp <- .01
##---------------------
.nu.hyp <- function(.data,.D) .D+2
nu.hyp <- cmpfun(.nu.hyp)
##---------------------
.psi.hyp <- function(.data,.D,.G)    {sum(diag(var(.data)))/.D}/{.G^(2/.D)}
psi.hyp <- cmpfun(.psi.hyp)
##---------------------
.lam.hyp <- function(.data,.D,.G)  tmp <- var(.data)/{.G^(2/.D)}
lam.hyp <- cmpfun(.lam.hyp)
##---------------------
.Hyper <- function(.data,.D,.G,.kp,.tol,.max.iter){
  mu.h <- mu.hyp(.data)
  kp.h <- kp.hyp
  nu.h <- nu.hyp(.data,.D)
  psi.h <- psi.hyp(.data,.D,.G)
  lam.h <- lam.hyp(.data,.D,.G)
  return(list(mu.h=mu.h,kp.h=kp.h,nu.h=nu.h,psi.h=psi.h,lam.h=lam.h,kp=.kp,tol=.tol,max.iter=max.iter))
}
Hyper <- cmpfun(.Hyper)
#mu.hyp(.data); nu.hyp(.data,4); psi.hyp(.data,4,2); lam.hyp(.data,4,2)
##----------------------------------------
##         Functions for estimators
##                NO PRIOR
##----------------------------------------
.n.g <- function(.z,.g)  sum(.z[,.g])
n.g <- cmpfun(.n.g)
##---------------------
.ybar.g <- function(.z,.data,.D,.N,.g,.n){
  jnk = .z[,.g]*.data
  tmp <- apply(jnk,2,sum)/.n
  return(tmp)
}
ybar.g <- cmpfun(.ybar.g)
##---------------------
.w.g <- function(.z,.data,.D,.N,.g,.n){
  #p3 <- vector("list",.N)
  p3 <- matrix(0,.D,.D)
  y.g <- ybar.g(.z,.data,.D,.N,.g,.n)
  .tmpdata = as.matrix(.data)
  for(j in 1:.N){
    p1 <- .tmpdata[j,]-y.g #as.numeric(.data[j,] - y.g)
    p2 <- p1 %*% t(p1)
    #p3[[j]] <- .z[j,.g] * p2
    p3 <- p3+ .z[j,.g]*p2
  }
  #tmp <- Reduce("+",p3)
  return(p3)
  #return(tmp)
}
w.g <- cmpfun(.w.g)
#n.g(z.init,1);ybar.g(z.init,.data,1);w.g(z.init,.data,1) 
##---------------------
.w.prior <- function(.z,.data,.D,.N,.g,.n,.hyp){
  y.g <- ybar.g(.z,.data,.D,.N,.g,.n)
  n.21 <- {.hyp$kp.h*.n} / {.hyp$kp.h +.n}
  jnk <- y.g - .hyp$mu.h
  n.22 <- jnk %*% t(jnk)
  n.23 <- w.g(.z,.data,.D,.N,.g,.n)
  tmp <- n.21*n.22 + n.23
  return(tmp)
}
w.prior <- cmpfun(.w.prior)
#w.prior(z.init,.data,1,.kp)
##---------------------
.mu.G <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  jnk <- tmp <- matrix(NA,.D,.G) 
  for(g in 1:.G) jnk[,g] <- ybar.g(.z,.data,.D,.N,g,.ng[g])

  if(length(.which.pri)==0){
    tmp <- jnk
  }else{
    if(.which.pri=="DEF" | .which.pri=="MOD"){
      for(g in 1:.G){
        num1 <- .ng[g]*jnk[,g]
        num1[is.na(num1)] <- 0
        tmp[,g] <- {num1+.hyp$kp.h*.hyp$mu.h} / {.ng[g] + .hyp$kp.h}
  }}}
  return(tmp)
}
mu.G <- cmpfun(.mu.G)
#mu.G(z.init,.data,2,150,2,ng,.which.pri,.hyp)
##----------------------------------------
##         SPHERICAL ESTIMATORS
##----------------------------------------
##--------         EII        ------------
##----------------------------------------
.L.EII <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
   w <- 0

  if(is.null(.which.pri)){
     for(g in 1:.G) w <- w +tr(w.g(.z,.data,.D.,.N,.g,.ng[g]))
     #W <- Reduce("+",w)
     dnm <- .N*.D
   }else{  
     for(g in 1:.G) w <- w+ tr(w.prior(.z,.data,.D,.N,g,.ng[g],.hyp))#w[[g]] <- tr(w.prior(.z,.data,.D,.N,g,.ng[g],.hyp))
     w <- w+.hyp$psi.h# W <- Reduce("+",w) + .hyp$psi.h
     dnm <- (.N+.G)*.D + .hyp$nu.h + 2
   }
  tmp <- w/dnm #as.numeric(W)/dnm
  return(tmp)
}
L.EII <-cmpfun(.L.EII)
##----------------------------------------
##--------         VII        ------------
##----------------------------------------
.L.VII <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  #w <- vector("list",.G) 
  w <- dnm <- rep(NA,.G)

  if(is.null(.which.pri)){
     for(g in 1:.G){
	    w[g]  <- tr(w.g(.z,.data,.D.,.N,.g,.ng[g])) #w[[g]] <- tr(w.g(.z,.data,.D.,.N,.g,.ng[g]))
	    dnm[g] <- .ng[g]*.D
 	   }
     	#W <- unlist(w)   
   }else{  
     for(g in 1:.G){
	    w[g] <- tr(w.prior(.z,.data,.D,.N,g,.ng[g],.hyp)) #w[[g]] <- tr(w.prior(.z,.data,.D,.N,g,.ng[g],.hyp))
	    dnm[g] <- (.ng[g]+1)*D + .hyp$nu.h + 2
     }
    w + .hyp$psi.h# W <- unlist(w) + .hyp$psi.h
   }
  tmp <- w/dnm #as.numeric(W)/dnm; 
  return(tmp)
}
L.vII <-cmpfun(.L.VII)
##----------------------------------------
##         DIAGONAL ESTIMATORS
##----------------------------------------
##--------         EEI        ------------
##----------------------------------------
.L.EEI <- function(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp){
  if(is.null(.which.pri)){
     W <- Reduce("+",.w)
     dnm <- .N
   }else{  
     if(.which.pri=="DEF")       pri <- .hyp$psi.h*diag(.D)
     if(.which.pri=="MOD")       pri <- diag(.hyp$lam.h)*diag(.D)      
     W <- Reduce("+",.w) + pri
     dnm <- .N + .hyp$nu.h + .G + 2
   }
  num <- det(diag(W)*diag(.D))^{1/.D} 
  tmp <- num/dnm
  return(tmp)
}
##---------------------
.A.EEI <- function(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp){
  if(is.null(.which.pri)){
    W <- Reduce("+",.w)
  }else{
    if(.which.pri=="DEF")       pri <- .hyp$psi.h*diag(.D)
    if(.which.pri=="MOD")       pri <- diag(.hyp$lam.h)*diag(.D) 
    W <- Reduce("+",.w) + pri
  } 
  tmp <- diag(W)/det(diag(W)*diag(.D))^{1/.D}
  return(tmp)
}
L.EEI <-cmpfun(.L.EEI)
A.EEI <-cmpfun(.A.EEI)

##----------------------------------------
##--------         EVI        ------------
##----------------------------------------
.L.EVI <-  function(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp){  
  dtr <- rep(NA,.G)  
 if(is.null(.which.pri)){  dnm <- .N
   }else{                  dnm <- .N + .hyp$nu.h + 2    }

  for(g in 1:.G) dtr[g] <- det(diag(.w[[g]])*diag(.D))^{1/.D}
  tmp <- sum(dtr) / dnm
  return(tmp)
}
##-----------
.A.EVI <-  function(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp){  
  tmp <- matrix(NA,.D,.G)
  for(g in 1:.G)    tmp[,g] <- diag(.w[[g]]) / det(diag(.w[[g]])*diag(.D))^{1/.D}
  return(tmp)
}
L.EVI <-cmpfun(.L.EVI)
A.EVI <-cmpfun(.A.EVI)
##----------------------------------------
##--------         VVI        ------------
##----------------------------------------
.L.VVI <- function(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp){  
  dnm <- tmp <- rep(NA,.G)
  
 if(is.null(.which.pri)){      dnm <- 0
   }else{                      dnm <- .hyp$nu.h + 3   }
  for(g in 1:.G) tmp[g] <- det(diag(.w[[g]])*diag(.D))^{1/.D} / {.ng[g]+dnm}
  return(tmp)
}
##-----------
.A.VVI <-  function(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp){
  tmp <- A.EVI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)
  return(tmp)
}
L.VVI <-cmpfun(.L.VVI)
A.VVI <-cmpfun(.A.VVI)
##----------------------------------------
##--------         VEI        ------------
##----------------------------------------
.L.VEI <- function(.z,.data,.D,.N,.G,.ng,.w,.A,.which.pri,.hyp){
  ##A should be a DxD matrix
  # w <- vector("list",.G)
  tmp <- rep(NA,.G)
  
  if(is.null(.which.pri)){ dnm = 0
  }else{  
    if(.which.pri=="DEF")       pri <- .hyp$psi.h*diag(.D)
    if(.which.pri=="MOD")       pri <- diag(.hyp$lam.h)*diag(.D) 
    dnm <-.hyp$nu.h + .G + 2
  }
  for(g in 1:.G){
    tmp.num <- sum(diag(.w[[g]] %*% solve(.A*diag(.D))))
    tmp[g] <- tmp.num / {.D*.ng[g]+dnm}
  }
  return(tmp)
}
##---------------------
.A.VEI <- function(.z,.data,.D,.N,.G,.ng,.w,.l,.which.pri,.hyp){
  ##l should be a vector of length G
  num <- diag(Reduce("+",.w))
  tmp <- num/det(num*diag(.D))^{1/.D}
  return(tmp)
}
##---------------------
.tst.VEI <- function(.z,.data,.D,.N,.G,.ng,.w,.l,.A,.which.pri,.hyp){
  ##A should be a Dx1 vector
  ##l should be a vector of length G
  p1 <- p2 <- rep(NA,.G)
  for(g in 1:.G){
    jnk <- .w[[g]] %*% solve(.A*diag(.D))
    p1[g] <- tr(jnk) #sum(diag(jnk))
    p2[g] <- .ng[g]*log(.l[g])
  }
  tmp <- sum(p1) + .D*sum(p2)
  return(tmp)
}
##-----------
.iter.VEI <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  tmp.L <- tmp.A <- min.F <- rep(list(),.hyp$max.iter)
  iter <- 1
  w <- vector("list",.G)
  
  if(is.null(.which.pri)){
    for(g in 1:.G) w[[g]] <- w.g(.z,.data,.D.,.N,g,.ng[g])
  }else{  
    if(.which.pri=="DEF" | .which.pri=="MOD"){
      for(g in 1:.G) w[[g]] <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp) #+ pri
    }}
  ##PROPER WAY TO INITIALIZE????
  tmp.L[[iter]] <- L.VVI(.z,.data,.D,.N,.G,.ng,w,.which.pri,.hyp)
  tmp.A[[iter]] <- A.VEI(.z,.data,.D,.N,.G,.ng,w,tmp.L[[iter]],.which.pri,.hyp)
  min.F[[iter]] <- tst.VEI(.z,.data,.D,.N,.G,.ng,w,tmp.L[[iter]],tmp.A[[iter]],.which.pri,.hyp)
  tst <- 1
  
  while(tst > .hyp$tol & iter < .hyp$max.iter){
    iter <- iter+1
    tmp.L[[iter]] <- L.VEI(.z,.data,.D,.N,.G,.ng,w,tmp.A[[iter-1]],.which.pri,.hyp)
    tmp.A[[iter]] <- A.VEI(.z,.data,.D,.N,.G,.ng,w,tmp.L[[iter]],.which.pri,.hyp)
    min.F[[iter]] <- tst.VEI(.z,.data,.D,.N,.G,.ng,w,tmp.L[[iter]],tmp.A[[iter]],.which.pri,.hyp)
    tst <- abs(min.F[[iter]]-min.F[[iter-1]]) ##Is F monotone? If so, abs() shouldn't be needed?
  }
  
  scale <- tmp.L[[iter]]
  shape <- tmp.A[[iter]]
  return(list(scale=scale,shape=shape))
}
L.VEI <-cmpfun(.L.VEI)
A.VEI <-cmpfun(.A.VEI)
tst.VEI <-cmpfun(.tst.VEI)
iter.VEI <-cmpfun(.iter.VEI)
##----------------------------------------
##         ELLIPSOIDAL ESTIMATORS
##----------------------------------------
##--------         EEE        ------------
##----------------------------------------
S.EEE <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  w <- vector("list",.G)

  if(is.null(.which.pri)){
     for(g in 1:.G) w[[g]] <- w.g(.z,.data,.D.,.N,g,.ng[g])
     W <- Reduce("+",w)
     dnm <- .N
   }else{  
     pri <- .hyp$lam.h     
     for(g in 1:.G) w[[g]] <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp)
     W <- Reduce("+",w) + pri
     dnm <- .N+.D+.G + .hyp$nu.h + 1
   }
  tmp <- W/dnm
  return(tmp)
}
#S.EEE(z,.data,g,which.pri,kp)
##----------------------------------------
##--------         EEV        ------------
##----------------------------------------
.L.EEV <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  omega <- vector("list",.G)
  dnm <- tmp <- rep(NA,.G)

  if(is.null(.which.pri)){
     for(g in 1:.G){
      .w <- w.g(.z,.data,.D.,.N,g,.ng[g])
	    omega[[g]] <- eigen(.w)$values*diag(.D)
      }
     dnm = .N
   }else{  
    pri <- .hyp$lam.h*diag(.D)   
    for(g in 1:.G){
      .w <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp) + pri
	    omega[[g]] <- eigen(.w)$values*diag(.D)
    }
     dnm <- .N + .hyp$nu.h + .G + 2
   }
   o = Reduce("+",omega) 
   tmp <- det(o)^{1/.D} / dnm
  return(tmp)
}
##---------------------
.AC.EEV <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  omega <- vector("list",.G)
  C <- array(NA,dim=c(.D,.D,.G))

 if(is.null(.which.pri)){
     for(g in 1:.G){
      .W <- w.g(.z,.data,.D.,.N,g,.ng[g])
 	    .eig <- eigen(.w[[g]])
	    omega[[g]] <- .eig$values*diag(.D)
	    C[,,g] <- .eig$vectors
     }
   }else{  
    pri <- .hyp$lam.h*diag(.D)   
    for(g in 1:.G){
      .w <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp) + pri
    	.eig <- eigen(.w)
    	omega[[g]] <- .eig$values*diag(.D)
    	C[,,g] <- .eig$vectors
   }}
   num <- Reduce("+",omega)
    A <- num / det(num)^{1/.D}
  return(list(A=A,C=C))
}
L.EEV <-cmpfun(.L.EEV)
AC.EEV <-cmpfun(.AC.EEV)
##----------------------------------------
##--------         VEV (IP)   ------------
##----------------------------------------
.L.VEV <- function(.z,.data,.D,.N,.G,.ng,.w,.A,.C,.which.pri,.hyp){
  ##A should be  DxD matrix 
  ##C a G length list of DxD matrix
  #w <- vector("list",.G)
  dnm <- tmp <- rep(NA,.G)
  
  if(is.null(.which.pri)){ 
    pri <- 0
    dnm <- 0
   }else{  
     pri <- .hyp$lam.h*diag(.D)     
     dnm <- .hyp$nu.h + .G + 2
   }
  for(g in 1:.G)   tmp[g] <- tr({.w[[g]]+pri}%*%.C[,,g] %*%solve(.A)%*% t(.C[,,g]))/ {.ng[g]*.D+dnm}
 
  return(tmp)
}
##---------------------
.AC.VEV <- function(.z,.data,.D,.N,.G,.w,.l,.which.pri,.hyp){
  omega <- tmpA <- vector("list",.G)
  C <- array(NA,dim=c(.D,.D,.G))
 
 if(is.null(.which.pri)){
     for(g in 1:.G){
       .eig <- eigen(.w)
	    C[,,g] <- .eig$vectors
	    omega[[g]] <- .eig$values*diag(.D)
     }
   }else{  
    pri <- .hyp$lam.h*diag(.D)   
    for(g in 1:.G){
      .eig <- eigen(.w)
	    C[,,g] <- .eig$vectors
	    omega[[g]] <- .eig$values*diag(.D)
   }}
  for(g in 1:.G) tmpA[[g]] <- omega[[g]]/.l[g]
  numA = Reduce("+",tmpA)
  A = numA/det(numA)^{1/.D}
return(list(C=C,A=A))
}
##---------------------
.tst.VEV <- function(.z,.data,.D,.N,.G,.w,.l,.A,.C,.which.pri,.hyp){
  ##A should be a DxD matrix
  ##C a G-list of Dxd matrix
  ##l should be a vector of length G
  p1 <- p2 <- rep(NA,.G)
 for(g in 1:.G){
   jnk <- .w[[g]]%*%.C[,,g]%*%solve(.A)%*%t(.C[,,g])
   p1[g] <- tr(jnk)/.l[g]
   p2[g] <- n.g(.z,g)*log(.l[g])
  }
  tmp <- sum(p1) + D*sum(p2)
  return(tmp)
}
##-----------
.iter.VEV <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  tmp.L <- tmp.A <- min.F <- rep(list(),.hyp$max.iter)
  iter <- 1
  .w <- vector("list",.G)
   if(is.null(.which.pri)){
       for(g in 1:.G) .w[[g]] <- w.g(.z,.data,.D.,.N,g,.ng[g])
    }else{  
       pri <- .hyp$lam.h*diag(.D)
      for(g in 1:.G) .w[[g]] <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp) + pri
    }
  ##PROPER WAY TO INITIALIZE????
  tmp.L[[iter]] <- L.VVI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)
  tmp.A[[iter]] <- AC.VEV(.z,.data,.D,.N,.G,.w,tmp.L[[iter]],.which.pri,.hyp)
  min.F[[iter]] <- tst.VEV(.z,.data,.D,.N,.G,.w,tmp.L[[iter]],tmp.A[[iter]]$A,tmp.A[[iter]]$C,.which.pri,.hyp)
  tst <- 1

  while(tst > .hyp$tol & iter < .hyp$max.iter){
   iter <- iter+1
   tmp.L[[iter]] <- L.VEV(.z,.data,.D,.N,.G,.w,tmp.A[[iter-1]]$A,tmp.A[[iter-1]]$C,.which.pri,.hyp)
   tmp.A[[iter]] <- AC.VEV(.z,.data,.D,.N,.G,.w,tmp.L[[iter]],.which.pri,.hyp)
   min.F[[iter]] <- tst.VEV(.z,.data,.D,.N,.G,.w,tmp.L[[iter]],tmp.A[[iter]]$A,tmp.A[[iter]]$C,.which.pri,.hyp)
   tst <- abs(min.F[[iter]]-min.F[[iter-1]]) ##Is F monotone? If so, abs() shouldn't be needed?
  }

 scale  <- tmp.L[[iter]]
 shape  <- tmp.A[[iter]]$A
 orient <- tmp.A[[iter]]$C
return(list(scale=scale,shape=shape,orientation=orient))
}
L.VEV <- cmpfun(.L.VEV)
AC.VEV <- cmpfun(.AC.VEV)
tst.VEV <- cmpfun(.tst.VEV)
iter.VEV <- cmpfun(.iter.VEV)
##----------------------------------------
##--------         VVV        ------------
##----------------------------------------
.S.VVV <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  w <- vector("list",.G)
  tmp <- array(NA,dim=c(.D,.D,.G))
  dnm <- rep(NA,.G)
 
 if(is.null(.which.pri)){
     for(g in 1:.G)       w[[g]] <- w.g(.z,.data,.D.,.N,g,.ng[g])
     dnm[g] = .ng[g]
   }else{  
    pri <- .hyp$lam.h*diag(.D)   
    for(g in 1:.G){
     w[[g]] <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp) + pri
     dnm[g] <- N + .hyp$nu.h + .G + 2
   }}
   for(g in 1:.G) tmp[,,g] <- w[[g]] / dnm[g]
  return(tmp)
}
S.VVV <- cmpfun(.S.VVV)
##----------------------------------------
##         EM - DIAGONAL MODELS
##----------------------------------------
##--------------
.m.EEI <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  .w <- vector("list",.G)    
  if(is.null(.which.pri)){
    for(g in 1:.G) .w[[g]] <- w.g(.z,.data,.D.,.N,g,.ng[g])
  }else{  
    if(.which.pri=="DEF")       pri <- .hyp$psi.h*diag(.D)
    if(.which.pri=="MOD")       pri <- diag(.hyp$lam.h)*diag(.D) 
    for(g in 1:.G) .w[[g]] <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp)
  }
  tmp.u <- mu.G(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp)
  tmp.L <- L.EEI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)
  tmp.A <- A.EEI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)

  pro <- apply(.z,2,sum)/.N
  mean <- tmp.u
  variance <- list()
  variance$d <- .D
  variance$G <- .G
  variance$scale <- tmp.L
  variance$shape <- tmp.A
  variance$Sigma <- tmp.L*tmp.A*diag(.D)
  variance$sigma <- array(NA,dim=c(.D,.D,.G))
  for(g in 1:.G) variance$sigma[,,g] <- variance$Sigma
  return(list(pro=pro,mean=mean,variance=variance))
}
#m.EEI(z.init,.data,2,NULL,.kp)
m.EEI <- cmpfun(.m.EEI)
##--
.m.VEI <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  tmp.u <- mu.G(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp)
  tmp <- iter.VEI(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp)

  pro <- apply(.z,2,sum)/.N
  mean <- tmp.u
  variance <- list()
  variance$d <- .D
  variance$G <- .G
  variance$scale <- tmp$scale
  variance$shape <- tmp$shape
  variance$sigma <- array(NA,dim=c(.D,.D,.G))
  for(g in 1:.G) variance$sigma[,,g] <- tmp$scale[g]*tmp$shape*diag(.D)
  return(list(pro=pro,mean=mean,variance=variance))
}
m.VEI <- cmpfun(.m.VEI)
#m.VEI(z.init,.data,2,NULL,.kp)
##--
.m.EVI <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  .w <- vector("list",.G)
  if(is.null(.which.pri)){
    for(g in 1:.G)       .w[[g]] <- w.g(.z,.data,.D.,.N,g,.ng[g])
  }else{  
    if(.which.pri=="DEF")       pri <- .hyp$psi.h*diag(.D)
    if(.which.pri=="MOD")       pri <- diag(.hyp$lam.h)*diag(.D) 
    for(g in 1:.G)       .w[[g]] <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp) + pri
  }
  
  tmp.u <- mu.G(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp)
  tmp.L <- L.EVI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)
  tmp.A <- A.EVI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)
  
  pro <- apply(.z,2,sum)/.N
  mean <- tmp.u
  variance <- list()
  variance$d <- .D
  variance$G <- .G
  variance$scale <- tmp.L
  variance$shape <- tmp.A
  variance$sigma <- array(NA,dim=c(.D,.D,.G))
  for(g in 1:.G) variance$sigma[,,g] <- tmp.L*tmp.A[,g]*diag(.D)
  return(list(pro=pro,mean=mean,variance=variance))
}
m.EVI <- cmpfun(.m.EVI)
#m.EVI(z.init,.data,2,NULL,.kp)
##--
.m.VVI <- function(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  .w <- vector("list",.G)
  if(is.null(.which.pri)){
    for(g in 1:.G)        .w[[g]] <- w.g(.z,.data,.D.,.N,g,.ng[g])
    dnm <- 0
  }else{  
    if(.which.pri=="DEF")       pri <- .hyp$psi.h*diag(.D)
    if(.which.pri=="MOD")       pri <- diag(.hyp$lam.h)*diag(.D)  
    for(g in 1:.G)    .w[[g]] <- w.prior(.z,.data,.D,.N,g,.ng[g],.hyp) + pri      
  }
  tmp.u <- mu.G(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp)
  tmp.L <- L.VVI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)
  tmp.A <- A.VVI(.z,.data,.D,.N,.G,.ng,.w,.which.pri,.hyp)
  
  pro <- apply(.z,2,sum)/.N
  mean <- tmp.u
  variance <- list()
  variance$d <- .D
  variance$G <- .G
  variance$scale <- tmp.L
  variance$shape <- tmp.A
  variance$sigma <- array(NA,dim=c(.D,.D,.G))
  for(g in 1:.G) variance$sigma[,,g] <- tmp.L[g]*tmp.A[,g]*diag(.D)
  return(list(pro=pro,mean=mean,variance=variance))
}
m.VVI <- cmpfun(.m.VVI)
#m.EVI(z.init,.data,2,NULL,.kp)
##--------
## For non-diagonal models
.gen.mstep <- function(.modelName,.z,.data,.D,.N,.G,.ng,.which.pri,.hyp){
  pro <- apply(.z,2,sum)/.N
  tmp.u <- mu.G(.z,.data,.D,.N,.G,.ng,.which.pri,.hyp)
  #tmp.var = list()                                    #  I add
  #variance = list()                                   #  I add
  if(!is.null(.which.pri)){ which.pri = priorControl()
  }else{ which.pri = .which.pri}
  #begin to change
  # original: if(length(.modelName)==1){tmp.var <- mstep(modelName = .modelName,data = .data,z = .z,prior=which.pri)     
  # standard: mstep(modelName = .modelName,data = .data,z = .z,prior=priorControl(functionName = which.pri) )
  #if(length(.modelName)==1){tmp.var <- my.m.step(modelName = .modelName,data = .data,z = .z ,prior=which.pri)   ORIGINAL
  if(length(.modelName)==1){tmp.var <- my.m.step(.modelName,.z,.data,.D,.N,.G,.which.pri=which.pri,.hyp)   # I delete ,prior=which.pri   change mstep to my.m.step
  }else{ jishu = length(.modelName)   
    for(i in 1:jishu)
      #   {tmp.var[[i]] <- my.m.step(modelName = .modelName[i],data = .data,.z,prior=which.pri )       ORIGINAL
    {tmp.var[[i]] <- my.m.step(.modelName = .modelName[i],.z,.data,.D,.N,.G,.which.pri=which.pri,.hyp)   # I delete ,prior=which.pri    change mstep to my.m.step
    variance[[i]] <- tmp.var[[i]]$parameters$variance    
    }
  }                                                                                                    # end   I add
# tmp.var <- mstep(modelName = .modelName,data = .data,.z,prior=which.pri)          #!!!!   modelName =
  
  return(list(pro=pro,G=.G,mean=tmp.u,variance=tmp.var$parameters$variance))
  #return(list(pro=pro,G=.G,mean=tmp.u,variance=variance))                 #  I add
}
gen.mstep <- cmpfun(.gen.mstep)
##--------
.my.m.step <- function(.modelName,.z,.data,.D,.N,.G,.which.pri,.hyp){
  ng <- rep(NA,.G)
  for(g in 1:.G) ng[g] <- n.g(.z,g) 
  if(all(.modelName %in% .diagModel) & !is.null(.which.pri)){         # I add all
   assign("tmp.m.func",eval(parse(text=paste("m.",.modelName,sep=""))))
   m.par <- tmp.m.func(.z,.data,.D,.N,.G,ng,.which.pri,.hyp)
  }else{
   m.par <- gen.mstep(.modelName,.z,.data,.D,.N,.G,ng,.which.pri,.hyp)        # ?
  }
  return(list(par=m.par))
}
my.m.step <- cmpfun(.my.m.step)
##--------
e.func <- function(.modelName,.data,.parameters){
 assign("tmp",eval(parse(text=paste("estep",.modelName,"(data = .data, parameters = .parameters)",sep=""))))      #?
 eFlag <- ifelse(is.na(tmp$loglik),"ERROR - Singluar",
                 ifelse(length(unique(map(tmp$z)))!=.parameters$variance$G,
                  "ERROR - Number of classes NEQ G",NA))
 return(list(model=.modelName,z=tmp$z,loglik=tmp$loglik,parameters=.parameters,	eFlag=eFlag))
}
##-- 
.nor.me <- function(.modelName,.data,.z.init,.g,.which.pri,.kp,.tol,.max.iter){
  D <- dim(.data)[2]
  N <- dim(.data)[1]
  
  if(!is.null(.which.pri)){  .hyp = Hyper(.data,D,.g,.kp,.tol,.max.iter)
   }else{                    .hyp = list(kp=.kp,tol=.tol,max.iter=.max.iter) }
  
  m.func <- my.m.step
  m <- e <- vector("list",.max.iter)
 
  m.init <- m.func(.modelName,.z.init,.data,D,N,.g,.which.pri,.hyp)        #?
  e.init <- e.func(.modelName,.data,m.init$par)          #?
  iter <- 1
 
  if(is.na(e.init$eFlag)){
   if(.g==1 & length(e.init$z) ==0) e.init$z <- matrix(1,N,1)
   m[[iter]] <- m.func(.modelName,e.init$z,.data,D,N,.g,.which.pri,.hyp)
   e[[iter]] <- e.func(.modelName,.data,m[[iter]]$par)
    tst <- e[[iter]]$loglik - e.init$loglik
  }else{ e[[1]] <- list(eFlag="ERROR - Singular") }
  
  tst = 1
  while(tst > .hyp$tol & iter < .hyp$max.iter & is.na(e[[iter]]$eFlag)){
   iter <- iter+1
   if(.g==1 & length(e[[iter-1]]$z) ==0) e[[iter-1]]$z <- matrix(1,N,1)
   m[[iter]] <- m.func(.modelName,e[[iter-1]]$z,.data,D,N,.g,.which.pri,.hyp)
   e[[iter]] <- e.func(.modelName,.data,m[[iter]]$par)
   tst <- e[[iter]]$loglik - e[[iter-1]]$loglik
  }

  res <- e[[iter]]
  z <- res$z
  if(!is.null(z)){ 
    class <- map(z)
   if(length(table(class)) != .g){
    loglik <- -Inf
    BIC <- NA
   }else{
    loglik <- res$loglik
    BIC <- bic(.modelName,loglik,N,D,.g)
   }}else{
    loglik <- NA
    BIC <- NA
  }
  if(is.na(e[[iter]]$eFlag)){
	return(list(z=z,class=class,loglik=loglik,bic=BIC,parameters=res$parameters,iter=iter))
  }else{
	return(list(z=z,class=rep(NA,N),loglik=loglik,bic=BIC,parameters=res$parameters,iter=iter))
  }
}
nor.me <- cmpfun(.nor.me)
##---------------------
##--------------------
## Consider more than 1 variance structure
.Mclustpri <- function(.data,.g1,.g2,.modelNames,.which.pri,.kp=.01,.tol=10^(-5),.max.iter=50){
  .D <- dim(.data)[2]
  if(is.null(.D)){.D = 0}          # I add
  if(.D==1){                           ##-- Beg: uni/multivariate
    NmodelName <- NULL			             
    if(any(.modelNames %in% .mclustModel[c(1,3,7)]))  NmodelName <- c(NmodelName,"E")
    if(any(.modelNames %in% .mclustModel[-c(1,3,7)])) NmodelName <- c(NmodelName,"V")
  }else{ NmodelName = .modelNames } 			##-- End: uni/multivariate
   mM = length(NmodelName)
   BIC <- matrix(NA,.g2-.g1+1,mM)
   dimnames(BIC) <- list(paste(.g1:.g2),NmodelName)
    
  if(.D==1){ tmphc <- hc(modelName="E",.data)
   }else{    tmphc <- hc(modelName="VII",.data)   }     # .data is NULL ?????
  tmp.cl <- hclass(tmphc,c(1:.g2))

  tmpMod = vector("list",length(NmodelName))
  names(tmpMod) <- NmodelName
  tmpG = vector("list",.g2)
  
  for(g in .g1:.g2){
    z.init <- unmap(tmp.cl[,g], groups=1:max(tmp.cl[,g]))
    tmpG[[g]] <- sapply(NmodelName, function(.mod) nor.me(.mod,.data,z.init,g,.which.pri,.kp,.tol,.max.iter),simplify=FALSE )    #  ?
    BIC[paste(g),] <- sapply(NmodelName, function(.mod) tmpG[[g]][[.mod]]$bic)
  }
  
  if(all(is.na(BIC))){
   cat(paste("No models could be fit for:",paste(.modelNames,collapse=",")," for",.g1,":",.g2,"\n",sep=""))
   mc=list(bic=NA)
  }else{
   BIC = matrix(as.numeric(BIC),dim(BIC))
   dimnames(BIC) <- list(paste(.g1:.g2),NmodelName)
   maxBIC = which(BIC==max(BIC,na.rm=TRUE),arr.ind=TRUE)
   maxG   = as.numeric(rownames(BIC)[maxBIC[1]])
   maxMod = colnames(BIC)[maxBIC[2]]
   mc     = tmpG[[maxG]][[maxMod]]
   mc$BIC = BIC
   mc$G   = maxG
   mc$model = maxMod
   mc$npar = {2 * mc$loglik - mc$bic}/log(dim(.data)[1])
   cat(paste("best norClust model:",maxMod,"with",maxG,"clusters \n"))
  }
  return(mc)
}
Mclustpri <- cmpfun(.Mclustpri)
##--------------------
#system.time(mc2 <-  Mclustpri(iris[,1:4],2,2,.mclustModel,"MOD",kp,tol,max.iter) )
#g=2;m="EII";.data=iris[,1:4];.g1=2;.g2=9;.modelName=.mclustModel;kp=.kp;tol=.tol;max.iter=.max.iter
#.modelName=m
#foo = Mclustpri(d[[1]][,1:2],2,2,"VEI","MOD",.kp=.01,.tol=10^(-5),.max.iter=50)
#load("d4_150.rdata"); .data=d[[1]][,1:2];.g1=.g2=2;.modelNames="VEI";.which.pri="MOD"