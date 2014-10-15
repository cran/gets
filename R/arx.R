arx <-
function(y, mc=NULL, ar=NULL, ewma=NULL, mxreg=NULL,
  arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL, zero.adj=0.1,
  vc.adj=TRUE, vcov.type=c("ordinary", "white"),
  qstat.options=NULL, tol=1e-07, LAPACK=FALSE, verbose=TRUE)
{
##arguments:
vcov.type <- match.arg(vcov.type)
aux <- list() #info for getsm function
if(NCOL(y)!=1) stop("Regressand must be a vector")
if(is.null(qstat.options)){
  if(is.null(ar)){ar.lag <- 1}else{ar.lag <- max(ar)+1}
  if(is.null(arch)){arch.lag <- 1}else{arch.lag <- max(arch)+1}
  qstat.options <- c(ar.lag, arch.lag)
}

##zoo:
y <- as.zoo(y)
y <- na.trim(y)
aux$y <- y
zoo.index.y <- index(y)
y <- coredata(y)
y.n <- length(y)
t1 <- zoo.index.y[1]
t2 <- zoo.index.y[y.n]
aux$t1 <- t1
aux$t2 <- t2

if(!is.null(mxreg)){
  mxreg <- as.zoo(mxreg)
  mxreg <- window(mxreg, start=t1, end=t2)
  mxreg <- cbind(coredata(mxreg))
}

if(!is.null(vxreg)){
  vxreg <- as.zoo(vxreg)
  vxreg <- window(vxreg, start=t1, end=t2)
  vxreg <- cbind(coredata(vxreg))
}

##aux:
aux$mc <- mc; aux$ar <- ar; aux$ewma <- ewma; aux$mxreg <- mxreg;
aux$arch <- arch; aux$asym <- asym; aux$log.ewma <- log.ewma;
aux$vxreg <- vxreg; aux$zero.adj <- zero.adj;
aux$vc.adj <- vc.adj; aux$vcov.type <- vcov.type;
aux$qstat.options <- qstat.options; aux$tol <- tol;
aux$LAPACK <- LAPACK;

### INITIALISE ##########

out <- list()
out$aux <- aux
out$call <- sys.call()

##regressors:
mX <- regs.mean(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg)
ewma.mean.chk <- if(is.null(ewma)){0}else{ifelse(is.null(ewma$lag),1,ewma$lag)}
max.ar <- if( is.null(ar)&&is.null(ewma) ){0}else{ max(ar,ewma.mean.chk) }
vYadj <- y[I(max.ar+1):y.n]
if(!is.null(mX)){mXadj <- cbind(mX[I(max.ar+1):y.n,])}

#### MEAN #######################

if(is.null(mX)){
  resids <- vYadj
  mean.results <- NULL
  vcov.mean <- NULL
}else{
  est.m <- ols(vYadj, mXadj, tol = tol, LAPACK=LAPACK,
    method=2)
  fit.m <- as.vector(mXadj%*%cbind(est.m$coefficients))
  resids <- vYadj - fit.m
  resids2 <- resids^2
  mXadj.n <- NROW(mXadj)
  mXadj.k <- NCOL(mXadj)
  d.f. <- mXadj.n - mXadj.k
  sigma2 <- sum(resids2)/d.f.

  #estimate s.e.; compute t-stats. and p-vals.:
  if(vcov.type == "ordinary"){
    varcovmat <- sigma2*est.m$xtxinv
    s.e. <- sqrt(as.vector(diag(varcovmat)))
  }
  if(vcov.type == "white"){
    matResids2 <- matrix(0, mXadj.n, mXadj.n)
    diag(matResids2) <- resids2
    omega.hat <- t(mXadj)%*%matResids2%*%mXadj
    varcovmat <- est.m$xtxinv%*%omega.hat%*%est.m$xtxinv
    coef.var <- as.vector(diag(varcovmat))
    s.e. <- sqrt(coef.var)
  }
  colnames(varcovmat) <- colnames(mXadj)
  rownames(varcovmat) <- colnames(mXadj)
  vcov.mean <- varcovmat
  t.stat <- est.m$coefficients/s.e.
  p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

  mean.results <- as.data.frame(cbind(est.m$coefficients, s.e., t.stat, p.val))
  colnames(mean.results) <- c("coef", "std.error", "t-stat", "p-value")
  rownames(mean.results) <- colnames(mX)
}

#### VOLATILITY ###################

e.n <- length(resids)

#adjust no. of rows of vxreg-matrix:
if(is.null(vxreg)){
  vxreg <- NULL
}else{
  if(y.n!=e.n){
    vxreg.names <- colnames(vxreg)
    vxreg <- cbind(vxreg[I(y.n-e.n+1):y.n,])
    colnames(vxreg) <- vxreg.names
  }
}

#estimate:
vX <- regs.var(resids, vc=TRUE, arch=arch, asym=asym,
  log.ewma=log.ewma, vxreg=vxreg, zero.adj=zero.adj)
ewma.vol.chk <- if(is.null(log.ewma)){0}else{
  ifelse(is.null(log.ewma$lag), 1, log.ewma$lag)
}
pstar <- max(arch, asym, ewma.vol.chk)
logepadj <- as.vector(vX[I(pstar+1):e.n, 1])
vXadj <- cbind(vX[c(pstar+1):e.n, -1])
est.v <- ols(logepadj, vXadj, tol=tol, LAPACK=LAPACK,
  method=2)
fit.v <- as.vector(vXadj%*%cbind(est.v$coefficients))
ustar <- logepadj-fit.v
ustar2 <- ustar^2
vXadj.n <- NROW(vXadj)
vXadj.k <- NCOL(vXadj)
d.f.v. <- vXadj.n - vXadj.k
sigma2.v <- sum(ustar2)/d.f.v.

varcovmat.v <- sigma2.v*est.v$xtxinv
s.e. <- sqrt(as.vector(diag(varcovmat.v)))
colnames(varcovmat.v) <- colnames(vX)[-1]
rownames(varcovmat.v) <- colnames(vX)[-1]
vcov.var <- varcovmat.v
t.stat <- est.v$coefficients/s.e.
p.val <- pt(abs(t.stat), d.f.v., lower.tail=FALSE)*2

if(vc.adj){
  Elnz2 <- -log(mean(exp(ustar)))
  t.stat[1] <- ((est.v$coefficients[1]-Elnz2)^2)/s.e.[1]^2
  p.val[1] <- pchisq(t.stat[1], 1, lower.tail=FALSE)
  est.v$coefficients[1] <- est.v$coefficients[1] - Elnz2
}
fit.v <- exp(fit.v - Elnz2)
resids.std <- resids[I(pstar+1):e.n]/sqrt(fit.v)

variance.results <- as.data.frame(cbind(est.v$coefficients, s.e., t.stat, p.val))
colnames(variance.results) <- c("coef", "std.error", "t-stat", "p-value")
rownames(variance.results) <- colnames(vX)[-1]

### DIAGNOSTICS #################

if(verbose){
  diagnostics <- matrix(NA, 5, 3)
  colnames(diagnostics) <- c("Chi-sq", "df", "p-value")
  rownames(diagnostics) <- c(paste("Ljung-Box AR(", qstat.options[1], ")", sep=""),
    paste("Ljung-Box ARCH(", qstat.options[2], ")", sep=""),
    "Skewness", "JB-test", "R-squared")
  ar.LjungBox <- Box.test(resids.std, lag = qstat.options[1], type="L")
  diagnostics[1,1] <- ar.LjungBox$statistic
  diagnostics[1,2] <- qstat.options[1]
  diagnostics[1,3] <- ar.LjungBox$p.value
  arch.LjungBox <- Box.test(resids.std^2, lag = qstat.options[2], type="L")
  diagnostics[2,1] <- arch.LjungBox$statistic
  diagnostics[2,2] <- qstat.options[2]
  diagnostics[2,3] <- arch.LjungBox$p.value
  skew.test <- skewness.test(resids.std)
  diagnostics[3,1] <- skew.test$statistic
  diagnostics[3,2] <- 1
  diagnostics[3,3] <- skew.test$p.value
  normality.test <- jb.test(resids.std)
  diagnostics[4,1] <- normality.test$statistic
  diagnostics[4,2] <- 2
  diagnostics[4,3] <- normality.test$p.value

  ##R-squared:
  TSS <- sum( (vYadj - mean(vYadj))^2 )
  RSS <- sum( (resids - mean(resids))^2 )
  Rsquared <- 1 - RSS/TSS
  diagnostics[5,1] <- Rsquared
}

### OUTPUT: ######################

#mean:
add.nas2mean <- length(y) - length(resids)
if(verbose){
  if(!is.null(mX)){
    out$mean.fit <- zoo(c(rep(NA, add.nas2mean), fit.m), order.by=zoo.index.y)
  }
}
resids <- zoo(c(rep(NA, add.nas2mean), resids), order.by=zoo.index.y)
out$resids <- resids

#volatility:
add.nas2var <- length(y) - length(ustar)
if(verbose){
  fit.v <- zoo(c(rep(NA, add.nas2var), fit.v), order.by=zoo.index.y)
  out$var.fit <- fit.v
  ustar <- zoo(c(rep(NA, add.nas2var), ustar), order.by=zoo.index.y)
  out$resids.ustar <- ustar
}
resids.std <- zoo(c(rep(NA, add.nas2var), resids.std), order.by=zoo.index.y)
out$resids.std <- resids.std
if(vc.adj){
  out$Elnz2 <- Elnz2
}

#results:
out$vcov.mean <- vcov.mean
out$vcov.var <- vcov.var
out$mean.results <- mean.results
out$variance.results <- variance.results
if(verbose){
  out$diagnostics <- diagnostics
}

out <- c(list(date=date()), out)
class(out) <- "arx"
return(out)
}
