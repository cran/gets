getsm <-
function(object, keep=NULL, vcov.type=NULL,
  t.pval=0.05, do.pet=TRUE, wald.pval=0.05,
  ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL,pval=0.025),
  info.method=c("sc", "aic", "hq"), include.empty=FALSE,
  zero.adj=NULL, vc.adj=NULL, tol=NULL, LAPACK=NULL,
  max.regs=100000, verbose=TRUE, alarm=FALSE)
{
### ARGUMENTS: ###########

#match arguments:
info.method <- match.arg(info.method)
keep.groups=NULL

y <- object$aux$y
zoo.index.y <- index(y)
y <- coredata(y)
y.n <- length(y)
t1 <- zoo.index.y[1]
t2 <- zoo.index.y[y.n]

mc <- object$aux$mc
ar <- object$aux$ar
ewma <- object$aux$ewma

mxreg <- object$aux$mxreg
if(!is.null(mxreg)){
  mxreg <- as.zoo(mxreg)
  mxreg <- window(mxreg, start=t1, end=t2)
  mxreg <- cbind(coredata(mxreg))
}

arch <- object$aux$arch
asym <- object$aux$asym
log.ewma <- object$aux$log.ewma

vxreg <- object$aux$vxreg
if(!is.null(vxreg)){
  vxreg <- as.zoo(vxreg)
  vxreg <- window(vxreg, start=t1, end=t2)
  vxreg <- cbind(coredata(vxreg))
}

if(is.null(vcov.type)){vcov.type <- object$aux$vcov.type}
if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
  ar.LjungB$lag <- object$aux$qstat.options[1]
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
}
if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
  arch.LjungB$lag <- object$aux$qstat.options[2]
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])
}

if(is.null(zero.adj)){ zero.adj <- object$aux$zero.adj }
if(is.null(vc.adj)){ vc.adj <- object$aux$vc.adj }
if(is.null(tol)){ tol <- object$aux$tol }
if(is.null(LAPACK)){ LAPACK <- object$aux$LAPACK }

### INITIALISE ##########

out <- list()
out$call <- sys.call()
out$aux <- object$aux
notes <- list()
spec <- list()
spec.results <- NULL

#check if variance equation:
if(!is.null(arch) || !is.null(asym) || !is.null(log.ewma) || !is.null(vxreg)){
  var.spec.chk <- TRUE
}else{var.spec.chk <- FALSE}

## GUM: ############################################

##mean regressors:
ewma.mean.chk <- if(is.null(ewma)){0}else{ifelse(is.null(ewma$lag),1,ewma$lag)}
max.ar <- if( is.null(ar) && is.null(ewma) ){0}else{ max(ar,ewma.mean.chk) }
mXunadj <- regs.mean(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg)
y.n <- length(y)
if(is.null(mXunadj)){
  stop("Mean equation empty")
}else{
  mX <- cbind(mXunadj[I(max.ar+1):y.n,])
}
yadj <- y[I(max.ar+1):y.n]

#misc.:
yadj.n <- length(yadj)
keep.n <- length(keep)
gum.n <- ncol(mX)
gum <- 1:gum.n
delete <- setdiff(gum, keep)
delete.n <- length(delete)

#deletable and non-deletable regressors:
if(delete.n > 0){mXdel <- cbind(mX[,delete])}else{mXdel <- NULL}
if(is.null(keep)){mXndel <- NULL}else{mXndel <- cbind(mX[,keep])}

#estimate GUM:
mXadj <- cbind(mXdel,mXndel)
est <- ols(yadj, mXadj, tol = tol, LAPACK=LAPACK,
  method=2)
fit <- as.vector(mXadj%*%cbind(est$coefficients))
resids <- as.vector(yadj) - fit
resids.n <- length(resids)
resids2 <- resids^2
mXadj.k <- NCOL(mXadj)
d.f. <- yadj.n - mXadj.k
sumResids2 <- sum(resids2)
sigma2 <- sumResids2/d.f.

#estimate s.e.; compute t-stats. and p-vals.:
if(vcov.type == "ordinary"){
  varcovmat <- sigma2*est$xtxinv
  coef.var <-as.vector(diag(varcovmat))
  s.e. <- sqrt(coef.var)
}
if(vcov.type == "white"){
  matResids2 <- matrix(0, yadj.n, yadj.n)
  diag(matResids2) <- resids2
  omega.hat <- t(mXadj)%*%matResids2%*%mXadj
  varcovmat <- est$xtxinv%*%omega.hat%*%est$xtxinv
  coef.var <- as.vector(diag(varcovmat))
  s.e. <- sqrt(coef.var)
}
t.stat <- est$coefficients/s.e.
p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

if(!is.null(vxreg)){
  vxreg.names <- colnames(vxreg)
  vxreg <- cbind(vxreg[c(yadj.n-resids.n+1):yadj.n,])
  colnames(vxreg) <- vxreg.names
}
est.var <- arx(resids, mc=NULL, arch=arch, asym=asym,
  log.ewma=log.ewma, vxreg=vxreg, zero.adj=zero.adj,
  vc.adj=vc.adj, tol=tol, LAPACK=LAPACK, verbose=TRUE)
sigma2.fit <- coredata(na.trim(est.var$var.fit))
resids.adj <- coredata(resids[c(resids.n-length(sigma2.fit)+1):resids.n])
resids.adj.n <- length(resids.adj)
zhat <- coredata(na.trim(est.var$resids.std))

#make diagnostics table:
diagnostics <- matrix(NA, 2, 3)
colnames(diagnostics) <- c("Chi-sq", "df", "p-value")
rownames(diagnostics) <- c(paste("Ljung-Box AR(", ar.LjungB[1],
  ")", sep=""), paste("Ljung-Box ARCH(", arch.LjungB[1], ")",
  sep=""))

#Ljung-Box test of zhat:
if(!is.null(ar.LjungB)){
  ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
  if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.gum.chk <- 0}else{ar.gum.chk <- 1}
  if(verbose == TRUE){
    diagnostics[1,1] <- ar.LjungBox$statistic
    diagnostics[1,2] <- ar.LjungB[1]
    diagnostics[1,3] <- ar.LjungBox$p.value
  }
}else{ar.gum.chk <- 1}

#Ljung-Box test of zhat^2:
if(!is.null(arch.LjungB)){
  arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
  if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.gum.chk <- 0}else{arch.gum.chk <- 1}
  if(verbose == TRUE){
    diagnostics[2,1] <- arch.LjungBox$statistic
    diagnostics[2,2] <- arch.LjungB[1]
    diagnostics[2,3] <- arch.LjungBox$p.value
  }
}else{arch.gum.chk <- 1}

#results:
if(verbose){
  keep.labels <- c(rep(0,delete.n), rep(1,keep.n))
  if(!is.null(keep.groups)){
    keep.groups.n <- length(keep.groups)
    if(keep.groups.n < 10){divisor <- 10}else{divisor <- 100}
    for(g in 1:length(keep.groups)){
      keep.groups.g <- keep.groups[[g]]
      keep.groups.g.n <- length(keep.groups.g)
      tmp.which <- which(c(delete,keep) %in% keep.groups.g[2:keep.groups.n])
      keep.labels[tmp.which] <- 1 + g/divisor
    }
  }

  results <- as.data.frame(cbind(c(delete,keep), keep.labels,
    est$coefficients, s.e., t.stat, p.val))
  colnames(results) <- c("reg.no", "keep", "coef", "std.error", "t-stat", "p-value")
  rownames(results) <- colnames(mXunadj[,c(delete,keep)]) #NULL
  out$gum.mean <- results
  out$gum.variance <- est.var$variance.results

  out$gum.diagnostics <- diagnostics
} #end if(verbose)

#if GUM passes diagnostic checks:
if((ar.gum.chk*arch.gum.chk) != 0){

  spec[[1]] <- gum

  #specification results
  logl <- -resids.adj.n*log(2*pi)/2 - sum(log(sigma2.fit))/2 - sum(resids.adj^2/sigma2.fit)/2
  info.results <- info.criterion(logl, resids.adj.n, mXadj.k, method = info.method)
  spec.results <- rbind( c(info.results$value, logl,
    info.results$n, info.results$k) )
  col.labels <- c(paste("info(", info.method, ")", sep=""),
    "LogL", "n", "k")
  row.labels <- c("spec1 (gum)")

  #record data for tests against gum:
  gum.regs <- c(delete, keep)
  gum.coefs <- est$coefficients
  gum.varcovmat <- varcovmat

}else{
  notes <- c(notes, c("MGUM does not pass one or more diagnostic checks"))
}

## EMPTY MODEL: ####################################

if( ar.gum.chk*arch.gum.chk!=0 && delete.n>0 && include.empty==TRUE ){

#DO NOT do pet in order to enable reality check:

  #benchmark model:
  if(keep.n==0){
    resids <- as.vector(yadj)
    resids.n <- length(resids)
    fit <- rep(0, yadj.n)
    mXndel.k <- 0
    resids2 <- resids^2; sumResids2 <- sum(resids2)
    sigma2 <- sumResids2/yadj.n
  }
  if(keep.n>0){
    est <- ols(yadj, mXndel, tol = tol, LAPACK=LAPACK,
      method=1)
    fit <- as.vector(mXndel%*%cbind(est$coefficients))
    resids <- as.vector(yadj) - fit
    mXndel.k <- NCOL(mXndel)
    resids2 <- resids^2; sumResids2 <- sum(resids2)
    sigma2 <- sumResids2/(yadj.n-mXndel.k)
  }

  #estimate sigma2.fit, resids.adj, zhat:
  est.var <- arx(resids, mc=NULL, arch=arch, asym=asym,
    log.ewma=log.ewma, vxreg=vxreg, zero.adj=zero.adj, vc.adj=vc.adj,
    tol=tol, LAPACK=LAPACK, verbose=TRUE)
  sigma2.fit <- coredata(na.trim(est.var$var.fit))
  resids.adj <- coredata(resids[c(resids.n-length(sigma2.fit)+1):resids.n])
  resids.adj.n <- length(resids.adj)
  zhat <- coredata(na.trim(est.var$resids.std))

  #Ljung-Box test for serial correlation in {z_t}:
  if(!is.null(ar.LjungB)){
    ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
    if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.bench.chk <- 0}else{ar.bench.chk <- 1}
  }else{ar.bench.chk <- 1}

  #Ljung-Box test for arch in {z_t}
  if(!is.null(arch.LjungB)){
    arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
    if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.bench.chk <- 0}else{arch.bench.chk <- 1}
  }else{arch.bench.chk <- 1}

  #check if empty model passes diagnostic checks:
  if((ar.bench.chk*arch.bench.chk) != 0){

    spec[[length(spec)+1]] <- if(is.null(keep)){0}else{keep}

    #specification results
    logl <- -resids.adj.n*log(2*pi)/2 - sum(log(sigma2.fit))/2 - sum(resids.adj^2/sigma2.fit)/2
    info.results <- info.criterion(logl, resids.adj.n, keep.n,
      method = info.method)

    #check pet result:
    spec.results <- rbind(spec.results,
      c(info.results$value, logl, info.results$n,
      info.results$k))
    row.labels <- c(row.labels, paste("spec", length(spec), " (empty)", sep=""))

  }else{
    notes <- c(notes, c("Empty mean model does not pass one or more diagnostic checks"))
  } #end if(empty passes diagnostics==TRUE){..}else{..}
} #end if(include empty model==TRUE)

## MULTI-PATH SEARCH: #################

insig.regs <- NULL
paths <- list()
if( ar.gum.chk*arch.gum.chk!=0 && delete.n>1 ){

  #paths:
  insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
  n.paths <- length(insig.regs)

  #if paths = 0:
  if(n.paths == 0){
    notes <- c(notes, c("No insignificant regressors in MGUM"))
  }

  #if paths > 0:
  if(n.paths > 0){

    #paths:
    for(i in 1:n.paths){

      #print path if verbose:
      if(verbose){
        cat(paste("Searching path no. ", i, " out of ",
          n.paths, "\n", sep=""))
#        writeLines(paste("Searching path no. ", i, " out of ", n.paths, sep=""))
#        print(paste("Searching path no. ", i, " out of ", n.paths, sep=""),
#          quote=FALSE)
      }

      #prepare single-path search:
      path <- insig.regs[i]
      delete.adj <- setdiff(delete, insig.regs[i])
      keep.adj <- as.numeric(keep)

      #single-path search:
      for(j in 1:max.regs){

        #matrices:
        mXdell <- if(length(delete.adj)==0){NULL}else{mX[,delete.adj]}
        mXndell <- if(is.null(keep.adj)){NULL}else{mX[,keep.adj]}

        #estimate model:
        mXadj <- cbind(mXdell,mXndell)
        if(!is.null(mXadj)){
          est <- ols(yadj, mXadj, tol = tol, LAPACK=LAPACK,
            method=2)
          fit <- as.vector(mXadj%*%cbind(est$coefficients))
        }else{
          fit <- rep(0, yadj.n)
        }
        resids <- as.vector(yadj) - fit
        resids2 <- resids^2
        mXadj.k <- ncol(mXadj)
        d.f. <- yadj.n - mXadj.k
        sumResids2 <- sum(resids2)
        sigma2 <- sumResids2/d.f.

        #make resids.adj, sigma2.fit, zhat:
        est.var <- arx(resids, mc=NULL, arch=arch, asym=asym,
          log.ewma=log.ewma, vxreg=vxreg, zero.adj=zero.adj,
          vc.adj=vc.adj, tol=tol, LAPACK=LAPACK, verbose=TRUE)
        sigma2.fit <- coredata(na.trim(est.var$var.fit))
        resids.adj <- coredata(resids[c(resids.n-length(sigma2.fit)+1):resids.n])
        resids.adj.n <- length(resids.adj)
        zhat <- coredata(na.trim(est.var$resids.std))

        #Ljung-Box test for serial correlation in {z_t}:
        if(!is.null(ar.LjungB)){
          ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
          if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.chk <- 0}else{ar.chk <- 1}
        }else{ar.chk <- 1}

        #Ljung-Box test for arch in {z_t}
        if(!is.null(arch.LjungB)){
          arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
          if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.chk <- 0}else{arch.chk <- 1}
        }else{arch.chk <- 1}

        #if either ar.chk or arch.chk is equal to 0, then move
        #path[length(path)] to keep.adj
        if((ar.chk*arch.chk) == 0){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*I(-1))
          next #next j
        }

        #check keep.groups:
        if(!is.null(keep.groups)){
          skip.chk <- FALSE
          for(g in 1:length(keep.groups)){
            group.g.n <- length(keep.groups[[g]])
            group.g <- keep.groups[[g]][2:group.g.n]
            vars.in.delete <- intersect(group.g, delete.adj)
            vars.in.delete.n <- length(vars.in.delete)
            vars.in.keep <- intersect(group.g, keep.adj)
            vars.in.keep.n <- length(vars.in.keep)
            if(vars.in.delete.n + vars.in.keep.n == keep.groups[[g]][1]){
              if(vars.in.delete.n > 0){
                keep.adj <- union(vars.in.delete, keep.adj)
                delete.adj <- setdiff(delete.adj, vars.in.delete)
                path <- union(path, vars.in.delete*I(-1))
                skip.chk <- TRUE
                break  ##break for loop with index g
              }
            }
          }
          if(skip.chk == TRUE){next}  ##next j
        } #end check keep.groups

        #if ar.chk*arch.chk==1:
        if(ar.chk*arch.chk == 1){

          #stop if no more deletable regressors:
          if(length(delete.adj)==0){
            spec.adj <- keep.adj
            break
          } #end if(length(..)==0)

          if(!is.null(mXadj)){

            #estimate s.e.; compute t-stats. and p-vals.:
            if(vcov.type == "ordinary"){
              coef.var <-as.vector(sigma2*diag(est$xtxinv))
              s.e. <- sqrt(coef.var)
            } #end "ordinary"

            if(vcov.type == "white"){
              matResids2 <- matrix(0, yadj.n, yadj.n)
              diag(matResids2) <- resids2
              omega.hat <- t(mXadj)%*%matResids2%*%mXadj
              varcovmat <- est$xtxinv%*%omega.hat%*%est$xtxinv
              coef.var <- as.vector(diag(varcovmat))
              s.e. <- sqrt(coef.var)
            } #end "white"

            #t-tests:
            t.stat <- est$coefficients/s.e.
            p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

          } #end if(!is.null(mXadj))

          #if any p-value > t.pval:
          if(sum(p.val[1:I(length(delete.adj))] > t.pval) > 0){
            reg.no <- which.max(p.val[1:I(length(delete.adj))])
            #do pet test:
            if(do.pet){
              deleted <- setdiff(delete, delete.adj[-reg.no])
              n.deleted <- length(deleted)
              mR <- NULL #initiate restriction matrix
              for(k in 1:gum.n){
                if(gum.regs[k] %in% deleted){
                  mR <- rbind(mR, c(rep(0,I(k-1)), 1, rep(0, I(gum.n-k) )))
                } #end if(gum.regs[k}..)
              } #end for(k in ..)
              mRestq <- mR%*%cbind(gum.coefs)
              wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=tol)%*%mRestq
              pet.chk <- as.logical(wald.pval < pchisq(wald.stat, n.deleted, lower.tail = FALSE))
            } #end if(do.pet)
            #if pet=TRUE then delete regressor, else move to keep:
            if(pet.chk){
              path <- union(path, delete.adj[reg.no])
              delete.adj <- delete.adj[-reg.no]
            }else{
              path <- union(path, delete.adj[reg.no]*I(-1))
              keep.adj <- union(delete.adj[reg.no], keep.adj)
              delete.adj <- delete.adj[-reg.no]
            } #end if(pet.chk)else{..}

          }else{
            spec.adj <- union(delete.adj, keep.adj)
            break
          } #end if..else.. any p-value > t.pval

        } #end if ar.chk*arch.chk==1

      } #### end single-path search

      #add path to the paths list:
      paths[[length(paths)+1]] <- path

      #check if spec.adj is already in spec:
      if(length(spec.adj)==0){spec.adj <- 0} #check if empty
      for(l in 1:length(spec)){
        chk.spec <- setequal(spec.adj, spec[[l]])
        if(chk.spec==TRUE){break}
      }

      #if spec.adj not in spec (among terminals):
      if(chk.spec==FALSE){

        #add spec.adj to spec:
        spec[[length(spec)+1]] <- spec.adj

        #specification results
        if(spec.adj[1]==0){
          n.spec.adj <- 0
        }else{
          n.spec.adj <- length(spec.adj)
        }
        logl <- -resids.adj.n*log(2*pi)/2 - sum(log(sigma2.fit))/2 - sum(resids.adj^2/sigma2.fit)/2
        info.results <- info.criterion(logl, resids.adj.n,
          n.spec.adj, method = info.method)

        #add terminal to spec:
        spec.results <- rbind(spec.results,
          c(info.results$value, logl, info.results$n,
          info.results$k))
        row.labels <- c(row.labels, paste("spec", length(spec), sep=""))

      } #end if(chk.spec==FALSE)
    } #end multi-path search: for(i in 1:n.paths) loop
  } #end if paths > 0
} #end if(ar/arch.chk and delete.n>1)

## FIND THE BEST MODEL: ########################

#if(verbose){
  if(!is.null(spec.results)){

    #if keep.groups != NULL:
    if(!is.null(keep.groups)){
      if((ar.bench.chk*arch.bench.chk) == 1){
        J <- I(1:NROW(spec.results))[-2]
        models <- cbind(J, spec.results[-2,])
        colnames(models) <- NULL
        notes <- c(notes, "keep.groups not NULL: Empty model (spec2) not included in selection")
      }else{
        J <- 1:nrow(spec.results)
        models <- cbind(J, spec.results)
        colnames(models) <- NULL
      }
    } #end if(!is.null(keep.groups))

    #if is.null(keep.groups):
    if(is.null(keep.groups)){
      J <- 1:NROW(spec.results)
      models <- cbind(J, spec.results)
      colnames(models) <- NULL
    } #end if(is.null(keep.groups))

    #find best model:
    where <- which.min(models[,2])
    best.spec <- spec[[where]] #winner

    #check for several minimums:
    min.models <- min(models[,2])
    wheres <- which(models[,2]==min.models)
    if(length(wheres)>1){notes <- c(notes, "Several terminal specifications attain the minimum information criterion")}

  } #end if(!is.null(spec.results))
#} #end if(verbose)

## ESTIMATE SPECIFIC MODEL: #############

specific.mean <- NULL
specific.diagnostics <- NULL
if(verbose){
  if(!is.null(spec.results)){

    if(best.spec[1]==0){
      resids <- yadj
      specific.mean <- "empty"
      Rsquared <- 0
      vcov.mean <- NULL
    }else{

      #specific model:
      specific <- sort(best.spec)  ##specific
      mXadj <- cbind(mX[,specific])
      est <- ols(yadj, mXadj, tol = tol, LAPACK=LAPACK,
        method=2)
      fit <- as.vector(mXadj%*%cbind(est$coefficients))
      resids <- as.vector(yadj) - fit
      resids2 <- resids^2
      mXadj.k <- ncol(mXadj)
      d.f. <- yadj.n - mXadj.k
      sumResids2 <- sum(resids2)
      sigma2 <- sumResids2/d.f.
      #R-squared:
      if(verbose){
        TSS <- sum( (yadj - mean(yadj))^2 )
        RSS <- sum( (resids - mean(resids))^2 )
        Rsquared <- 1 - RSS/TSS
      }

      #estimate s.e.; compute t-stats. and p-vals.:
      if(vcov.type == "ordinary"){
        varcovmat <- sigma2*est$xtxinv
        s.e. <- sqrt(as.vector(diag(varcovmat)))
#OLD:
#        coef.var <-as.vector(sigma2*diag(est$xtxinv))
#        s.e. <- sqrt(coef.var)
      }
      if(vcov.type == "white"){
        matResids2 <- matrix(0, yadj.n, yadj.n)
        diag(matResids2) <- resids2
        omega.hat <- t(mXadj)%*%matResids2%*%mXadj
        varcovmat <- est$xtxinv%*%omega.hat%*%est$xtxinv
        coef.var <- as.vector(diag(varcovmat))
        s.e. <- sqrt(coef.var)
      }
      colnames(varcovmat) <- colnames(mXadj)
      rownames(varcovmat) <- colnames(mXadj)
      vcov.mean <- varcovmat
      t.stat <- est$coefficients/s.e.
      p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

      #make results table:
      specific.mean <- as.data.frame(cbind(specific, est$coefficients, s.e., t.stat, p.val))
      colnames(specific.mean) <- c("reg.no", "coef", "std.error", "t-stat", "p-value")
      rownames(specific.mean) <- colnames(mXunadj)[specific.mean[,1]] #NULL

    } #end if(best.spec[1]==0)else{..}

    #mean.fit, resids, sigma2.fit, zhat:
    est.var <- arx(resids, mc=NULL, arch=arch, asym=asym,
      log.ewma=log.ewma, vxreg=vxreg, zero.adj=zero.adj,
      vc.adj=vc.adj, tol=tol, LAPACK=LAPACK, verbose=TRUE)
    vcov.var <- vcov.arx(est.var, spec="variance")
    specific.variance <- cbind(est.var$variance.results)
    var.fit <- coredata(est.var$var.fit)
    var.fit <- zoo(c(rep(NA, c(y.n-length(var.fit))), var.fit), order.by=zoo.index.y)
    zhat <- coredata(est.var$resids.std)
    zhat <- zoo(c(rep(NA, c(y.n-length(zhat))), zhat), order.by=zoo.index.y)
    resids <- zoo(c(rep(NA, c(y.n-length(resids))), resids), order.by=zoo.index.y)
    fit <- zoo(c(rep(NA, c(y.n-length(fit))), fit), order.by=zoo.index.y)
    out <- c(list(mean.fit=fit, resids=resids, var.fit=var.fit,
      resids.std=zhat, Elnz2=est.var$Elnz2,
      vcov.mean=vcov.mean, vcov.var=vcov.var), out)

    #make diagnostics table:
    if(!is.null(ar.LjungB) || !is.null(arch.LjungB)){
      specific.diagnostics <- matrix(NA, 3, 3)
      colnames(specific.diagnostics) <- c("Chi-sq", "df", "p-value")
      rownames(specific.diagnostics) <- c(paste("Ljung-Box AR(", ar.LjungB[1],
        ")", sep=""), paste("Ljung-Box ARCH(", arch.LjungB[1], ")",
        sep=""), "R-squared")
      if(!is.null(ar.LjungB)){
        ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
        specific.diagnostics[1,1] <- ar.LjungBox$statistic
        specific.diagnostics[1,2] <- ar.LjungB[1]
        specific.diagnostics[1,3] <- ar.LjungBox$p.value
      }
      if(!is.null(arch.LjungB)){
        arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
        specific.diagnostics[2,1] <- arch.LjungBox$statistic
        specific.diagnostics[2,2] <- arch.LjungB[1]
        specific.diagnostics[2,3] <- arch.LjungBox$p.value
      }
      ##R-squared:
      specific.diagnostics[3,1] <- Rsquared

    } #end make diagnostics table
  } #end if(!is.null(spec.results))
} #end if(verbose)

## CREATE OUT LIST ################################

out$keep <- keep
out$insigs.in.gum <- insig.regs

if((ar.gum.chk*arch.gum.chk) != 0){
  if(verbose){ if(length(paths)==0){out$paths <- NULL}else{out$paths <- paths} }
  out$terminals <- spec
  colnames(spec.results) <- col.labels
  where.empty <- which(spec.results[,"k"]==0)
  if(include.empty==FALSE && length(where.empty) > 0){
    row.labels[where.empty] <- paste(row.labels[where.empty],
      " (empty)", sep="")
  }
  rownames(spec.results) <- row.labels
  out$terminals.results <- spec.results
  out$specific.mean <- specific.mean
  out$specific.variance <- specific.variance
  out$specific.diagnostics <- specific.diagnostics
} #end if((ar.gum.chk*arch.gum.chk) != 0)

if(length(notes) > 0){
  out$notes <- notes
}
out <- c(list(date=date()), out)
class(out) <- "gets"

if(alarm){ alarm() }
return(out)
}
