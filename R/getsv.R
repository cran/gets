getsv <-
function(object, keep=c(1), t.pval=0.05, do.pet=TRUE,
  wald.pval=0.05, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025),
  info.method=c("sc", "aic", "hq"), include.empty=FALSE,
  zero.adj=NULL, vc.adj=NULL, tol=NULL, LAPACK=NULL,
  max.regs=1e+05, verbose=TRUE, alarm=FALSE)
{
### ARGUMENTS ###########

info.method <- match.arg(info.method)
vc=TRUE #obligatory
keep.groups=NULL

##zoo:
#if(class(object)=="arx"){ e <- na.trim(object$resids) }
#if(class(object)=="gets"){ e <- na.trim(object$resids) }
e <- na.trim(object$resids)
index.e <- index(e)
e <- coredata(e)
e.n <- length(e)

arch <- object$aux$arch
asym <- object$aux$asym
log.ewma <- object$aux$log.ewma

vxreg <- object$aux$vxreg
if(!is.null(vxreg)){
  vxreg <- as.zoo(vxreg)
  vxreg <- window(vxreg, start=index.e[1],
    end=index.e[e.n])
  vxreg <- cbind(coredata(vxreg))
}

#variance-covariance matrix:
vcov.type <- "ordinary"

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
notes <- list()
spec <- list()
spec.results <- NULL

## CREATE REGRESSOR MATRIX: ########################

mXorig <- regs.var(e, vc=1, arch=arch, log.ewma=log.ewma,
  asym=asym, vxreg=vxreg, zero.adj=zero.adj)
mX.colnames <- colnames(mXorig[,-1])
if(!is.null(log.ewma)){ewma.chk <- 1}else{ewma.chk <- 0}
t1 <- max(arch, asym, ewma.chk) + 1
eadj <- e[t1:e.n]
logep <- mXorig[t1:e.n, 1]
mX <- cbind(mXorig[t1:e.n, -1])

## GUM: ############################################

if(!(1 %in% keep)){
  notes <- c(notes, "Regressor 1 set to 'keep'")
  keep <- union(keep, 1)
}
keep.n <- length(keep)
gum.n <- NCOL(mX)
gum <- 1:gum.n
delete <- setdiff(gum, keep)
delete.n <- length(delete)

#deletable and non-deletable regressors:
if(delete.n > 0){mXdel <- cbind(mX[,delete])}else{mXdel <- NULL}
if(is.null(keep)){mXndel <- NULL}else{mXndel <- cbind(mX[,keep])}

#estimate GUM:
mXadj <- cbind(mXdel,mXndel)
est <- ols(logep, mXadj, tol = tol, LAPACK=LAPACK,
  method=2)
fit <- as.vector(mXadj%*%cbind(est$coefficients))
ustar <- as.vector(logep) - fit
ustar2 <- ustar^2
mXadj.n <- NROW(mXadj)
mXadj.k <- NCOL(mXadj)
d.f. <- mXadj.n - mXadj.k
sumUstar2 <- sum(ustar2)
sigma2 <- sumUstar2/d.f.

#estimate s.e.; compute t-stats. and p-vals.:
if(vcov.type == "ordinary"){
  varcovmat <- sigma2*est$xtxinv
  coef.var <-as.vector(diag(varcovmat))
  s.e. <- sqrt(coef.var)
}
t.stat <- est$coefficients/s.e.
p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

#make standardised residuals {z_hat}:
Elnz2 <- -log(mean(exp(ustar)))
vconstadj <- -Elnz2
sigma2adj <- exp(fit + vconstadj)
zhat <- eadj/sqrt(sigma2adj)
out$gum.resids.std <- zoo(c(rep(NA, c(e.n-length(zhat))), zhat), order.by=index.e)

#adjust variance constant and its p-value:
if(vc.adj==TRUE){
  where.const <- which(union(delete,keep)==1)
  waldstat <- ((est$coefficients[where.const] + vconstadj)^2)/s.e.[where.const]^2
  p.val[where.const] <- pchisq(waldstat, 1, lower.tail = FALSE)
  est$coefficients[where.const] <- est$coefficients[where.const] + vconstadj
}

#make diagnostics table:
if(verbose == TRUE){
  diagnostics <- matrix(NA, 2, 2)
  colnames(diagnostics) <- c("Chi-sq", "p-value")
  rownames(diagnostics) <- c(paste("Ljung-Box AR(", ar.LjungB[1],
    ")", sep=""), paste("Ljung-Box ARCH(", arch.LjungB[1], ")",
    sep=""))
}

#Ljung-Box test for serial correlation in {z_t}:
if(!is.null(ar.LjungB)){
  ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
  if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.gum.chk <- 0}else{ar.gum.chk <- 1}
  if(verbose == TRUE){
    diagnostics[1,1] <- ar.LjungBox$statistic
    diagnostics[1,2] <- ar.LjungBox$p.value}
}else{ar.gum.chk <- 1}

#Ljung-Box test for arch in {z_t^2}
if(!is.null(arch.LjungB)){
  arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
  if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.gum.chk <- 0}else{arch.gum.chk <- 1}
  if(verbose == TRUE){
    diagnostics[2,1] <- arch.LjungBox$statistic
    diagnostics[2,2] <- arch.LjungBox$p.value}
}else{arch.gum.chk <- 1}

#results:
if(verbose){
  keep.labels <- c(rep(0,delete.n), rep(1,keep.n))
#  if(!is.null(keep.groups)){
#    keep.groups.n <- length(keep.groups)
#    if(keep.groups.n < 10){divisor <- 10}else{divisor <- 100}
#    for(g in 1:length(keep.groups)){
#      keep.groups.g <- keep.groups[[g]]
#      keep.groups.g.n <- length(keep.groups.g)
#      tmp.which <- which(c(delete,keep) %in% keep.groups.g[2:keep.groups.n])
#      keep.labels[tmp.which] <- 1 + g/divisor
#    }
#  }

  results <- as.data.frame(cbind(c(delete,keep), keep.labels,
    est$coefficients, s.e., t.stat, p.val))
  colnames(results) <- c("reg.no", "keep", "coef", "std.error", "t-stat", "p-value")
  rownames(results) <- mX.colnames[c(delete,keep)] #NULL
  out$gum.variance <- results
  out$gum.diagnostics <- diagnostics
} #end if(verbose)

#if GUM passes diagnostic checks:
if((ar.gum.chk*arch.gum.chk) != 0){

  spec[[1]] <- gum

  #specification results
  logl <- -mXadj.n*log(2*pi)/2 - sum(log(sigma2adj))/2 - sum(eadj^2/sigma2adj)/2
  info.results <- info.criterion(logl, mXadj.n, mXadj.k, method=info.method)
  spec.results <- rbind( c(info.results$value, logl,
    info.results$n, info.results$k) )
  col.labels <- c(paste("info(", info.method, ")", sep=""),
    "LogL", "n", "k")
#OLD:
#  spec.results <- rbind( c(info.results$value, info.results$n, info.results$k) )
#  col.labels <- c(paste("info(", info.method, ")", sep=""), "n", "k")
  row.labels <- c("spec1 (gum)")

  #record data for Wald-tests against gum:
  gum.regs <- c(delete, keep)
  gum.coefs <- est$coefficients
  gum.varcovmat <- varcovmat

}else{
  notes <- c(notes, c("GUM does not pass one or more diagnostic checks"))
}

##adjust variance constant and its p-value:
#if(vc.adj==TRUE){
#  where.const <- which(union(delete,keep)==1)
#  waldstat <- ((est$coefficients[where.const] + vconstadj)^2)/s.e.[where.const]^2
#  p.val[where.const] <- pchisq(waldstat, 1, lower.tail = FALSE)
#  est$coefficients[where.const] <- est$coefficients[where.const] + vconstadj
#}

## EMPTY MODEL: ####################################

if( ar.gum.chk*arch.gum.chk!=0 && delete.n>0 && include.empty==TRUE && keep.n>0 ){

  #estimate benchmark model:
  est <- ols(logep, mXndel, tol = tol, LAPACK=LAPACK,
    method=2)
  fit <- as.vector(mXndel%*%cbind(est$coefficients))
  ustar <- as.vector(logep) - fit
  ustar.n <- length(ustar)
  mXndel.k <- NCOL(mXndel)
  ustar2 <- ustar^2; sumUstar2 <- sum(ustar2)
  sigma2 <- sumUstar2/(ustar.n-mXndel.k)

  #make standardised residuals {z_hat}:
  Elnz2 <- -log(mean(exp(ustar)))
  vconstadj <- -Elnz2
  sigma2adj <- exp(fit + vconstadj)
  zhat <- eadj/sqrt(sigma2adj)

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

    spec[[length(spec)+1]] <- keep

    #specification results
    logl <- -ustar.n*log(2*pi)/2 - sum(log(sigma2adj))/2 - sum(eadj^2/sigma2adj)/2
    info.results <- info.criterion(logl, ustar.n, mXndel.k, method=info.method)

    #compute wald-test:
#    if(do.pet){
#      mR <- NULL
#      for(k in 1:gum.n){
#        if(gum.regs[k] %in% delete){
#          mR <- rbind(mR, c(rep(0,I(k-1)), 1, rep(0, I(gum.n-k) )))
#        }
#      }
#      mRestq <- mR%*%cbind(gum.coefs)
#      wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=tol)%*%mRestq
#      pet.chk <- as.logical(wald.pval < pchisq(wald.stat, delete.n, lower.tail = FALSE))
#    }else{pet.chk <- TRUE} #end if(do.pet)

    #check pet result:
#    if(pet.chk){
      spec.results <- rbind(spec.results, c(info.results$value,
        logl, info.results$n, info.results$k))
      row.labels <- c(row.labels, paste("spec", length(spec), " (empty)", sep=""))
#    }else{notes <- c(notes, c("Empty model fails PET against VGUM"))}

  }else{
    notes <- c(notes, c("Empty log-variance model does not pass one or more diagnostic checks"))
  } #end if(empty passes diagnostics==TRUE){..}else{..}
} #end if(include empty model==TRUE)

## MULTI-PATH SEARCH: #################

insig.regs <- NULL
paths <- list()
if( ar.gum.chk*arch.gum.chk!=0 && delete.n>1 ){

  #paths:
  insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
  n.paths <- length(insig.regs)
  paths <- list()

  #if paths = 0:
  if(n.paths == 0){
    notes <- c(notes, c("No insignificant regressors in VGUM"))
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
        mXndell <- mX[,keep.adj]

        #estimate model:
        mXadj <- cbind(mXdell,mXndell)
        est <- ols(logep, mXadj, tol = tol, LAPACK=LAPACK,
          method=2)
        fit <- as.vector(mXadj%*%cbind(est$coefficients))
        ustar <- as.vector(logep) - fit
        ustar2 <- ustar^2
        mXadj.k <- ncol(mXadj)
        d.f. <- mXadj.n - mXadj.k
        sumUstar2 <- sum(ustar2)
        sigma2 <- sumUstar2/d.f.

        #make standardised residuals {z_hat}:
        Elnz2 <- -log(mean(exp(ustar)))
        vconstadj <- -Elnz2
        sigma2adj <- exp(fit + vconstadj)
        zhat <- eadj/sqrt(sigma2adj)

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

        #if either ar.chk or arch.chk is unequal to 1, then move
        #path[length(path)] to keep.adj
        if((ar.chk*arch.chk) == 0){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*I(-1))
          next #next j
        }

#        #check keep.groups:
#        if(!is.null(keep.groups)){
#          skip.chk <- FALSE
#          for(g in 1:length(keep.groups)){
#            group.g.n <- length(keep.groups[[g]])
#            group.g <- keep.groups[[g]][2:group.g.n]
#            vars.in.delete <- intersect(group.g, delete.adj)
#            vars.in.delete.n <- length(vars.in.delete)
#            vars.in.keep <- intersect(group.g, keep.adj)
#            vars.in.keep.n <- length(vars.in.keep)
#            if(vars.in.delete.n + vars.in.keep.n == keep.groups[[g]][1]){
#              if(vars.in.delete.n > 0){
#                keep.adj <- union(vars.in.delete, keep.adj)
#                delete.adj <- setdiff(delete.adj, vars.in.delete)
#                path <- union(path, vars.in.delete*I(-1))
#                skip.chk <- TRUE
#                break  ##break for loop with index g
#              }
#            }
#          }
#          if(skip.chk == TRUE){next}  ##next j
#        } #end check keep.groups

        #if ar.chk*arch.chk==1:
        if(ar.chk*arch.chk == 1){

          #stop if no more deletable regressors:
          if(length(delete.adj)==0){
            spec.adj <- keep.adj
            break
          } #end if(length(..)==0)

          #estimate s.e.; compute t-stats. and p-vals.:
          if(vcov.type == "ordinary"){
            coef.var <-as.vector(sigma2*diag(est$xtxinv))
            s.e. <- sqrt(coef.var)
          }
          t.stat <- est$coefficients/s.e.
          p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

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
          } #end if..else.. p-value > t.pval, then delete regressor

        } #end if ar.chk*arch.chk==1

      } #### end single-path search

      #add path to the paths list:
      paths[[length(paths)+1]] <- path

      #check if spec.adj is already in spec:
      n.spec <- length(spec)
      for(l in 1:n.spec){
        chk.spec <- setequal(spec.adj, spec[[l]])
        if(chk.spec==TRUE){break}
      }

      #if spec.adj not in spec (among terminals):
      if(chk.spec==FALSE){

        #add spec.adj to spec:
        spec[[length(spec)+1]] <- spec.adj

        #specification results
        n.spec.adj <- length(spec.adj)
        logl <- -mXadj.n*log(2*pi)/2 - sum(log(sigma2adj))/2 - sum(eadj^2/sigma2adj)/2
        info.results <- info.criterion(logl, mXadj.n, n.spec.adj, method=info.method)

        #add terminal to spec:
        spec.results <- rbind(spec.results, c(info.results$value,
          logl, info.results$n, info.results$k))
        row.labels <- c(row.labels, paste("spec", length(spec), sep=""))

      } #end if(chk.spec==FALSE)
    } #end multi-path search: for(i in 1:n.paths) loop
  } #end if(n.paths > 0)
} #end if(ar/arch.chk and delete.n>1)

## FIND THE BEST MODEL: ########################

if(verbose){
  if(!is.null(spec.results)){

#    #if keep.groups != NULL:
#    if(!is.null(keep.groups)){
#      if((ar.bench.chk*arch.bench.chk) == 1){
#        J <- I(1:nrow(spec.results))[-2]
#        models <- cbind(J, spec.results[-2,])
#        colnames(models) <- NULL
#        notes <- c(notes, "keep.groups not NULL: Empty model (spec2) not included in selection")
#      }else{
#        J <- 1:nrow(spec.results)
#        models <- cbind(J, spec.results)
#        colnames(models) <- NULL
#      }
#    } #end if(!is.null(keep.groups))

    #if keep.groups == NULL
    if(is.null(keep.groups)){
      J <- 1:nrow(spec.results)
      models <- cbind(J, spec.results)
      colnames(models) <- NULL
    } #end if(is.null(keep.groups))

    #find the best model:
    where <- which.min(models[,2])
    best.spec <- spec[[where]] #winner

    #check for several minimums:
    min.models <- min(models[,2])
    wheres <- which(models[,2]==min.models)
    if(length(wheres)>1){notes <- c(notes, "Several terminal specifications attain the minimum information criterion")}

  } #end if(!is.null(spec.results))
} #end if(verbose)

## ESTIMATE SPECIFIC MODEL #############

specific.results <- NULL
specific.diagnostics <- NULL
#OLD: if(verbose){..
if(!is.null(spec.results)){

  #estimate specific:
  specific <- sort(best.spec)  ##specific
  mXadj <- cbind(mX[,specific])
  est <- ols(logep, mXadj, tol = tol, LAPACK=LAPACK,
    method=2)
  fit <- as.vector(mXadj%*%cbind(est$coefficients))
  ustar <- as.vector(logep) - fit
  ustar2 <- ustar^2
  mXadj.k <- ncol(mXadj)
  d.f. <- mXadj.n - mXadj.k
  sumUstar2 <- sum(ustar2)
  sigma2 <- sumUstar2/d.f.

  #estimate s.e.; compute t-stats. and p-vals.:
  if(vcov.type == "ordinary"){
    varcovmat <- sigma2*est$xtxinv
    s.e. <- sqrt(as.vector(diag(varcovmat)))
#OLD:
#      coef.var <-as.vector(sigma2*diag(est$xtxinv))
#      s.e. <- sqrt(coef.var)
  }
  colnames(varcovmat) <- colnames(mXadj)
  rownames(varcovmat) <- colnames(mXadj)
  vcov.var <- varcovmat
  t.stat <- est$coefficients/s.e.
  p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

  #resids, sigma2.fit, zhat:
  Elnz2 <- -log(mean(exp(ustar)))
  vconstadj <- -Elnz2
  sigma2adj <- exp(fit + vconstadj)
  zhat <- eadj/sqrt(sigma2adj)
  resids.std <- zoo(c(rep(NA, I(e.n-length(zhat))), zhat), order.by=index.e)
  volatility <- zoo(c(rep(NA, I(e.n-length(zhat))), sigma2adj), order.by=index.e)
#OLD:
#  resids <- zoo(eadj, order.by=index.e)
  ustar <- zoo(c(rep(NA, I(e.n-length(ustar))), ustar), order.by=index.e)
  out <- c(list(var.fit=volatility,
    resids.ustar=ustar, resids.std=resids.std, Elnz2=Elnz2,
    vcov.var=vcov.var), out)
#OLD:
#  out <- c(list(resids=resids, var.fit=volatility,
#    resids.ustar=ustar, resids.std=resids.std, Elnz2=Elnz2,
#    vcov.var=vcov.var), out)

  #adjust log-variance constant and its p-value:
  if( vc.adj==TRUE && 1%in%specific ){
    where.const <- which(specific==1)
    waldstat <- ((est$coefficients[where.const] + vconstadj)^2)/s.e.[where.const]^2
    p.val[where.const] <- pchisq(waldstat, 1, lower.tail = FALSE)
    est$coefficients[where.const] <- est$coefficients[where.const] + vconstadj
  }

  #make results table:
  specific.results <- as.data.frame(cbind(specific, est$coefficients, s.e., t.stat, p.val))
  colnames(specific.results) <- c("reg.no", "coef", "std.error", "t-stat", "p-value")
  rownames(specific.results) <- mX.colnames[specific] #NULL

  #make diagnostics table:
  if(!is.null(ar.LjungB) || !is.null(arch.LjungB)){
    specific.diagnostics <- matrix(NA, 2, 2)
    colnames(specific.diagnostics) <- c("Chi-sq", "p-value")
    rownames(specific.diagnostics) <- c(paste("Ljung-Box AR(", ar.LjungB[1],
      ")", sep=""), paste("Ljung-Box ARCH(", arch.LjungB[1], ")",
      sep=""))
    if(!is.null(ar.LjungB)){
      ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
      specific.diagnostics[1,1] <- ar.LjungBox$statistic
      specific.diagnostics[1,2] <- ar.LjungBox$p.value
    }
    if(!is.null(arch.LjungB)){
      arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
      specific.diagnostics[2,1] <- arch.LjungBox$statistic
      specific.diagnostics[2,2] <- arch.LjungBox$p.value
    }
  } #end make diagnostics table
} #end if(!is.null(spec.results))

## OUT LIST: ################################

out$keep <- keep
out$insigs.in.gum <- insig.regs
out$resids <- zoo(eadj, order.by=index.e)

if( ar.gum.chk*arch.gum.chk != 0){
  if(verbose){ if(length(paths)==0){out$paths <- NULL}else{out$paths <- paths} }
  out$terminals <- spec
  colnames(spec.results) <- col.labels
  rownames(spec.results) <- row.labels
  out$terminals.results <- spec.results
  if(verbose){
    out$specific.variance <- specific.results
    out$specific.diagnostics <- specific.diagnostics
  }
} #end if((ar.gum.chk*arch.gum.chk) != 0)

if(length(notes) > 0){
  out$notes <- notes
}

out <- c(list(date=date()), out)
class(out) <- "gets"

if(alarm){ alarm() }
return(out)
}
