getsm <-
function(object, t.pval=0.05, wald.pval=t.pval,
  vcov.type=NULL, do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc", "aic", "hq"),
  keep=NULL, include.gum=FALSE, include.empty=FALSE,
  max.paths=NULL, max.regs=NULL, zero.adj=NULL, vc.adj=NULL,
  verbose=TRUE, print.searchinfo=TRUE, estimate.specific=TRUE,
  plot=NULL, alarm=FALSE)
{
  ### ARGUMENTS: ###########

  info.method <- match.arg(info.method)

  ##determine no. of paths:
  if( !is.null(max.paths) && max.paths < 1){
    stop("'max.paths' cannot be smaller than 1")
  }

  ##determine vcov.type:
  if(is.null(vcov.type)){
    vcov.type <- object$aux$vcov.type
  }else{
    vcovTypes <- c("ordinary", "white", "newey-west")
    which.type <- charmatch(vcov.type, vcovTypes)
    vcov.type <- vcovTypes[which.type]
  }

  ##determine ar and arch lags:
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])

  ## max.regs, tol, LAPACK:
  if(is.null(max.regs)){ max.regs <- 10*object$aux$y.n }
  tol <- object$aux$tol
  LAPACK <- object$aux$LAPACK

  ##check if mean equation:
  if( is.null(object$aux$mX) ){ stop("Mean equation empty") }

  ##check if variance equation:
  ##in the future, check object$aux$vX instead?:
  if( is.null(object$variance.results )){
    var.spec.chk <- FALSE
  }else{
    var.spec.chk <- TRUE
  }

  ## check if estimator user-defined:
  if( !is.null(object$aux$user.estimator) ){
    if(estimate.specific){
      estimate.specific <- FALSE
      message("  'estimate.specific' set to FALSE")
    }
    if( is.null(plot) || identical(plot,TRUE) ){
      plot <- FALSE
      message("  'plot' set to FALSE")
    }
  }

  ### INITIALISE ##########

  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  messages <- NULL
  spec <- list() #terminal specifications
  spec.results <- NULL

  ##deletable, non-deletable regressors, re-organise:
  keep.n <- length(keep)
  gum <- 1:object$aux$mXncol
  delete <- setdiff(gum, keep) #integer(0) if empty
  delete.n <- length(delete)
  if(delete.n > 0){mXdel <- cbind(object$aux$mX[,delete])}else{mXdel <- NULL}
  if(is.null(keep)){mXndel <- NULL}else{mXndel <- cbind(object$aux$mX[,keep])}
  mXadj <- cbind(mXdel,mXndel)

  ## add keep column to results:
  tmp <- rep(0,object$aux$mXncol)
  if(!is.null(keep)){ tmp[keep] <- 1 }
  tmpdf <- cbind(tmp, object$mean.results)
  tmp <- 1:object$aux$mXncol
  tmpdf <- cbind(tmp, tmpdf)
  colnames(tmpdf)[1:2] <- c("reg.no", "keep")
  out$gum.mean <- tmpdf
  out$gum.variance <- object$variance.results
  out$gum.diagnostics <- object$diagnostics

  ## GUM: ##################

  ##estimate GUM:
  if( is.null(object$aux$user.estimator) ){

    ##usual ols:
    estMethod <- which( vcov.type==c("none", "none", "ordinary",
      "white", "newey-west") )
    est <- ols(object$aux$y, mXadj, tol=object$aux$tol,
      LAPACK=object$aux$LAPACK, method=estMethod,
      user.fun=NULL, user.options=NULL)
    est$std.residuals <- coredata(na.trim(object$std.residuals))
    est$logl <- object$logl
    if( !is.null(object$aux$loge2.n) ){
      est$n <- object$aux$loge2.n
    }

  }else{

    ##user-defined estimator:
    est <- do.call(object$aux$user.estimator$name,
      list(object$aux$y, mXadj), envir=.GlobalEnv)
    #delete?:
    if( is.null(est$vcov) && !is.null(est$vcov.mean) ){
      est$vcov <- est$vcov.mean
    }

  } #end if( is.null(..) )

  ##diagnostics:
  if( !is.null(est$residuals) ){
    gum.chk <- diagnostics(est, ar.LjungB=ar.LjungB,
      arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
      verbose=FALSE, user.fun=user.diagnostics)
  }else{
    gum.chk <- TRUE
  }

  ## if GUM passes diagnostic checks:
  if(gum.chk){

    spec[[1]] <- spec.gum <- gum #add gum to list of terminal specs

    ## vcov, t-stats, p-vals:
    stderrs <- sqrt(diag(est$vcov))
    t.stat <- est$coefficients/stderrs
    p.val <- pt(abs(t.stat), est$df, lower.tail=FALSE)*2

    ##specification results
    info.results <- info.criterion(est$logl, est$n, est$k,
      method=info.method)
    spec.results <- rbind( c(info.results$value, est$logl,
      info.results$n, info.results$k) )
    col.labels <- c(paste("info(", info.method, ")", sep=""),
      "logl", "n", "k")
    row.labels <- c("spec 1 (gum):")

    ##record data for tests against gum:
    gum.regs <- c(delete, keep)
    gum.coefs <- object$mean.results[gum.regs,1]
    gum.varcovmat <- est$vcov #OLD: varcovmat

  }else{
    messages <- paste(messages,
      "- MGUM does not pass one or more diagnostic checks", sep="")
  }

  ## FUTURE: ADD 1-CUT MODEL #########

  ## EMPTY MODEL: ################

  if( gum.chk && delete.n>0 && include.empty ){

    ## DO NOT do pet in order to enable reality check!

    ## estimate model:
    if( is.null(object$aux$user.estimator) ){
      est <- ols(object$aux$y, mXndel, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, method=estMethod)
    }else{
      est <- do.call(object$aux$user.estimator$name,
        list(object$aux$y, mXndel), envir=.GlobalEnv)
      #delete?:
      if( is.null(est$vcov) && !is.null(est$vcov.mean) ){
        est$vcov <- est$vcov.mean
      }
    } #end if( is.null(..) )

    #add var.fit, modify resids, resids.std, logl and n:
    if(var.spec.chk){
      residsAdj <- zoo(est$residuals, order.by=object$aux$y.index)
      est.var <- arx(residsAdj, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg,
        zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
        tol=object$aux$tol, LAPACK=object$aux$LAPACK, plot=FALSE)
      est$std.residuals <- coredata(na.trim(est.var$std.residuals))
      est$var.fit <- coredata(na.trim(est.var$var.fit))
      est$logl <- est.var$logl
      est$residuals <- est$residuals[c(object$aux$y.n-object$aux$loge2.n+1):object$aux$y.n]
      est$n <- length(est$residuals)
    }

    diagnostics.chk <- diagnostics(est, ar.LjungB=ar.LjungB,
      arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
      verbose=FALSE, user.fun=user.diagnostics)

    ##if empty model passes diagnostic checks:
    if(diagnostics.chk){

      ##add empty to spec:
      spec[[length(spec)+1]] <- if(is.null(keep)){0}else{keep}

      ##specification results
      info.results <- info.criterion(est$logl, est$n, est$k,
        method = info.method)
      spec.results <- rbind(spec.results,
        c(info.results$value, est$logl, info.results$n,
        info.results$k))
      row.labels <- c(row.labels,
        paste("spec ", length(spec), " (empty):", sep=""))

    }else{
        messages <- paste(messages,
          "- Empty mean model does not pass one or more diagnostic checks",
          sep="")
    } #end if(empty passes diagnostics==TRUE){..}else{..}

  } #end if(include empty model==TRUE)

## MULTI-PATH SEARCH: #################

#future:
#regs <- list()
#regs.info <- list()

#the following should probably be moved to just before
#"single path search", that is, just beneath "prepare
#single-path search":
#regs.current.path <- list() #the regressions of the current path
#regs.info.current.path <- list() #the regression info of current path
#regression info: list(which.path=??, where.in.path=??)

insig.regs <- NULL
paths <- list()
if( gum.chk && delete.n>1 ){

  ## no. of paths:
  insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
  if( !is.null(max.paths) ){
    if(max.paths < length(insig.regs)){
      pvalRanksInv <- rank( 1-p.val[insig.regs] )
      insig.regs <- delete[ pvalRanksInv <= max.paths ]
    }
  }
  n.paths <- length(insig.regs)

  ## if paths = 0:
  if(n.paths == 0){
    messages <- paste(messages,
      "- All regressors significant in GUM mean equation", sep="")
  }

  ## if paths > 0:
  if(n.paths > 0){

    if(print.searchinfo){
      message(n.paths, " path(s) to search")
      message("Searching: ", appendLF=FALSE)
    }

    ## paths:
    for(i in 1:n.paths){

      ## print searchinfo:
      if(print.searchinfo){
        newLine <- ifelse(i==n.paths, TRUE, FALSE)
        message(i, " ", appendLF=newLine)
      }

      ## prepare single-path search:
      path <- insig.regs[i]
      delete.adj <- setdiff(delete, insig.regs[i])
      keep.adj <- as.numeric(keep)

      ## single-path search of path i:
      for(j in 1:max.regs){

        ## matrices:
        mXdell <- if(length(delete.adj)==0){NULL}else{object$aux$mX[,delete.adj]}
        mXndell <- if(is.null(keep.adj)){NULL}else{object$aux$mX[,keep.adj]}
        mXadj <- cbind(mXdell,mXndell)
        mXadj.k <- NCOL(mXadj)

        ## estimate model:
        if( is.null(object$aux$user.estimator) ){
          est <- ols(object$aux$y, mXadj, tol=object$aux$tol,
            LAPACK=object$aux$LAPACK, method=estMethod,
            user.fun=NULL, user.options=NULL)
        }else{
          ##user-defined estimator:
          est <- do.call(object$aux$user.estimator$name,
            list(object$aux$y, mXadj), envir=.GlobalEnv)
          #delete?:
          if( is.null(est$vcov) && !is.null(est$vcov.mean) ){
            est$vcov <- est$vcov.mean
          }
        } #end if(is.null(..))else(..)

        #add var.fit, modify resids, resids.std, logl and n:
        if(var.spec.chk){
          residsAdj <- zoo(est$residuals, order.by=object$aux$y.index)
          est.var <- arx(residsAdj, vc=object$aux$vc,
            arch=object$aux$arch, asym=object$aux$asym,
            log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg,
            zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
            tol=object$aux$tol, LAPACK=object$aux$LAPACK, plot=FALSE)
          est$std.residuals <- coredata(na.trim(est.var$std.residuals))
          est$var.fit <- coredata(na.trim(est.var$var.fit))
          est$logl <- est.var$logl
          est$residuals <- est$residuals[c(object$aux$y.n-object$aux$loge2.n+1):object$aux$y.n]
          est$n <- length(est$residuals)
        }

        ##diagnostics:
        diagnostics.chk <- diagnostics(est, ar.LjungB=ar.LjungB,
          arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
          verbose=FALSE, user.fun=user.diagnostics)

        ## if diagnostics fails (FALSE),
        ## then move path[length(path)] to keep.adj
        if(!diagnostics.chk){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*c(-1))
          next #next j
        }

        #if diagnostics are ok (TRUE):
        if(diagnostics.chk){

          ## stop if no more deletable regressors:
          if( length(delete.adj)==0 ){
            spec.adj <- keep.adj
            break
          } #end if(length(..)==0)

          ## vcov, t-stats, p-vals:
#for the future:
#          if( is.null(est$vcov) ){
#            est$vcov <- vcovFun(est, method="ordinary")
#          }
          stderrs <- sqrt(diag(est$vcov))
          t.stat <- est$coefficients/stderrs
          p.val <- pt(abs(t.stat), est$df, lower.tail=FALSE)*2

# for new version: here is where I need to identify the coefficient,
# among those not in keep, with highest p-value. sketch?:
# gumRegs # regressors in gum, e.g. 1:NCOL(x)
# keep # the regressors to be kept, e.g. 1:2
# keepAdj # keep + those put in keep
# regsAdj # regressors left
# pvalRanks <- rank(p.val)
# pvalRanksAdj <- pvalRanks[-keepAdj]
# pvalMin <- which.min(pvalRanksAdj)
# reg.no <- keepAdj[ which(pvalRanks == pvalMin) ]

          ## if any p-value > t.pval:
          if(sum(p.val[1:c(length(delete.adj))] > t.pval) > 0){

            reg.no <- which.max(p.val[1:I(length(delete.adj))])

            ## do pet test (i.e. wald-test):
            if(do.pet){
              deleted <- setdiff(delete, delete.adj[-reg.no])
              n.deleted <- length(deleted)
#for the new version:
#              mR <- matrix(0, length(gum.regs), length(gum.regs))
#              diag(mR) <- 1
#              mR <- rbind(mR[deleted,])
              mR <- NULL #initiate restriction matrix
              for(k in 1:length(gum.regs)){
              #for(k in 1:object$aux$mXncol){
                if(gum.regs[k] %in% deleted){
                  mR <- rbind(mR, c(rep(0,I(k-1)), 1, rep(0, I(object$aux$mXncol-k) )))
                } #end if(gum.regs[k}..)
              } #end for(k in ..)

              mRestq <- mR%*%cbind(gum.coefs)
              wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=object$aux$tol)%*%mRestq
              pet.chk <- as.logical(wald.pval < pchisq(wald.stat, n.deleted, lower.tail = FALSE))
            }else{
              pet.chk <- TRUE
            } #end if(do.pet)else..

            ## delete regressor if(pet.chk), else move to keep:
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

        } #end if diagnostics are ok

      } #### end single-path search: for(j in..

      #it is probably at this point that I should introduce
      #the 'bookkeeping' with respect to regs, regs.info,
      #regs.current.path and regs.info.current.path

      #add path to the paths list:
      paths[[length(paths)+1]] <- path

      #check if spec.adj (terminal) is already in spec:
      if(length(spec.adj)==0){spec.adj <- 0} #check if completely empty
      for(l in 1:length(spec)){
        chk.spec <- setequal(spec.adj, spec[[l]])
        if(chk.spec==TRUE){break} #stop for(l in..)
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
#        if( is.null(object$aux$user.estimator) ){
#          est$logl <- -resids.adj.n*log(2*pi)/2 - sum(log(est$var.fit))/2 - sum(est$std.residuals^2)/2
#        }
        info.results <- info.criterion(est$logl, est$n,
          n.spec.adj, method=info.method)

        #add terminal to spec.results:
        spec.results <- rbind(spec.results,
          c(info.results$value, est$logl, info.results$n,
          info.results$k))
        row.labels <- c(row.labels, paste("spec ", length(spec), ":", sep=""))

      } #end if(chk.spec==FALSE)

    } #end multi-path search: for(i in 1:n.paths) loop

  } #end if paths > 0
} #end if( gum.chk!=0 && delete.n>1 )

  ## FIND THE BEST MODEL: ########################

  #future?: if include.gum=FALSE (default), then first check if
  #spec results is empty, then add gum to it if it is

  if(!is.null(spec.results)){

    J <- 1:NROW(spec.results)
    models <- cbind(J, spec.results)
    colnames(models) <- NULL

    #find best model and check for several minimums:
    if(include.gum){
      min.value <- min(models[,2])
      where <- which(min.value==models[,2])
    }else{
      if(length(spec)==1){
        where <- 1
      }else{
        min.value <- min(models[-1,2])
        where <- which(min.value==models[-1,2]) + 1
      } #end if(length(spec)==1)
    } #end if(include.gum)..
    if(length(where)>1){
      messages <- paste(messages,
        "- Several terminal specifications attain the minimum information criterion",
        sep="")
      }
    best.spec <- spec[[where[1]]] #winner

  } #end if(!is.null(spec.results))

  ## OUTPUT ################################

  out$keep <- keep
  #out$insigs.in.gum <- insig.regs #don't think I need this one...

  ##if no search has been undertaken:
  if(is.null(spec.results)){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }

  ##if search has been undertaken:
  if(!is.null(spec.results)){

    ##terminals results:
    if(length(paths)==0){
      out$paths <- NULL
    }else{ out$paths <- paths }
    out$terminals <- spec
    colnames(spec.results) <- col.labels
    where.empty <- which(spec.results[,"k"]==0)
    if(include.empty==FALSE && length(where.empty) > 0){
      row.labels[where.empty] <- paste("spec ", where.empty,
        " (empty):", sep="")
    }
    rownames(spec.results) <- row.labels
    out$terminals.results <- spec.results

    if(!estimate.specific){
      if(best.spec==0 || is.na(best.spec) || length(best.spec)==0 ){
        out$specific.spec <- NULL
      }else{
        specific <- sort(best.spec)
        names(specific) <- object$aux$mXnames[specific]
        out$specific.spec <- specific
      }
    } #end if(!estimate.specific)

    if(estimate.specific){

      ##prepare for estimation:
      yadj <- zoo(cbind(object$aux$y),
        order.by=object$aux$y.index)
      colnames(yadj) <- object$aux$y.name
      specific <- sort(best.spec)
      if(specific[1]==0){
        mXadj <- NULL
      }else{
        mXadj <- cbind(object$aux$mX[,specific])
        colnames(mXadj) <- object$aux$mXnames[specific]
        mXadj <- zoo(mXadj, order.by=object$aux$y.index)
      }
      if(is.null(object$aux$vxreg)){
        vxregAdj <- NULL
      }else{
        vxregAdj <- zoo(object$aux$vxreg,
          order.by=object$aux$y.index)
      }
      if(is.null(ar.LjungB)){
        ar.LjungB <- object$aux$qstat.options[1]
      }
      if(is.null(arch.LjungB)){
        arch.LjungB <- object$aux$qstat.options[2]
      }

      ##estimate specific model:
      est <- arx(yadj, mxreg=mXadj, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=vxregAdj,
        zero.adj=object$aux$zero.adj,
        vc.adj=object$aux$vc.adj, vcov.type=vcov.type,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        user.diagnostics=user.diagnostics, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, plot=FALSE)

      ##rename various stuff:
      est$call <- est$date <- NULL
      where.diagnostics <- which(names(est)=="diagnostics")
      if(length(where.diagnostics)>0){
        names(est)[where.diagnostics] <- "specific.diagnostics"
      }
      est$aux$y.name <- object$aux$y.name
      est <- unclass(est)
      names(specific) <- colnames(mXadj)
      out$specific.spec <- specific
      out <- c(out,est)
      ##solution to issue raised by J-dog?:
#      if(specific[1]==0){
#        out$aux$mX <- out$aux$mXnames <- NULL
#      }

    } #end if(estimate.specific)

  } #end if(!is.null(spec.results))

  if(!is.null(messages)){ out$messages <- messages }
  if( is.null(object$aux$user.estimator) ){
    out$aux$y <- object$aux$y
    out$aux$y.index <- object$aux$y.index
    out$aux$y.n <- object$aux$y.n
    out$aux$y.name <- object$aux$y.name
    out$aux$mXnames.gum <- object$aux$mXnames
    out$aux$call.gum <- object$call
    if(is.null(out$aux$vcov.type)){ out$aux$vcov.type <- vcov.type }
  }
  out <- c(list(date=date(), gets.type="getsm"), out)
  out$time.finished <- date()
  class(out) <- "gets"

  if(alarm){ alarm() }
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.gets(out) }
  return(out)
}
