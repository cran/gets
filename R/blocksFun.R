blocksFun <-
function(y, x, untransformed.residuals=NULL,
  blocks=NULL, no.of.blocks=NULL, max.block.size=30,
  ratio.threshold=0.8, gets.of.union=TRUE, force.invertibility=FALSE,
  user.estimator=list(name="ols"), t.pval=0.001, wald.pval=t.pval,
  do.pet=FALSE, ar.LjungB=NULL, arch.LjungB=NULL,
  normality.JarqueB=NULL, user.diagnostics=NULL,
  gof.function=list(name="infocrit"), gof.method=c("min","max"),
  keep=NULL, include.gum=FALSE, include.1cut=FALSE,
  include.empty=FALSE, max.paths=NULL, turbo=FALSE,
  parallel.options=NULL, tol=1e-07, LAPACK=FALSE, max.regs=NULL,
  print.searchinfo=TRUE, alarm=FALSE)
{
  ## contents:
  ## 1 initiate
  ## 2 x and blocks arguments
  ## 3 loop on x matrices
  ## 4 make return object
    
  ##-------------------------------
  ## 1 initiate
  ##-------------------------------

  gof.method <- match.arg(gof.method)

  ##make result list, add to list:
  result <- list()
  result$call <- sys.call()
  result$time.started <- date()
  result$time.finished <- NA
  result$messages <- NULL

  ##parallel.options argument:
  if(!is.null(parallel.options)){

    ##if(numeric):
    if(is.numeric(parallel.options)){
      clusterSpec <- parallel.options
      OScores <- detectCores()
      if(parallel.options > OScores){
        stop("parallel.options > number of cores/threads")
      }
    }

    ##varlist for clusterExport:
    if(is.list(parallel.options)){
      clusterVarlist <- parallel.options$varlist
    }else{
      clusterVarlist <- NULL
    }
    clusterVarlist <- c(clusterVarlist,
      "dropvar", "getsFun", "ols", "infocrit", "diagnostics")
    if(!is.null(user.diagnostics)){
      clusterVarlist <- c(clusterVarlist, user.diagnostics$name)
    }
    if(!is.null(user.estimator)){
      clusterVarlist <- c(clusterVarlist, user.estimator$name)
    }
    if(!is.null(gof.function)){
      clusterVarlist <- c(clusterVarlist, gof.function$name)
    }

    #for the future?: add memory.limit()/memory.size() = max cores check?

  } #end if(!is.null(parallel.options))


  ##-------------------------------
  ## 2 x argument
  ##-------------------------------
  
  ##x is a matrix:
  if( is.matrix(x) ){
    xMatrixName <- deparse(substitute(x))
    x <- list(x=x) #convert to list
    names(x) <- xMatrixName
  }

  ##x is a list of matrices:
  if( is.list(x) ){

    ##ensure matrices are named:
    xMatrixNames <- paste0("X", 1:length(x))
    if( is.null(names(x)) ){
      names(x) <- xMatrixNames
    }else{
      for(i in 1:length(x)){
        if( names(x)[i] %in% c("", NA) ){
          names(x)[i] <- xMatrixNames[i]
        }
      } #close for..loop
    }

    ##handle colnames:
    for(i in 1:length(x)){

      xColNames <- colnames(x[[i]])
      if( is.null(xColNames) ){
        xColNames <- paste0("X", i, ".xreg", 1:NCOL(x[[i]]))
      }
      if( any(xColNames=="") ){
        missing.colnames <- which(xColNames == "")
        for(j in 1:length(missing.colnames)){
          #fixed by Jonas: 
          xColNames[ missing.colnames[j] ] <-
            paste0("X", i, ".xreg", missing.colnames[j])
        }
      }
      xColNames <- make.unique(xColNames)
      colnames(x[[i]]) <- xColNames

    } #end for(..) loop
    
    ##do NOT check that colnames() are unique across matrices,
    ##so that matrices can contain the same regressors
    
  } #end if( is.list(x) )

  ##-------------------------------
  ## 3 blocks and keep arguments
  ##-------------------------------

  ##check blocks:
  if( is.list(blocks) ){
    if( length(x)!=length(blocks) ){
      stop("No. of matrices unequal to length(blocks)")
    }
    blocks.is.list <- TRUE
  }else{
    blocks.is.list <- FALSE
    blocks <- list()
  }

  ##keep is vector:
  if( !is.null(keep) && !is.list(keep) && is.vector(keep) ){
    keeptmp <- keep
    keep <- list()
    keep[[1]] <- keeptmp
    if( length(x)>1 ){
      for(i in 2:length(x)){ keep[[i]] <- integer(0) }
    }
  }

  ##check keep argument, name keep items:
  if( is.list(keep) ){
  
    ##check keep argument:
    if( length(x)!=length(keep) ){
      stop("Length(keep) unequal to no. of matrices in 'x'")
    }
    
    ##name keep items:
    for(i in 1:length(keep)){
      if( length(keep[[i]])>0 ){
        names(keep[[i]]) <- colnames(x[[i]])[ keep[[i]] ]
      } 
    }

    ##name entries in keep list:
    names(keep) <- names(x)
        
  } #end if( is.list(keep) )
 
  
  ##-------------------------------
  ## 4 loop on x matrices
  ##-------------------------------

  ##create list w/union of retained regressors from
  ##each x matrix:
  xUnionOfModels <- list() 
  
  ##loop on x:
  for(i in 1:length(x)){

    ##add entry i to list:
    xUnionOfModels[[i]] <- integer(0)
    
    ##blocks:
    if( !blocks.is.list ){

      y.n <- NROW(y)
      ncol.adj <- NCOL(x[[i]])

      ##determine no. of blocks:
      if( is.null(no.of.blocks) ){
        blockratio.value <- ncol.adj/(ratio.threshold*ncol.adj)
        blocksize.value <-
          ncol.adj/min(y.n*ratio.threshold, max.block.size)
        no.of.blocks <- max(2,blockratio.value,blocksize.value)
        no.of.blocks <- ceiling(no.of.blocks)
        no.of.blocks <- min(ncol.adj, no.of.blocks) #ensure blocks < NCOL
      }

      ##make partitions:
      blocksize <- ceiling(ncol.adj/no.of.blocks)
      partitions.t2 <- blocksize
      for(j in 1:no.of.blocks){
        if( blocksize*j <= ncol.adj ){
          partitions.t2[j] <- blocksize*j
        }
      }      
      ##check if last block contains last regressor:
      if( partitions.t2[ length(partitions.t2) ] < ncol.adj ){
        partitions.t2 <- c(partitions.t2, ncol.adj)
      }
      blocksadj <- length(partitions.t2)
      partitions.t1 <- partitions.t2 + 1
      partitions.t1 <- c(1, partitions.t1[ -blocksadj ])

      ##finalise:
      tmp <- list()
      for(j in 1:blocksadj){
        tmp[[j]] <- partitions.t1[j]:partitions.t2[j]
      }
      blocks[[i]] <- tmp

    } #end if(!blocks.is.list)

#    ##name the blocks:
#    names(blocks[[i]]) <- paste0("block", 1:length(blocks[[i]]))

    ##add keep entries to blocks:
    if( length(keep[[i]])>0 ){
      for(j in 1:length(blocks[[i]])){
        blocks[[i]][[j]] <- sort(union(keep[[i]], blocks[[i]][[j]]))
      }
    }
            
    ##xkeep argument:
    if( length(keep[[i]])==0 ){
      xkeep <- NULL
    }else{ xkeep <- names(keep[[i]]) }
    
    ##make blocks function for lapply/parLapply:
    XblocksFun <- function(j, i, x, blocks,
      parallel.options, y, untransformed.residuals, user.estimator,
      t.pval, wald.pval, do.pet, ar.LjungB, arch.LjungB,
      normality.JarqueB, user.diagnostics, gof.function, gof.method,
      xkeep, include.gum, include.1cut, include.empty, max.paths,
      turbo, force.invertibility, tol, LAPACK, max.regs,
      print.searchinfo){

      ##check if block contains 1 regressor:
      if( length(blocks[[i]][[j]])==1 ){
        tmp <- colnames(x[[i]])[ blocks[[i]][[j]] ]
        mX <- cbind(x[[i]][, blocks[[i]][[j]] ])
        colnames(mX) <- tmp
      }else{
        mX <- cbind(x[[i]][, blocks[[i]][[j]] ])
      }

      ##apply dropvar:
      if( force.invertibility ){
        mX <- dropvar(mX, tol=tol, LAPACK=LAPACK, silent=TRUE)
      }

      ##set xkeep argument:
      if( !is.null(xkeep) ){
        xkeep <- which( colnames(mX) %in% xkeep )
      }
      
      ##print info:
      if( is.null(parallel.options) ){
        if(print.searchinfo){
          message("\n", appendLF=FALSE)
          message(names(x)[i],
            " block ", j, " of ", length(blocks[[i]]), ":",
            appendLF=TRUE)
        }
      }

      ##do gets inside XblocksFun:
      getsx <- getsFun(y, mX,
        untransformed.residuals=untransformed.residuals,
        user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, gof.function=gof.function,
        gof.method=gof.method, keep=xkeep, include.gum=include.gum,
        include.1cut=include.1cut, include.empty=include.empty,
        max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
        max.regs=max.regs, print.searchinfo=print.searchinfo,
        alarm=FALSE)
      if( !is.null(getsx$messages) ){ message(getsx$messages) }

      if( is.null(getsx$specific.spec) ){
        xSpecificmodels <- NULL
      }else{
        xSpecificmodels <- names(getsx$specific.spec)
      }

      ##return
      return(xSpecificmodels)
      ##for the future?: return(getsx) - careful! this would change
      ##the subsequent code substantially

    } #close XblocksFun

    ##call XblocksFun/do gets on each block: NO parallel computing
    if( is.null(parallel.options) ){
      xSpecificmodels <- lapply(1:length(blocks[[i]]),
        XblocksFun, i, x, blocks, parallel.options,
        y, untransformed.residuals, user.estimator, t.pval, wald.pval,
        do.pet, ar.LjungB, arch.LjungB, normality.JarqueB,
        user.diagnostics, gof.function, gof.method, xkeep, include.gum,
        include.1cut, include.empty, max.paths, turbo,
        force.invertibility, tol, LAPACK, max.regs, print.searchinfo)
    }
      
    ##call XblocksFun/do gets on each block: WITH parallel computing
    if( !is.null(parallel.options) ){

      ##print info:
      if(print.searchinfo){
        message("\n", appendLF=FALSE)
        message("Preparing parallel computing...",
          appendLF=TRUE)
        message(names(x)[i],
          " blocks to search in parallel: ", length(blocks[[i]]),
          appendLF=TRUE)
        message("Searching...", appendLF=TRUE)
      }

      blocksClust <- makeCluster(clusterSpec, outfile="") #make cluster
      clusterExport(blocksClust, clusterVarlist,
        envir=.GlobalEnv) #idea for the future?: envir=clusterEnvir
      xSpecificmodels <- parLapply(blocksClust,
        1:length(blocks[[i]]), XblocksFun, i, x,
        blocks, parallel.options, y, untransformed.residuals,
        user.estimator, t.pval, wald.pval, do.pet, ar.LjungB,
        arch.LjungB, normality.JarqueB, user.diagnostics, gof.function,
        gof.method, xkeep, include.gum, include.1cut, include.empty,
        max.paths, turbo, force.invertibility, tol, LAPACK, max.regs,
        print.searchinfo)
      stopCluster(blocksClust)

    } #end if( parallel computing )

    ##union of retained variables:
    ##------------------------------------

    ##union of retained variables (names):
    xNames <- NULL
    if( length(xSpecificmodels)>0 ){
      #which variables retained?:
      for(j in 1:length(xSpecificmodels)){
        #check if non-empty:
        if( !is.null(xSpecificmodels[[j]]) ){
          xNames <- union(xNames, xSpecificmodels[[j]])
        }
      }
    }

    ##NOT do gets of union:
    if( gets.of.union==FALSE ){ xUnionOfModels[[i]] <- xNames }
    
    ##DO gets of union:
    if( gets.of.union==TRUE && length(xNames)>0 ){

      if( print.searchinfo ){
        message("\n", appendLF=FALSE)
        message("GETS of union of retained ",
          names(x)[i], " variables... ",
          appendLF=TRUE)
        message("\n", appendLF=FALSE)
      }
  
      ##build regressor matrix:
      mX <- cbind(x[[i]][,xNames])
      colnames(mX) <- xNames
      if( force.invertibility ){
        mX <- dropvar(mX, tol=tol, LAPACK=LAPACK, silent=TRUE)
      }
      
      ##build xkeep:
      if( !is.null(keep) ){
        xkeep <- NULL
        for(j in 1:length(keep)){
          xkeep <- union(xkeep, names(keep[[j]]))
        }
        xkeep <- which( colnames(mX) %in% xkeep)
      }
      
      ##do gets of union:
      getsx <- getsFun(y, mX,
        untransformed.residuals=untransformed.residuals,
        user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, gof.function=gof.function,
        gof.method=gof.method, keep=xkeep, include.gum=include.gum,
        include.1cut=include.1cut, include.empty=include.empty,
        max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
        max.regs=max.regs, print.searchinfo=print.searchinfo,
        alarm=FALSE)
      if( !is.null(names(getsx$specific.spec)) ){
        xUnionOfModels[[i]] <- names(getsx$specific.spec)
      }

    } #end if( do gets )

  } #end for(i) loop (on x matrices)

  ##add names:
  names(blocks) <- names(x)
  names(xUnionOfModels) <- names(x)


  ##-------------------------------
  ## 5 make return object:
  ##-------------------------------

  result$y <- y
  result$x <- list()
  for(i in 1:length(x)){ result$x[[i]] <- colnames(x[[i]]) }
  names(result$x) <- names(x)  
  result$blocks <- blocks
  result$keep <- keep
  result$specific.spec <- xUnionOfModels
  result$time.finished <- date()
  if(alarm){ alarm() }
  return(result)
  
}
