print.gets <-
function(x, ...)
{
  ##determine spec:
  specType <- switch(as.character(x$call)[1],
    getsm="mean", getsv="variance")

  ##header - first part:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(specType=="mean"){
    cat("Dependent var.:", x$aux$y.name, "\n")
  }
  cat("Method: Ordinary Least Squares (OLS)\n")

  ##header - if mean:
  if(specType=="mean"){
    cat("Variance-Covariance:", switch(x$aux$vcov.type,
      ordinary = "Ordinary", white = "White (1980)",
      "newey-west" = "Newey and West (1987)"), "\n")
    if(!is.null(x$aux$y.n)){
      cat("No. of observations (mean eq.):", x$aux$y.n, "\n") }
  }

  ##header - if variance:
  if(specType=="variance"){
    if(!is.null(x$aux$loge2.n)){
      cat("No. of observations (variance eq.):",
        x$aux$loge2.n, "\n") }
  }

  ##header - sample info:
  if(!is.null(x$resids)){
    indexTrimmed <- index(na.trim(x$resids))
    isRegular <- is.regular(x$resids, strict=TRUE)
    isCyclical <- frequency(x$resids) > 1
    if(isRegular && isCyclical){
      cycleTrimmed <- cycle(na.trim(x$resids))
      startYear <- floor(as.numeric(indexTrimmed[1]))
      startAsChar <- paste(startYear,
        "(", cycleTrimmed[1], ")", sep="")
      endYear <- floor(as.numeric(indexTrimmed[length(indexTrimmed)]))
      endAsChar <- paste(endYear,
        "(", cycleTrimmed[length(indexTrimmed)], ")", sep="")
    }else{
      startAsChar <- as.character(indexTrimmed[1])
      endAsChar <- as.character(indexTrimmed[length(indexTrimmed)])
    }
    cat("Sample:", startAsChar, "to", endAsChar, "\n")
  } #end if(!is.null..)

  ##gum:
  if(specType=="mean"){
    cat("\n")
    cat("GUM mean equation:\n")
    cat("\n")
    printCoefmat(x$gum.mean, dig.tst=0, tst.ind=c(1,2),
      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }
  if(!is.null(x$gum.variance)){
    cat("\n")
    cat("GUM log-variance equation:\n")
    cat("\n")
    if(specType=="mean"){
      printCoefmat(x$gum.variance, signif.stars=FALSE)
    }
    if(specType=="variance"){
      printCoefmat(x$gum.variance, dig.tst=0, tst.ind=c(1,2),
      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
    }
  }
  if(!is.null(x$gum.diagnostics)){
    cat("\n")
    cat("Diagnostics:\n")
    cat("\n")
    printCoefmat(x$gum.diagnostics, dig.tst=0, tst.ind=2,
      signif.stars=FALSE)
  }

  ##paths:
  cat("\n")
  cat("Paths searched: \n")
  cat("\n")
  if(is.null(x$paths)){
    print(NULL)
  }else{
    for(i in 1:length(x$paths)){
      cat("path",i,":",x$paths[[i]],"\n")
    }
  } #end if(is.null(x$paths))

  ##terminal models and results:
  cat("\n")
  cat("Terminal models: \n")
  if(!is.null(x$terminals)){
    cat("\n")
    for(i in 1:length(x$terminals)){
      cat("spec",i,":",x$terminals[[i]],"\n")
    }
  }
  cat("\n")
  printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4),
    signif.stars=FALSE)

  ##specific model:
  if(specType=="mean" && !is.null(x$specific.spec)){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if(!is.null(x$mean.results)){
      printCoefmat(x$mean.results, signif.stars=FALSE)
    }
    if(x$specific.spec[1]==0){
      cat("empty\n")
    }
##in the future: use estimate.specific=FALSE more directly?
    if(x$specific.spec[1]!=0 && is.null(x$mean.results)){
      cat("Not estimated, since estimate.specific=FALSE\n")
    }
  }
  if(!is.null(x$variance.results)){
    cat("\n")
    cat("SPECIFIC log-variance equation:\n")
    cat("\n")
    printCoefmat(x$variance.results, signif.stars=FALSE)
#    printCoefmat(x$variance.results, dig.tst=0, tst.ind=c(1,2),
#      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }

  ##diagnostics and fit:
  if(!is.null(x$specific.diagnostics)){

    #fit-measures:
    mGOF <- matrix(NA, 3, 1)
    rownames(mGOF) <- c("SE of regression", "R-squared",
      paste("Log-lik.(n=", length(na.trim(x$resids.std)), ")", sep=""))
    colnames(mGOF) <- ""
    mGOF[1,1] <- sigma.gets(x) #OLD: sqrt( RSS/(nobs-DFs) )
    mGOF[2,1] <- rsquared(x) #OLD: x$specific.diagnostics[4,1]
    mGOF[3,1] <- as.numeric(logLik.arx(x))

    cat("\n")
    cat("Diagnostics:\n")
    cat("\n")
    printCoefmat(x$specific.diagnostics, dig.tst=0, tst.ind=2,
      signif.stars=FALSE)
    printCoefmat(mGOF, digits=6, signif.stars=FALSE)
  #OLD: print(mGOF)

  }

  ##messages:
  if(!is.null(x$messages)){
    message("\n", appendLF=FALSE)
    message(x$messages)
  }

}
