print.isat <-
function(x, ...)
{
  ##specification type:
  specType <- "mean"

  ##header:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Dependent var.:", x$aux$y.name, "\n")
  cat("Method: Ordinary Least Squares (OLS)\n")
  cat("Variance-Covariance:", switch(x$aux$vcov.type,
    ordinary = "Ordinary", white = "White (1980)",
    "newey-west" = "Newey and West (1987)"), "\n")

  ##header - sample info:
  cat("No. of observations (mean eq.):", x$aux$y.n, "\n")
  tmp <- zoo(x$aux$y, order.by=x$aux$y.index)

  indexTrimmed <- index(na.trim(tmp))
  isRegular <- is.regular(tmp, strict=TRUE)
  isCyclical <- frequency(tmp) > 1
  if(isRegular && isCyclical){
    cycleTrimmed <- cycle(na.trim(tmp))
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
    printCoefmat(x$gum.variance, signif.stars=FALSE)
  }
  cat("\n")
  cat("Diagnostics and fit:\n")
  cat("\n")
  printCoefmat(x$gum.diagnostics, dig.tst=0, tst.ind=2,
    signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)

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
  if(!is.null(x$terminals)){
    cat("\n")
    cat("Terminal models: \n")
    cat("\n")
    for(i in 1:length(x$terminals)){
      cat("spec",i,":",x$terminals[[i]],"\n")
    }
  }
  if(!is.null(x$terminals.results)){
    cat("\n")
    printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4),
      signif.stars=FALSE)
  }

  ##specific model:
  if(specType=="mean" && !is.null(x$specific.spec)){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if(!is.null(x$mean.results)){
      print(x$mean.results)
#OLD: DOES NOT WORK IN A PREDICTABLE WAY!
#      printCoefmat(x$mean.results, signif.stars=FALSE,
#        P.values=FALSE, has.Pvalues=FALSE)
    }
    if(x$specific.spec[1]==0){
      cat("empty\n")
    }
##in the future: use estimate.specific=FALSE more directly?
    if(x$specific.spec[1]!=0 && is.null(x$mean.results)){
      cat("Not estimated\n")
    }
  }
  
  ##diagnostics and fit:
  if(!is.null(x$specific.diagnostics)){

    #fit-measures:
    mGOF <- matrix(NA, 3, 1)
    rownames(mGOF) <- c("SE of regression", "R-squared",
      paste("Log-lik.(n=", length(na.trim(x$std.residuals)), ")", sep=""))
    colnames(mGOF) <- ""
    mGOF[1,1] <- sigma.isat(x) #OLD: sqrt( RSS/(nobs-DFs) )
    mGOF[2,1] <- rsquared(x) #OLD: x$specific.diagnostics[4,1]
    mGOF[3,1] <- as.numeric(logLik.arx(x))

    cat("\n")
    cat("Diagnostics and fit:\n")
    cat("\n")
    printCoefmat(x$specific.diagnostics, dig.tst=0, tst.ind=2,
      signif.stars=FALSE)
    printCoefmat(mGOF, digits=6, signif.stars=FALSE)

  }

  ##messages:
  if(!is.null(x$messages)){
    message("\n", appendLF=FALSE)
    message(x$messages)
  }

}
