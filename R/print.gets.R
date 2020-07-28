print.gets <-
function(x, signif.stars=TRUE, ...)
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
  estType <- ifelse(is.null(x$aux$user.estimator),
      "Ordinary Least Squares (OLS)", "User defined")
  cat("Method:", estType, "\n")

  ##header - if mean:
  if( specType=="mean" ){
    vcovType <- "Unknown"
    if( !is.null(x$aux$vcov.type) ){
      vcovType <- switch(x$aux$vcov.type,
        ordinary = "Ordinary", white = "White (1980)",
        "newey-west" = "Newey and West (1987)")
    }
    cat("Variance-Covariance:", vcovType, "\n")
    if(!is.null(x$aux$y.n)){
      cat("No. of observations (mean eq.):", x$aux$y.n, "\n") }
  }

  ##header - if variance:
  if( specType=="variance" ){
    if(!is.null(x$aux$loge2.n)){
      cat("No. of observations (variance eq.):",
        x$aux$loge2.n, "\n") }
  }

  ##header - sample info:
  if( !is.null(x$residuals) ){
    indexTrimmed <- index(na.trim(x$residuals))
    isRegular <- is.regular(x$residuals, strict=TRUE)
    isCyclical <- frequency(x$residuals) > 1
    if(isRegular && isCyclical){
      cycleTrimmed <- cycle(na.trim(x$residuals))
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
  if( specType=="mean" && !is.null(x$gum.mean) ){
    cat("\n")
    cat("GUM mean equation:\n")
    cat("\n")
    printCoefmat(x$gum.mean, tst.ind=c(1,2),
      signif.stars=signif.stars)
#OLD:
#    printCoefmat(x$gum.mean, dig.tst=0, tst.ind=c(1,2),
#      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }
  if( !is.null(x$gum.variance) ){
    cat("\n")
    cat("GUM log-variance equation:\n")
    cat("\n")
    if(specType=="mean"){
      printCoefmat(x$gum.variance, signif.stars=FALSE)
    }
    if(specType=="variance"){
      printCoefmat(x$gum.variance, tst.ind=c(1,2),
        signif.stars=signif.stars)
#OLD:
#      printCoefmat(x$gum.variance, dig.tst=0, tst.ind=c(1,2),
#        signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
    }
  }
  if( !is.null(x$gum.diagnostics) ){
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
  if(!is.null(x$terminals)){
    cat("\n")
    cat("Terminal models: \n")
    if(!is.null(x$terminals)){
      cat("\n")
      for(i in 1:length(x$terminals)){
        cat("spec",i,":",x$terminals[[i]],"\n")
      }
    }
  }
  if(!is.null(x$terminals.results)){
    cat("\n")
    printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4),
      signif.stars=FALSE)
  }
  
  ##specific mean model:
  if( specType=="mean" && !is.null(x$terminals.results) ){
#OLD: if( specType=="mean" && !is.null(x$specific.spec) ){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if( !is.null(x$mean.results) ){
      printCoefmat(x$mean.results, signif.stars=signif.stars)
#      printCoefmat(x$mean.results, signif.stars=FALSE)
    }
    if( length(x$specific.spec)==0 ){
#OLD: if( x$specific.spec[1]==0 ){
      cat("  the empty model\n")
    }
  }

  ##specific log-variance model:
  if( !is.null(x$variance.results) ){
    cat("\n")
    cat("SPECIFIC log-variance equation:\n")
    cat("\n")
    printCoefmat(x$variance.results, signif.stars=signif.stars)
#    printCoefmat(x$variance.results, signif.stars=FALSE)
#    printCoefmat(x$variance.results, dig.tst=0, tst.ind=c(1,2),
#      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }

  ##diagnostics and fit:
  if( !is.null(x$specific.diagnostics) ){

    #fit-measures:
    mGOFnames <- "SE of regression"
    mGOF <- sigma.gets(x) 
    if( specType == "mean" ){
      mGOFnames <- c(mGOFnames, "R-squared")
      mGOF <- rbind(mGOF, rsquared(x))
    }
    logl <- logLik.arx(x)
    mGOFnames <- c(mGOFnames,
      paste0("Log-lik.(n=", attr(logl,"n"), ")") )
    mGOF <- rbind(mGOF, as.numeric(logl))
    rownames(mGOF) <- mGOFnames
    colnames(mGOF) <- ""

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
