print.arx <-
function(x, signif.stars=TRUE, ...)
#print.arx <- function(x, signif.stars=FALSE, ...)
{
  ##check if mean and variance have been fitted:
  xNames <- names(x)
  meanResults <- ifelse("mean.results" %in% xNames, TRUE, FALSE)
  varianceResults <- ifelse("variance.results" %in% xNames, TRUE, FALSE)

  ##header - first part:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(meanResults || varianceResults){
    estType <- ifelse(is.null(x$aux$user.estimator),
      "Ordinary Least Squares (OLS)", "User defined")
    cat("Dependent var.:", x$aux$y.name, "\n")
    cat("Method:", estType, "\n")
  }

  ##header - if mean results:
  if(meanResults){
    if(is.null(x$aux$user.estimator)){
      cat("Variance-Covariance:", switch(x$aux$vcov.type,
        ordinary = "Ordinary", white = "White (1980)",
        "newey-west" = "Newey and West (1987)"), "\n")
    }
    if("residuals" %in% xNames){
      cat("No. of observations (mean eq.):",
        length(na.trim(x$residuals)), "\n")
    }
  }

  ##header - if variance results:
  if( varianceResults && "resids.std" %in% xNames ){
    cat("No. of observations (variance eq.):",
      length(na.trim(x$std.residuals)), "\n")
  }

  ##header - sample info:
  if( "residuals" %in% xNames ){
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
  } #end if( "residuals" %in% xNames )

  ##print mean results:
  if(meanResults){
    cat("\n")
    cat("Mean equation:\n")
    cat("\n")
    printCoefmat(x$mean.results, signif.stars=signif.stars)
  }

  ##print variance results:
  if(varianceResults){
    cat("\n")
    cat("Log-variance equation:\n")
    cat("\n")
    printCoefmat(x$variance.results, signif.stars=signif.stars)
  }

  ##print if no results:
  if( !meanResults && !varianceResults ){
    cat("\n")
    cat("   No estimation results\n")
  }

  ##create goodness-of-fit matrix:
  if( !"gof" %in% xNames && is.null(x$aux$user.estimator) ){
    gof <- matrix(NA, 3, 1)
    rownames(gof) <- c("SE of regression", "R-squared",
      paste("Log-lik.(n=", length(na.trim(x$std.residuals)), ")", sep=""))
    colnames(gof) <- ""
    gof[1,1] <- sigma.arx(x)
    gof[2,1] <- rsquared(x)
    gof[3,1] <- as.numeric(logLik.arx(x))
    x$gof <- gof
  }

  ##print diagnostics and fit:
  if( !is.null(x$diagnostics) ) {
    cat("\n")
    cat("Diagnostics and fit:\n")
    cat("\n")
    printCoefmat(x$diagnostics, dig.tst=0, tst.ind=2,
      signif.stars=FALSE)
#NEW (suggested by Moritz)?:
#    printCoefmat(x$diagnostics, tst.ind=2,
#      signif.stars=signif.stars, has.Pvalue=TRUE)
    if( !is.null(x$gof) ){
      printCoefmat(x$gof, digits=6, signif.stars=FALSE)
    }
  }

}
