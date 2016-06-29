print.isat <-
function(x, ...)
{
  ##determine type:
  spec <- "mean"

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
  if(is.regular(tmp, strict=TRUE)){
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

#OLD:
#  cat("Sample (mean eq.):",
#    as.character(x$aux$y.index[1]), "to",
#    as.character(x$aux$y.index[x$aux$y.n]), "\n")

  ##gum:
  if(spec=="mean"){
    cat("\n")
    cat("GUM mean equation:\n")
    cat("\n")
    printCoefmat(x$gum.mean, dig.tst=0, tst.ind=c(1,2))
#OLD:    print(x$gum.mean)
  }
  if(!is.null(x$gum.variance)){
    cat("\n")
    cat("GUM log-variance equation:\n")
    cat("\n")
    printCoefmat(x$gum.variance)
#OLD:    print(x$gum.variance)
  }
  cat("\n")
  cat("Diagnostics:\n")
  cat("\n")
  printCoefmat(x$gum.diagnostics, dig.tst=0, tst.ind=2)
#OLD:  print(x$gum.diagnostics)

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
    printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4))
    #OLD: print(x$terminals.results)
  }


  ##specific model:
  if(spec=="mean" && !is.null(x$specific.spec)){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if(!is.null(x$mean.results)){
      printCoefmat(x$mean.results)
#OLD:      print(x$mean.results)
    }
    if(x$specific.spec[1]==0){
      cat("empty\n")
    }
##in the future: use estimate.specific=FALSE more directly?
    if(x$specific.spec[1]!=0 && is.null(x$mean.results)){
      cat("Not estimated\n")
    }
  }
  if(!is.null(x$specific.diagnostics)){
    cat("\n")
    cat("Diagnostics:\n")
    cat("\n")
    printCoefmat(x$specific.diagnostics, dig.tst=0, tst.ind=2)
#OLD:    print(x$specific.diagnostics)
  }

  ##delete this one? or at least change?:
  ##notes:
#  if(!is.null(x$notes)){
#    cat("\n")
#    cat("Notes:\n")
#    cat("\n")
#    for(i in 1:length(x$notes)){
#      cat("-",x$notes[[i]],"\n")
#    }
#    cat("\n")
#  }
}
