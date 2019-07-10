printtex <-
function(x, fitted.name=NULL, xreg.names=NULL,
  digits=3, intercept=TRUE, gof=TRUE, diagnostics=TRUE)
{

  ##is class(x)="arx"/"gets"/"isat"?:
  ##---------------------------------
  
  xName <- deparse(substitute(x))
  xClass <- class(x)
  if( xClass %in% c("arx","gets","isat") ){
    yName <- ifelse(is.null(fitted.name), x$aux$y.name, fitted.name)
  }else{
    yName <- ifelse(is.null(fitted.name), "y", fitted.name)
    message(paste0("\n '", xName, "'", " is not of class 'arx', ",
      "'gets' or 'isat', LaTeX code may contain errors:\n"))
  }
  yName <- paste0("\\widehat{", yName, "}")

  ##equation:
  ##---------

  coefs <- coef(x)
  if(is.null(xreg.names)){
    coefsNames <- names(coefs)
  }else{
    coefsNames <- xreg.names
    if( length(coefs) != length(xreg.names) ){
      message(paste0("\n length of 'xreg.names' does not match",
        " length of 'coef(x)'\n"))
    }
  }
  intercept <- as.numeric(intercept)
  if( intercept > 0 ){ coefsNames[ intercept ] <- "" }
  coefs <- as.numeric(coefs)
  stderrs <- as.numeric(sqrt(diag(vcov(x))))

  eqtxt <- NULL
  if(length(coefs) > 0){
    for(i in 1:length(coefs) ){
      ifpluss <- ifelse(i==1, "", " + ")
      eqtxt <- paste(eqtxt,
        ifelse(coefs[i] < 0, " - ", ifpluss), "\\underset{(",
        format(round(stderrs[i], digits=digits), nsmall=digits), ")}{",
        format(round(abs(coefs[i]), digits=digits), nsmall=digits), "}",
        coefsNames[i], sep="")
    }
  }

  txtAddEq <- ifelse(gof+diagnostics>0, " \\\\[2mm]", "")
  eqtxt <- paste0("  ", yName, " &=& ", eqtxt, "", txtAddEq, " \n")

  ##goodness of fit:
  ##----------------

  goftxt <- NULL
  if(gof){
    goftxt <- "   &&"
    iT <- ""
    if(xClass %in% c("arx","gets","isat") ){
      goftxt <- paste(goftxt, " R^2=",
        format(round(rsquared(x), digits=digits), nsmall=digits),
        " \\qquad \\widehat{\\sigma}=",
        format(round(sigma(x), digits=digits), nsmall=digits),
        sep="")
      iT <- x$aux$y.n
    }
    goftxt <- paste(goftxt, " \\qquad LogL=",
      format(round(as.numeric(logLik(x)), digits=digits), nsmall=digits),
        "\\qquad T = ", iT, " \\nonumber \\\\ \n", sep="")
  }
  
  ##diagnostics:
  ##------------

  diagtxt <- NULL
  if(xClass=="arx" && diagnostics==TRUE){
    dfDiags <- diagnostics(x,
      ar.LjungB=c(ar.LjungB=x$aux$qstat.options[1],1),
      arch.LjungB=c(ar.LjungB=x$aux$qstat.options[2],1),
      normality.JarqueB=TRUE, verbose=TRUE)
    diagtxt <- paste("  ", " && \\underset{[p-val]}{ AR(",
      x$aux$qstat.options[1], ") }:", " \\underset{[",
      format(round(dfDiags[1,3], digits=digits), nsmall=digits), "]}{",
      format(round(dfDiags[1,1], digits=digits), nsmall=digits), "}",
      "\\qquad \\underset{[p-val]}{ ARCH(",
      x$aux$qstat.options[2], ")}:", "\\underset{[",
      format(round(dfDiags[2,3], digits=digits), nsmall=digits), "]}{",
      format(round(dfDiags[2,1], digits=digits), nsmall=digits), "}",
      "\\qquad \\underset{[p-val]}{ Normality }:", "\\underset{[",
      format(round(dfDiags[3,3], digits=digits), nsmall=digits), "]}{",
      format(round(dfDiags[3,1], digits=digits), nsmall=digits), "}",
      " \\nonumber \n", sep="")
  }
  
  ##print code:
  ##-----------

  cat("\\begin{eqnarray}\n")
  cat(eqtxt)
  cat(goftxt)
  cat(diagtxt)
  cat("\\end{eqnarray}\n")

}
