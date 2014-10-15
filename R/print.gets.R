print.gets <-
function(x, ...)
{
  #determine type:
  if(as.character(x$call)[1]=="getsm"){ spec <- "mean" }
  if(as.character(x$call)[1]=="getsv"){ spec <- "variance" }

  #header:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Method: Ordinary Least Squares (OLS)\n")
  if(spec=="mean"){
    cat("No. of observations (mean eq.):", length(na.trim(x$resids)), "\n")
    cat("Sample (mean eq.):",
      as.character(index(na.trim(x$resids))[1]), "to",
      as.character(index(na.trim(x$resids))[length(na.trim(x$resids))]), "\n")
  }else{
    cat("No. of observations (variance eq.):", length(na.trim(x$gum.resids.std)), "\n")
    cat("Sample (variance eq.):",
      as.character(index(na.trim(x$gum.resids.std))[1]), "to",
      as.character(index(na.trim(x$gum.resids.std))[length(na.trim(x$gum.resids.std))]), "\n")
  }

  #gum:
  if(spec=="mean"){
    cat("\n")
    cat("GUM mean equation:\n")
    cat("\n")
    print(x$gum.mean)
  }
  cat("\n")
  cat("GUM log-variance equation:\n")
  cat("\n")
  print(x$gum.variance)
  cat("\n")
  cat("Diagnostics:\n")
  cat("\n")
  print(x$gum.diagnostics)

  #paths:
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

  #terminal models and results:
  cat("\n")
  cat("Terminal models: \n")
  if(!is.null(x$terminals)){
    cat("\n")
    for(i in 1:length(x$terminals)){
      cat("spec",i,":",x$terminals[[i]],"\n")
    }
  }
  cat("\n")
  print(x$terminals.results)

  #specific model:
  if(spec=="mean"){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    print(x$specific.mean)
  }
  cat("\n")
  cat("SPECIFIC log-variance equation:\n")
  cat("\n")
  print(x$specific.variance)
  cat("\n")
  cat("Diagnostics:\n")
  cat("\n")
  print(x$specific.diagnostics)

  #notes:
  if(!is.null(x$notes)){
    cat("\n")
    cat("Notes:\n")
    cat("\n")
    for(i in 1:length(x$notes)){
      cat("-",x$notes[[i]],"\n")
    }
    cat("\n")
  }
}
