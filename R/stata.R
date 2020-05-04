stata <-
function(object, file=NULL, print=TRUE,
  return=FALSE)
{
  out <- list()
  out$object.name <- deparse(substitute(object))

  ##index, data, names:
  out$index <- object$aux$y.index
  out$data <- cbind(object$aux$y, object$aux$mX)
  out$data <- as.data.frame(out$data)
  out$data <- cbind(as.character(out$index), out$data)
  out$names <- gsub("[.]","",tolower(c("index", object$aux$y.name, object$aux$mXnames)))
  where.mconst <- which(out$names=="mconst")
  if(length(where.mconst) > 0){
    out$data <- out$data[-where.mconst]
    out$names <- out$names[-where.mconst]
    noConstant <- FALSE
  }else{
    noConstant <- TRUE
  }
  colnames(out$data) <- out$names

  ##Stata code to estimate the model:
  outNames <- out$names
  outNames[1] <- "regress"
  out$regress <- paste(outNames, collapse=" ")
  if( noConstant==TRUE || object$aux$vcov.type!="ordinary" ){

    cmdOptions <- NULL
    if(noConstant){ cmdOptions <- c(cmdOptions, "noconstant") }
    if(object$aux$vcov.type!="ordinary"){
      cmdOptions <- c(cmdOptions, "vce(robust)")
    }
    cmdOptions <- paste(cmdOptions, collapse=" ")
    out$regress <- paste(out$regress, ",", cmdOptions, collapse="")

  }

  ##if print=TRUE and is.null(file):
  if(print && is.null(file)){

    ##Stata code to estimate the model:
    message("STATA code to estimate the model:\\n")
    message(" ", out$regress, "\\n")

    ##R code to export the data:
    message("R code (example) to export the data of the model:\\n")
    message(paste("  stata(", out$object.name, ", file='C:/Users/myname/Documents/getsdata.csv')\\n", sep=""))

  } #close if(print && is.null(file))

  ##if save data:
  if(!is.null(file)){
    write.csv(out$data, file, row.names=FALSE)
    ##if print=TRUE:
    if(print){
      message("Data saved in:\\n")
      message("  ", file, "\\n", sep="")
      message("STATA code to estimate the model:\\n")
      message(" ", out$regress, "\\n")
    }
  } #end if(!is.null(file))

  ##out:
  if(return){ return(out) }

}
