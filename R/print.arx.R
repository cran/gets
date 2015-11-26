print.arx <-
function(x, ...)
{
  ##check if mean and variance have been fitted:
  if(is.null(x$mean.results)){
    meanResults <- FALSE
  }else{
    meanResults <- TRUE
  }
  if(is.null(x$variance.results)){
    varianceResults <- FALSE
  }else{
    varianceResults <- TRUE
  }

  ##header:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(meanResults || varianceResults){
    cat("Method: Ordinary Least Squares (OLS)\n")
  }
  if(meanResults){
#  if(!is.null(x$mean.results)){
    cat("Variance-Covariance:", switch(x$aux$vcov.type,
      ordinary = "Ordinary", white = "White (1980)",
      "newey-west" = "Newey and West (1987)"), "\n")
    cat("No. of observations (mean eq.):",
      length(na.trim(x$resids)), "\n")
  }
  if(varianceResults){
#  if(!is.null(x$variance.results)){
    cat("No. of observations (variance eq.):",
      length(na.trim(x$resids.std)), "\n")
  }
  cat("Sample:",
    as.character(index(na.trim(x$resids))[1]), "to",
    as.character(index(na.trim(x$resids))[length(na.trim(x$resids))]), "\n")

  ##mean results:
  if(meanResults){
#  if(!is.null(x$mean.results)){
    cat("\n")
    cat("Mean equation:\n")
    cat("\n")
    print(x$mean.results)
  }

  ##variance results:
  if(varianceResults){
#  if(!is.null(x$variance.results)){
    cat("\n")
    cat("Log-variance equation:\n")
    cat("\n")
    print(x$variance.results)
  }

  ##Fit and diagnostics:
  mGOF <- matrix(NA, 2, 1)
  rownames(mGOF) <- c("R-squared",
    paste("Log-lik.(n=", length(na.trim(x$resids.std)), ")", sep=""))
  colnames(mGOF) <- ""
  mGOF[1,1] <- x$diagnostics[4,1]
  mGOF[2,1] <- as.numeric(logLik.arx(x))
  cat("\n")
  cat("Diagnostics:\n")
  cat("\n")
  print(x$diagnostics[1:3,])
  print(mGOF)

}
