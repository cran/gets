print.arx <-
function(x, ...)
{
  #header:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Method: Ordinary Least Squares (OLS)\n")
  cat("No. of observations (mean eq.):", length(na.trim(x$resids)), "\n")
  cat("No. of observations (variance eq.):", length(na.trim(x$resids.std)), "\n")
  cat("Sample (mean eq.):",
    as.character(index(na.trim(x$resids))[1]), "to",
    as.character(index(na.trim(x$resids))[length(na.trim(x$resids))]), "\n")
  if(!is.null(x$mean.results)){
    cat("\n")
    cat("Mean equation:\n")
    cat("\n")
    print(x$mean.results)
  }
  cat("\n")
  cat("Log-variance equation:\n")
  cat("\n")
  print(x$variance.results)
  cat("\n")
  cat("Diagnostics:\n")
  cat("\n")
  print(x$diagnostics)
}
