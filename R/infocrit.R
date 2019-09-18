infocrit <-
function(x, method=c("sc", "aic", "aicc", "hq"))
{
  method <- match.arg(method)
  ##Schwarch criterion:
  if(method == "sc") infoVal <- -2*x$logl/x$n + x$k*log(x$n)/x$n
  ##Akaike criterion:
  if(method == "aic") infoVal <- -2*x$logl/x$n + 2*x$k/x$n
  ##AICC of Hurvich and Tsai (1989):
  if(method == "aicc") infoVal <- -2*x$logl/x$n + 2*x$k/(x$n-x$k-1)
  ##Hannan-Quinn criterion:
  if(method == "hq") infoVal <- -2*x$logl/x$n + 2*x$k*log(log(x$n))/x$n
  ##return:
  return(infoVal)
}
