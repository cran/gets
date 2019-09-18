leqwma <-
function(x, length=5, k=1, p=2, as.vector=FALSE,
  lag=NULL, start=NULL)
{
  eqwma(x, length=length, k=k, p=p, log=TRUE, abs=TRUE,
    as.vector=as.vector, lag=lag, start=start)
}
