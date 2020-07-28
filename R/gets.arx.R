gets.arx <-
function(x, spec=NULL, ...)
{
  ##determine spec:
  if(is.null(spec)){
    if( !is.null(x$mean.results) ){ spec <- "mean" }
    if( is.null(x$mean.results)
      && !is.null(x$variance.results) ){ spec <- "variance" }
  }else{
    specType <- c("mean", "variance")
    whichType <- charmatch(spec, specType)
    spec <- specType[ whichType ]  
  }
  
  ##do the gets modelling:
  if( spec=="mean" ){
    result <- getsm(x, ...)
  }else{
    result <- getsv(x, ...)
  }

  ##return result:
  return(result)
    
}
