coef.gets <-
function(object, spec=NULL, ...)
{
  if( is.null(spec) ){
    spec <- switch(object$gets.type, getsm="mean", getsv="variance")
  }
  coef.arx(object, spec=spec)
}
