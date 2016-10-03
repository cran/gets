predict.isat <-
function(object, n.ahead=12,
  newmxreg=NULL, newindex=NULL, return=TRUE, plot=TRUE,
  plot.options=list(), ...)
{

  ##create new object to add stuff to in order to use predict.arx()
  object.new <- object

  ##------------
  ## arguments:
  ##------------

  if("mX" %in% names(object$aux)) {

    ##check if constant is retained:
    if("mconst" %in% object$aux$mXnames){
      object.new$call$mc <- TRUE
    }else{
      object.new$call$mc <- NULL
    }

    ##what dynamics specified in gum?
    gum.ar <- eval(object$call$ar)
    ##what dynamics remain in spec?
    spec.ar <- as.numeric(gsub("ar(\\d+)","\\1",object$aux$mXnames[grep("^ar\\d+$",object$aux$mXnames)]))
    if(NROW(spec.ar)==0) {
      object.new$call$ar <- NULL
    } else { ##check that dynamics in specific are subset of those in gum
      object.new$call$ar <- spec.ar[spec.ar %in% gum.ar]
    }

    ##"mxreg" argument:
    mc.and.ar.length <- length(object.new$call$mc)+length(object.new$call$ar)
    if(NCOL(object$aux$mX) > mc.and.ar.length){
      object.new$call$mxreg <- "mxreg"
    }else{
      object.new$call$mxreg <- NULL
    }

  } else {

    object.new$call$mc <- NULL
    object.new$call$ar <- NULL
    ##to do: log.ewma
    object.new$call$mxreg <- NULL

  }


  ##-----------------------------------
  ## pass on arguments to predict.arx:
  ##-----------------------------------

  out <- predict.arx(object.new, spec="mean", n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=NULL, newindex=newindex,
    return=return, plot=plot, plot.options=plot.options)

  ##-------------------
  ## return forecasts:
  ##-------------------

  if(return){ return(out) }

}
