predict.gets <-
function(object, spec=NULL, n.ahead=12,
  newmxreg=NULL, newvxreg=NULL, newindex=NULL,
  n.sim=1000, innov=NULL, return=TRUE, plot=TRUE,
  plot.options=list(), ...)
{
  ##create new object to add stuff to in order to use predict.arx()
  object.new <- object

  ##spec:
  if(is.null(spec)){
    if(object$gets.type=="getsm"){ spec <- "mean" }
    if(object$gets.type=="getsv"){ spec <- "variance" }
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  } #end if(..)else(..)
  if( spec %in% c("variance","both") ){ stop("variance-forecasts with predict.gets not available yet") }


  ##-----------------------------------
  ## arguments mean-equation:
  ##-----------------------------------

  if("mX" %in% names(object$aux)) {
    ##what dynamics specified in gum?
    gum.ar <- eval(object$aux$call.gum$ar)
    ##what dynamics remain in spec?
    spec.ar <- as.numeric(gsub("ar(\\d+)","\\1",object$aux$mXnames[grep("^ar\\d+$",object$aux$mXnames)]))

    if(NROW(spec.ar)==0) {
      object.new$call$ar <- NULL
    } else { ##check that dynamics in specific are subset of those in gum
      object.new$call$ar <- spec.ar[spec.ar %in% gum.ar]
    }

    ##need to check for constant
    if(!is.null(object$aux$call.gum$mc)) {
      ##but the constant needs to have been retained...
      if ( "mconst" %in% object$aux$mXnames ) { object.new$call$mc <- TRUE }
    }

    if(is.null(object$aux$call.gum$mc) && "mconst" %in% object$aux$mXnames) {
      ##ensure that when we call predict.arx we set mc=TRUE
      object.new$call$mc <- !is.null(object$aux$call.gum$mc)
      ##need to also get rid of mconst in mX
      #object.new$mX <- object$aux$mX[,-grep("mconst",object$aux$mXnames)]
    }

    ##mxreg needs to be null if only thing in it is dynamics (i.e. ar variables)
    if(NROW(grep("^ar\\d+$",object$aux$mXnames))==object$aux$mXncol) {
      object.new$call$mxreg <- NULL
    } else {
      object.new$call$mxreg <- object$aux$call.gum$mxreg
    } ##if there is non-arch variables in mXreg, user needs to specify newmxreg

  } else {

    object.new$call$mc <- NULL
    object.new$call$ar <- NULL
    ##to do: log.ewma
    object.new$call$mxreg <- NULL

  }

  ##-----------------------------------
  ## arguments variance-equation:
  ##-----------------------------------

#  if("vX" %in% names(object$aux)) {
#
#    ##what dynamics specified in gum?
#    gum.arch <- eval(object$aux$call.gum$arch)
#    ##what dynamics remain in spec?
#    spec.arch <- as.numeric(gsub("arch(\\d+)","\\1",object$aux$vXnames[grep("^arch\\d+$",object$aux$vXnames)]))
#    if(NROW(spec.arch)==0) {
#      object.new$call$arch <- NULL
#    } else { ##check that dynamics in specific are subset of those in gum
#      object.new$call$arch <- spec.arch[spec.arch %in% gum.arch]
#    }
#
#    ##which asyms specified in gum?
#    gum.asym <- eval(object$aux$call.gum$asym)
#    ##which asyms remain in spec?
#    spec.asym <- as.numeric(gsub("asym(\\d+)","\\1",object$aux$vXnames[grep("^asym\\d+$",object$aux$vXnames)]))
#    if(NROW(spec.asym)==0) {
#      object.new$call$asym <- NULL
#    } else { ##check that dynamics in specific are subset of those in gum
#      object.new$call$asym <- spec.asym[spec.asym %in% gum.asym]
#    }
#
#    ##which asyms specified in gum?
#    gum.log.ewma <- eval(object$aux$call.gum$log.ewma)
#    if(!is.null(gum.log.ewma)){ stop("'log.ewma' not available yet for 'predict.gets' ") }
#
#    ##vxreg needs to be null if only thing in it is dynamics (i.e. arch variables)
#    if(NROW(grep("^arch\\d+$",object$aux$vXnames))==object$aux$vXncol) {
#      object.new$call$vxreg <- NULL
#    } else {
#      object.new$call$vxreg <- object$aux$call.gum$vxreg
#    } ##if there is non-arch variables in vXreg, user needs to specify newvxreg
#
#  } else {

    object.new$call$arch <- NULL
    object.new$call$asym <- NULL
    object.new$call$log.ewma <- NULL
    object.new$call$vxreg <- NULL

#  }

  ##-----------------------------------
  ## pass on arguments to predict.arx:
  ##-----------------------------------

  out <- predict.arx(object.new, spec=spec, n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=newvxreg, newindex=newindex,
    n.sim=n.sim, innov=innov, return=return, plot=plot,
    plot.options=plot.options)

  ##-------------------
  ## return forecasts:
  ##-------------------

  if(return){ return(out) }

}
