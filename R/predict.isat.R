predict.isat <-
function(object, n.ahead=12,
  newmxreg=NULL, newindex=NULL, return=TRUE, plot=NULL,
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
    ##what dynamics remain in specific?
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

    ##if sis and tis terms retained need to adapt mxreg call ...
    if(!is.null(object$ISnames)) {

      ##J-dog, add your code here??
      if(is.null(object$call$mxreg)) { #need to ensure predict.arx knows there are mx variables
        object.new$call$mxreg <- "mXis"
      }

      ##... and need to specify newmxreg of right dimension:
      ##if we're here it means the isat call did not specify any mxregs
      ##hence what is in mX in the object are terms retained by isat
      ##we can automatically create isat terms into sample period (exception uis)
      if(is.null(newmxreg)) {
        ##if no newmxreg specified we add iis/sis/tis from scratch

        ##first check that there shouldn't be something in newmxreg...
        if(!is.null(object$call$mxreg)){ stop("'newmxreg' is NULL") }

        ##assuming not, then we start from scratch adding the indicators...
        newmxreg <- c()
      }

      if(any(regexpr("^iis",object$ISnames)>-1)){##isat retained some iis terms
        for(i in object$ISnames[grep("^iis",object$ISnames)]) {
          newmxreg <- cbind(newmxreg,rep(0,n.ahead))
        }
      }
      if(any(regexpr("^sis",object$ISnames)>-1)){##isat retained some sis terms
        for(i in object$ISnames[grep("^sis",object$ISnames)]) {
          newmxreg <- cbind(newmxreg,rep(1,n.ahead))
        }
      }
      if(any(regexpr("^tis",object$ISnames)>-1)){##isat retained some tis terms
        for(i in object$ISnames[grep("^tis",object$ISnames)]) {
          newmxreg <- cbind(newmxreg,
                            seq(1,n.ahead)+object$aux$mX[NROW(object$aux$mX),i])
        }
      }

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
