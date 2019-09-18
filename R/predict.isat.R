predict.isat <-
function(object, n.ahead=12, newmxreg=NULL,
  newindex=NULL, n.sim=2000, probs=NULL, ci.levels=NULL,
  quantile.type=7, return=TRUE, verbose=FALSE, plot=NULL,
  plot.options=list(), ...)

{
  ## 1 arguments of mean-equation
  ## 2 make plot.options argument
  ## 3 pass arguments on to predict.arx
  ## 4 return forecasts

  ##create new object to add stuff to in order to use predict.arx()
  objectNew <- object

  ##-----------------------------------
  ## 1 arguments of mean-equation:
  ##-----------------------------------

  ##coefficients of mean spec in final model:
  coefsMean <- coef.arx(objectNew, spec="mean")

  ##there is no mean equation:
  if( length(coefsMean)==0 ){

    objectNew$call$mc <- NULL
    objectNew$call$ar <- NULL
    objectNew$call$ewma <- NULL
    objectNew$call$mxreg <- NULL

  }

  ##there is a mean equation:
  if( length(coefsMean)>0 ){

    ##initiate index counter (used for mxreg):
    indxCounter <- 0

    ##mc argument:
    mconstRetained <- "mconst" %in% names(coefsMean)
    if( mconstRetained ){
      objectNew$call$mc <- TRUE
      indxCounter <- indxCounter + 1
    }else{
      objectNew$call$mc <- NULL
    }

    ##ar argument:
    gumTerms <- eval(object$call$ar)
#OLD:
#    gumTerms <- eval(object$aux$call.gum$ar)
    gumNamesAr <- paste0("ar", gumTerms)
    whichRetained <- which( gumNamesAr %in% names(coefsMean) )
    if( length(whichRetained)==0 ){
      objectNew$call$ar <- NULL
    }else{
      objectNew$call$ar <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }

    ##ewma argument:
    gumTerms <- eval(object$call$ewma)
#OLD:
#    gumTerms <- eval(object$aux$call.gum$ewma)
    gumNamesEwma <- paste0("EqWMA(", gumTerms$length, ")")
    whichRetained <- which( gumNamesEwma %in% names(coefsMean) )
    if( length(whichRetained)==0 ){
      objectNew$call$ewma <- NULL
    }else{
      objectNew$call$ewma <-
        list( length=gumTerms$length[ whichRetained ] )
      indxCounter <- indxCounter + length(whichRetained)
    }

    ##mxreg argument:
    if(indxCounter==0){ whichRetainedCoefs <- coefsMean }
    if(indxCounter>0){ whichRetainedCoefs <- coefsMean[ -c(1:indxCounter) ] }
    if( length(whichRetainedCoefs)==0 ){
      objectNew$call$mxreg <- NULL
    }else{
      whichRetainedNames <- names(whichRetainedCoefs)
      objectNew$call$mxreg <- whichRetainedNames
#more correct (but not needed, since mxreg only needs to be non-NULL)?:
#      whichRetained <- which( object$aux$mXnames %in% whichRetainedNames )
#      mxreg <- cbind(object$aux$mX[, whichRetained ])
#      colnames(mxreg) <- whichRetainedNames
#      objectNew$call$mxreg <- mxreg
    }

  } #end if( length(coefsMean)>0 )

  ##here: introduce the automated detection of indicators and how they
  ##should modify newmxreg?

  ##----------------------------------
  ## 2 make plot.options argument:
  ##----------------------------------
  
  if(is.null(plot.options$start.at.origin)){
    plot.options$start.at.origin <- FALSE
  }
  if(is.null(plot.options$line.at.origin)){
    plot.options$line.at.origin <- TRUE
  }
  if(is.null(plot.options$fitted)){
    plot.options$fitted <- TRUE
  }
  
  
  ##-------------------------------------
  ## 3 pass arguments on to predict.arx:
  ##-------------------------------------

  innov <- rnorm(n.ahead*n.sim) #force normal errors
  result <- predict.arx(objectNew, spec="mean", n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=NULL, newindex=newindex,
    n.sim=n.sim, innov=innov, probs=probs, ci.levels=ci.levels,
    quantile.type=quantile.type, return=return, verbose=verbose,
    plot=plot, plot.options=plot.options)

  ##---------------------
  ## 4 return forecasts:
  ##---------------------

  if(return){ return(result) }

}
