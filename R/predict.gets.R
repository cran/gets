predict.gets <-
function(object, spec=NULL, n.ahead=12,
  newmxreg=NULL, newvxreg=NULL, newindex=NULL,
  n.sim=5000, innov=NULL, probs=NULL, ci.levels=NULL, 
  quantile.type=7, return=TRUE, verbose=FALSE, plot=NULL,
  plot.options=list(), ...)  
{

  ##create new object to add stuff to in order to use predict.arx()
  objectNew <- object


  ##-----------------------------------
  ## arguments mean-equation:
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
    gumTerms <- eval(object$aux$call.gum$ar)
    gumNamesAr <- paste0("ar", gumTerms)
    whichRetained <- which( gumNamesAr %in% names(coefsMean) )
    if( length(whichRetained)==0 ){
      objectNew$call$ar <- NULL
    }else{
      objectNew$call$ar <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
        
    ##ewma argument:
    gumTerms <- eval(object$aux$call.gum$ewma)
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
#      objectNew$call$mxreg <- xreg
    }

  } #end if( length(coefsMean)>0 )
  

  ##-----------------------------------
  ## arguments variance-equation:
  ##-----------------------------------

  ##coefficients of variance spec in final model:
  coefsVar <- coef.arx(objectNew, spec="variance")
  if( length(coefsVar)>0 ){ #remove Elnz2 estimate:
    coefsVar <- coefsVar[ -length(coefsVar) ]  
  }

  ##there is no variance equation:
  if( length(coefsVar)==0 ){

    objectNew$call$vc <- NULL
    objectNew$call$arch <- NULL
    objectNew$call$asym <- NULL
    objectNew$call$log.ewma <- NULL
    objectNew$call$vxreg <- NULL

  }

  ##there is a variance equation:
  if( length(coefsVar)>0 ){

    ##vc argument (always present in variance equations):
    objectNew$call$vc <- TRUE
    indxCounter <- 1 #used for vxreg
    
    ##arch argument:
    gumTerms <- eval(object$aux$call.gum$arch)
    gumNamesArch <- paste0("arch", gumTerms)
    whichRetained <- which( gumNamesArch %in% names(coefsVar) )
    if( length(whichRetained)==0 ){
      objectNew$call$arch <- NULL
    }else{
      objectNew$call$arch <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
    
    ##asym argument:
    gumTerms <- eval(object$aux$call.gum$asym)
    gumNamesAsym <- paste0("asym", gumTerms)
    whichRetained <- which( gumNamesAsym %in% names(coefsVar) )
    if( length(whichRetained)==0 ){
      objectNew$call$asym <- NULL
    }else{
      objectNew$call$asym <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
    
    ##log.ewma argument:
    gumTerms <- eval(object$aux$call.gum$log.ewma)
    gumNamesLogEwma <- paste0("logEqWMA(", gumTerms$length, ")")
    whichRetained <- which( gumNamesLogEwma %in% names(coefsVar) )
    if( length(whichRetained)==0 ){
      objectNew$call$log.ewma <- NULL
    }else{
      objectNew$call$log.ewma <-
        list( length=gumTerms$length[ whichRetained ] )
      indxCounter <- indxCounter + length(whichRetained)
    }

    ##vxreg argument:
    whichRetainedCoefs <- coefsVar[ -c(1:indxCounter) ]
    if( length(whichRetainedCoefs)==0 ){
      objectNew$call$vxreg <- NULL
    }else{
      whichRetainedNames <- names(whichRetainedCoefs)
      objectNew$call$vxreg <- whichRetainedNames
#more correct (but not needed, since vxreg only needs to be non-NULL)?:
#      whichRetained <- which( object$aux$vXnames %in% whichRetainedNames )
#      vxreg <- cbind(object$aux$vX[, whichRetained ])
#      colnames(vxreg) <- whichRetainedNames
#      objectNew$call$vxreg <- vxreg
    }

  } #end if( length(coefsVar)>0 )


  ##----------------------------------
  ## pass arguments on to predict.arx:
  ##----------------------------------

  result <- predict.arx(objectNew, spec=spec, n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=newvxreg, newindex=newindex,
    n.sim=n.sim, innov=innov, probs=probs, ci.levels=ci.levels,
    quantile.type=quantile.type, return=return, verbose=verbose,
    plot=plot, plot.options=plot.options)

  ##-------------------
  ## return forecasts:
  ##-------------------

  if(return){ return(result) }

}
