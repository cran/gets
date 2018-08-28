diagnostics <-
function(x, ar.LjungB=c(1,0.025), arch.LjungB=c(1,0.025),
  normality.JarqueB=NULL, verbose=TRUE, user.fun=NULL, ...)
{
  ##initiate:
  if( is.null(x$std.residuals) ){
    zhat <- x$residuals
  }else{
    zhat <- x$std.residuals
  }
  diagnosticsGood <- TRUE

  ##serial correlation:
  if( !is.null(ar.LjungB) ){
    ar.LjungBox <- Box.test(zhat, lag=ar.LjungB[1], type="L")
    if( ar.LjungBox$p.value <= ar.LjungB[2] ){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end serial correlation

  ##arch:
  if( diagnosticsGood && !is.null(arch.LjungB) ){
    zhat2 <- zhat^2
    arch.LjungBox <- Box.test(zhat2, lag=arch.LjungB[1], type="L")
    if(arch.LjungBox$p.value <= arch.LjungB[2]){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end arch

  ##normality:
  if( diagnosticsGood && !is.null(normality.JarqueB)){
    n <- length(zhat)
    avgzhat <- mean(zhat) #do I really need this?
    zhat.avgzhat <- zhat-avgzhat #do I really need this?
    zhat.avgzhat2 <- zhat.avgzhat^2
    K <- n*sum(zhat.avgzhat^4)/(sum(zhat.avgzhat2)^2)
    S <- (sum(zhat.avgzhat^3)/n)/(sum(zhat.avgzhat2)/n)^(3/2)
    JB <- (n/6)*(S^2 + 0.25*((K-3)^2))
    JBpval <- pchisq(JB, df = 2, lower.tail=FALSE)
    if(JBpval <= normality.JarqueB){
        diagnosticsGood <- FALSE
        diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end normality

  ##user.fun:
  if( diagnosticsGood && !is.null(user.fun) ){
    userVals <- do.call(user.fun$name, x, envir=.GlobalEnv)
    userVals <- rbind(userVals)
    if( length(user.fun)>=2 ){
      userFunPval <- as.numeric(userVals[,3])
      if( any(userFunPval <= user.fun[[2]]) ){
        diagnosticsGood <- FALSE
      }
    }
  }

  ##if(verbose): return diagnostics table
  if(verbose){
    result <- NULL
    resultRowNames <- NULL
    if(exists("ar.LjungBox")){
      tmp <- as.numeric(ar.LjungBox[1:3])
      resultRowNames <- c(resultRowNames,
        paste("Ljung-Box AR(", tmp[2], ")", sep=""))
      result <- rbind(result, tmp)
    }
    if(exists("arch.LjungBox")){
      tmp <- as.numeric(arch.LjungBox[1:3])
      resultRowNames <- c(resultRowNames,
        paste("Ljung-Box ARCH(", tmp[2], ")", sep=""))
      result <- rbind(result, tmp)
    }
    if(exists("JBpval")){
      tmp <- c(JB, 2, JBpval)
      resultRowNames <- c(resultRowNames,
        paste("Jarque-Bera", sep=""))
      result <- rbind(result, tmp)
    }
    if(exists("userVals")){
      result <- rbind(result, userVals)
      userValsNames <- rownames(userVals)
      if( identical(userValsNames, "userVals") ){
        userValsNames <- user.fun$name
      }
      if( is.null(userValsNames) ){
        userValsNames <- rep(user.fun$name, NROW(userVals))
      }
      resultRowNames <- c(resultRowNames, userValsNames)
    }
    if(!is.null(result)){
      rownames(result) <- resultRowNames
      colnames(result) <- c("Chi-sq", "df", "p-value")
      if(!is.null(user.fun)){ colnames(result)[1] <- "Statistic" }
    }
  } #end verbose

  ##if(!verbose): return logical
  if(!verbose){ result <- diagnosticsGood }

  ##return result:
  return(result)

}
