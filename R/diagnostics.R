diagnostics <-
function(x, s2=1, y=NULL, xreg=NULL,
  ar.LjungB=c(1,0.025), arch.LjungB=c(1,0.025),
  normality.JarqueB=NULL, verbose=TRUE, user.fun=NULL)
{
  diagnosticsGood <- TRUE
  if(s2 == 1){ zhat <- x }else{ zhat <- x/sqrt(s2) }

  ##serial correlation:
  if(!is.null(ar.LjungB)){
    ar.LjungBox <- Box.test(zhat, lag=ar.LjungB[1], type="L")
    if(ar.LjungBox$p.value <= ar.LjungB[2]){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end if(!is.null(..))

  ##arch:
  if(diagnosticsGood && !is.null(arch.LjungB)){
    zhat2 <- zhat^2
    arch.LjungBox <- Box.test(zhat2, lag=arch.LjungB[1], type="L")
    if(arch.LjungBox$p.value <= arch.LjungB[2]){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end arch

  ##normality:
  if(diagnosticsGood && !is.null(normality.JarqueB)){
    n <- length(zhat)
    avgzhat <- mean(zhat)
    zhat.avgzhat <- zhat-avgzhat
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
  if(diagnosticsGood && !is.null(user.fun)){
    userVals <- do.call(user.fun$name, list(y, xreg, x, zhat),
      envir=.GlobalEnv)
    if( length(user.fun)>=2 ){
      userFunPval <- as.numeric(rbind(userVals)[,3])
      if( any(userFunPval <= user.fun[[2]]) ){
        diagnosticsGood <- FALSE
      }
    }
  }

  ##if(verbose): diagnostics table
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
      result <- rbind(result, as.numeric(userVals))
      resultRowNames <- c(resultRowNames, user.fun$name)
    }
    if(!is.null(result)){
      rownames(result) <- resultRowNames
      colnames(result) <- c("Chi-sq", "df", "p-value")
      if(!is.null(user.fun)){ colnames(result)[1] <- "Statistic" }
    }
  }

  ##if(!verbose): logical
  if(!verbose){ result <- diagnosticsGood }

  ##return result:
  return(result)

}
