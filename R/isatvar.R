isatvar <-
function(x, lr=FALSE, conscorr=FALSE, effcorr=FALSE, mcor = 1,  mxfull = NULL, mxbreak=NULL)
{
  
  
  
  if (lr == TRUE && !is.null(mxfull)){
    print("Warning: LR currently not defined with user-specified break variables. LR set to FALSE")
    lr <- FALSE
  }
  
  
  
  if(!is.null(as.list(x$call)$uis) && is.null(mxfull)){
    
    if (!is.null(as.list(x$call)$iis) || !is.null(as.list(x$call)$sis) ){
      print("Warning: uis specified but no mxfull variable given. Using mconst instead.")
      mxfull <- "mconst"
    } else {
      stop("uis specified but no mxfull variable given")
    }
    
    
  }
  
  
  if (is.null(mxfull)){   ##if no full variable is specified, then the constant is used
    mxfull <- "mconst"
  }
  
  
  if (!is.null(x$mean.fit)){
    
    if (!is.null(mxbreak))
    {
      
      ISnames <- x$ISnames[grep(mxbreak, x$ISnames)]
      ####################################
    } else {
      
      mxbreak_iis <- "iis"
      mxbreak_sis <- "sis"

      ISnames <- c(x$ISnames[grep(mxbreak_iis, x$ISnames)], x$ISnames[grep(mxbreak_sis, x$ISnames)])
    }
    
    
    if (!is.null(ISnames)){
      
      
      var.rel <- c( which(substr(x$aux$mXnames,1,6) %in% mxfull), which(x$aux$mXnames %in% ISnames))
      
      
      if (conscorr==TRUE){
        x$vcov.mean <- x$vcov.mean*as.numeric(isvarcor(x$aux$t.pval,1)[2]^2)
      }
      
      if (effcorr==TRUE){
        
        if (!is.null(x$keep)) {
          x$vcov.mean[x$keep, x$keep] <- x$vcov.mean[x$keep, x$keep] * as.numeric(isvareffcor(x$aux$t.pval, 1, mcor)[2]^2)
        }
      }
      
      #coefficient path function
      is.x <- cbind(x$aux$mX[, x$aux$mXnames %in% ISnames])
      is.coef.ests <- coef.isat(x)[ISnames]
      
      if(is.null(mxfull)){
        coef.path <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index) 
      } else {
        is.x[is.x!=0] <- 1 #the break variable is not necessarily an indicator variable
        coef.path <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
      }
      
      
      
      
      colnames(coef.path) <- "coef.path"
      
      const <- coef.isat(x)[mxfull]
      const.path <- coef.path + const
      colnames(const.path) <- "const.path"
      
    } else {
      var.rel <- which(substr(x$aux$mXnames,1,6) %in% mxfull)
      
      #coefficient path function
      const <- coef.isat(x)[mxfull]
      
      coef.path <- zoo(0, order.by = x$aux$y.index)
      
      const.path <- coef.path + const
      
    }
    
    vcov.rel <- x$vcov.mean[var.rel,var.rel]
    
    
    dim.var <- NCOL(vcov.rel)
    const.var <- matrix(NA, dim.var, 1 )
    
    #construct a matrix to multiply by the variances
    dim.in <- NCOL(x$aux$mX[,var.rel])
    
    const.mat <- matrix(NA, NROW(x$aux$mX[,var.rel]), dim.in)
    
    if(is.null(mxfull)){
      indic.mat <- x$aux$mX[,var.rel]
      
    } else {
      indic.mat <- x$aux$mX[,var.rel]
      indic.mat[indic.mat!=0] <- 1
      
    }
    
    
    
    if (dim.var > 1) #if there are indicators retained
    {
      
      #order the breaks in correct order so the covariance matrix and relevance order is correct
      
      if(is.null(mxfull)){
        order.mat <- apply(indic.mat[,], 2, function(x) min(which(x==1)))  #find where each indicator is first one, for sorting
        
      } else { #if using a custom break variable
        order.mat <- apply(indic.mat[,], 2, function(x) min(which(x!=0)))
      }
      
      
      
      
      indic.mat <- indic.mat[,order(order.mat)]  #order the indicators
      vcov.rel <- vcov.rel[order(order.mat),order(order.mat)] #order the covariance matrix
      
      for (i in 1:dim.var)
      {
        const.var[i] <- sum(vcov.rel[1:i,1:i])   #sum over the expanding variance covariance matrix to get the variance of the sums of coefficients
      }
      
      const.mat <- indic.mat
      
      for (j in 2:(dim.in))
      {
        
        const.mat[which(rowSums(as.matrix(const.mat[,(j):(dim.in)]))>0),j-1] <- 0     #puts zeros in the appropriate places to make sure the correct s.e. is applied for each point in time
        
      }
      
      ind.var.mat <- const.mat %*% const.var   #the variance as it applies to each subsection
      
    } else { #just the constant remains
      
      const.var <- vcov.rel
      const.mat <- indic.mat
      ind.var.mat <- const.mat * const.var
      
    }
    
    ind.se.mat <- sqrt(ind.var.mat)   #the standard errors of the coefficient as it changes with the SIS breaks
    
    ####the coefficient path of the LR mean
    
    if(lr){
      if (!is.null(as.list(x$call)$ar)){
        
        if (!is.null(x$mean.fit)){
          
          vcov.rel.tot <- x$vcov.mean
          
          
          
          coef.rel <- coef.isat(x)
          
          ###coefficient path of LR mean
          
          arcall <- as.list(x$call)$ar
          arnames <- paste("ar",eval(as.expression(arcall)), sep="")
          
          ar.coefs <- coef.isat(x)[arnames]
          ar.sum <- sum(ar.coefs)
          lr.path <- const.path/(1-ar.sum)
          
          ###variance of the coefficient path
          
          ar.var <- vcov.rel.tot[arnames,arnames] #variance of ar() terms
          armu.cov <- vcov.rel.tot[c(mxfull, ISnames), arnames] #covariance of ar(1) and const + sis
          
          
          if (!is.null(ISnames)) {
            var.rel <- c(which(substr(x$aux$mXnames, 1, 6) %in%
                                 mxfull), which(x$aux$mXnames %in% ISnames))
          } else {
            var.rel <- which(substr(x$aux$mXnames, 1, 6) %in%
                               mxfull)
          }
          
          vcov.rel <- x$vcov.mean[var.rel, var.rel]
          dim.var <- NCOL(vcov.rel)
          const.var <- matrix(NA, dim.var, 1)
          dim.in <- NCOL(x$aux$mX[, var.rel])
          indic.mat <- x$aux$mX[, var.rel]
          const.mat <- matrix(NA, NROW(x$aux$mX[, var.rel]), dim.in)
          armu.cov.sum <- matrix(NA, dim.var, 1)
          
          dim.ar <- NCOL(ar.var)
          
          if (dim.var > 1) {
            order.mat <- apply(indic.mat[, ], 2, function(x) min(which(x ==
                                                                         1)))
            indic.mat <- indic.mat[, order(order.mat)]
            vcov.rel <- vcov.rel[order(order.mat), order(order.mat)]   #ordering the variables for cumulative summing for covariances of sis
            
            #do the same for the autoregressive terms
            
            if (dim.ar > 1) #if more than one ar term
            {
              armu.cov <- armu.cov[order(order.mat), ]
              
            } else {
              armu.cov <- armu.cov[order(order.mat)]
            }
            
            for (i in 1:dim.var) {
              const.var[i] <- sum(vcov.rel[1:i, 1:i])
              if (dim.ar > 1) #if more than one ar term
              {
                armu.cov.sum[i] <- sum(armu.cov[1:i,])
              } else {
                armu.cov.sum[i] <- sum(armu.cov[1:i])
              }
            }
            
            const.mat <- indic.mat
            for (j in 2:(dim.in)) {
              const.mat[which(rowSums(as.matrix(const.mat[,
                                                          (j):(dim.in)])) > 0), j - 1] <- 0
            }
            
            #variance of ar
            ar.var.sum <- sum(ar.var)
            ar.var.part <- ((const.path^2)/((1-ar.sum)^4))*ar.var.sum
            
            #covariance and variance part of sis
            ind.var.mat <- const.mat %*% const.var # this one should be ok
            covar.is.part <- (1/(1-ar.sum)^2)*ind.var.mat
            
            #covariance part of ar and sis
            ind.ar.covar.mat <- const.mat %*% armu.cov.sum
            covar.aris.part <- 2*((const.path)/(1-ar.sum)^2)*(1/(1-ar.sum))*ind.ar.covar.mat
            
            lr.mean.var <- ar.var.part + covar.is.part + covar.aris.part
            lr.mean.se <- sqrt(lr.mean.var)
            
          } else { #if no steps retained
            
            ar.var.sum <- sum(ar.var)
            ar.var.part <- ((const.path^2)/((1-ar.sum)^4))*ar.var.sum
            
            ind.var.mat <-  vcov.rel
            covar.is.part <- (1/(1-ar.sum)^2)*ind.var.mat
            
            #covariance part of ar and sis
            
            armu.cov.sum <- sum(armu.cov)
            ind.ar.covar.mat <- armu.cov.sum
            covar.aris.part <- 2*((const.path)/(1-ar.sum)^2)*(1/(1-ar.sum))*ind.ar.covar.mat
            
            lr.mean.var <- ar.var.part + covar.is.part + covar.aris.part
            lr.mean.se <- sqrt(lr.mean.var)
            
          }
          
        } #if isnullxmean closed
        
        if(!is.null(as.list(x$call)$tis)){
          if (as.list(x$call)$tis) { #if TIS
            
            ind.var.mat <- NA
            ind.se.mat <- NA
            lr.path <- NA
            lr.mean.var <- NA
            lr.mean.se <- NA
          }
        }
        const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat, lr.path, lr.mean.var, lr.mean.se)
        colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se", "lr.path", "lr.var", "lr.se")
        
      }  else { #if there are no ar
        
        if(!is.null(as.list(x$call)$tis)){
          if (as.list(x$call)$tis) { #if TIS
            ind.var.mat <- NA
            ind.se.mat <- NA
          }
        }
        const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat)
        colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se")
        
      }  #if no ar closed
    } else {   #if no lr
      
      if(!is.null(as.list(x$call)$tis)){
        if (as.list(x$call)$tis) {
          ind.var.mat <- NA
          ind.se.mat <- NA
        }
      }
      const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat)
      colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se")
      
      
      
    }#if lr closed
    
    const.varse <- zoo(const.varse , order.by=x$aux$y.index)
    return(const.varse)
    
  } ##if (is null) closed
  
}
