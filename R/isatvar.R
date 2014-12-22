isatvar <-
function(x, ...)
{

  if (!is.null(x$mean.fit)){
    if (!is.null(x$ISnames)){
      var.rel <- c( which(substr(x$aux$mXnames,1,6) %in% "mconst"), which(x$aux$mXnames %in% x$ISnames))  #vector of where the constant and is terms are
    } else {
      var.rel <- which(substr(x$aux$mXnames,1,6) %in% "mconst")
    }

    vcov.rel <- x$vcov.mean[var.rel,var.rel]
    dim.var <- NCOL(vcov.rel)
    const.var <- matrix(NA, dim.var, 1 )

    #construct a matrix to multiply by the variances
    dim.in <- NCOL(x$aux$mX[,var.rel])
    indic.mat <- x$aux$mX[,var.rel]
    const.mat <- matrix(NA, NROW(x$aux$mX[,var.rel]), dim.in)


    if (dim.var > 1) #if there are indicators retained
    {

      #order the breaks in correct order so the covariance matrix and relevance order is correct
      order.mat <- apply(indic.mat[,], 2, function(x) min(which(x==1)))  #find where each indicator is first one, for sorting

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
    const.varse <- cbind(ind.var.mat, ind.se.mat)
    colnames(const.varse) <- c("const.var", "const.se")
    const.varse <- zoo(const.varse , order.by=x$aux$y.index)
    return(const.varse)

  } ##if (is null) closed

}
