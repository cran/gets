isatdates <-
function(x){
  
  mxbreak_iis <- "iis"
  mxbreak_sis <- "sis"
  mxbreak_tis <- "tis"
  
  iis.names <- c(x$ISnames[grep(mxbreak_iis, x$ISnames)])
  sis.names <- c(x$ISnames[grep(mxbreak_sis, x$ISnames)])
  tis.names <- c(x$ISnames[grep(mxbreak_tis, x$ISnames)])
  
  ##### iis
  if(length(iis.names) != 0){
    iis.breaks <- data.frame(matrix(NA, nrow=NROW(iis.names), ncol=1))
    names(iis.breaks) <- c("breaks")  
    is.m <- as.matrix(x$aux$mX[,iis.names])
    is.index <- which(x$aux$mXnames %in% iis.names)
    colnames(is.m) <- iis.names
    iis.date <- apply(is.m,2, function(x) (which(x>0))[1])  
    iis.date.index <- iis.date
    iis.date <- x$aux$y.index[iis.date]
    
    iis.breaks$breaks <- iis.names
    iis.breaks$date <- iis.date
    iis.breaks$index <- iis.date.index
    iis.breaks$coef <- x$mean.results$coef[is.index]
    iis.breaks$coef.se <- x$mean.results$std.error[is.index]
    iis.breaks$coef.t <- x$mean.results$`t-stat`[is.index]
    iis.breaks$coef.p <- x$mean.results$`p-value`[is.index]
    
  } else {
    iis.breaks <- NULL 
  }
  
  ##### sis
  if(length(sis.names) != 0){
    sis.breaks <- data.frame(matrix(NA, nrow=NROW(sis.names), ncol=1))
    names(sis.breaks) <- c("breaks")  
    sis.m <- as.matrix(x$aux$mX[,sis.names])
    sis.index <- which(x$aux$mXnames %in% sis.names)
    colnames(sis.m) <- sis.names
    sis.date <- apply(sis.m,2, function(x) (which(x>0))[1])  
    sis.date.index <- sis.date
    sis.date <- x$aux$y.index[sis.date]
    
    sis.breaks$breaks <- sis.names
    sis.breaks$date <- sis.date
    sis.breaks$index <- sis.date.index
    sis.breaks$coef <- x$mean.results$coef[sis.index]
    sis.breaks$coef.se <- x$mean.results$std.error[sis.index]
    sis.breaks$coef.t <- x$mean.results$`t-stat`[sis.index]
    sis.breaks$coef.p <- x$mean.results$`p-value`[sis.index]
    
  } else {
    sis.breaks <- NULL 
  }
  
  ##### tis
  if(length(tis.names) != 0){
    tis.breaks <- data.frame(matrix(NA, nrow=NROW(tis.names), ncol=1))
    names(tis.breaks) <- c("breaks")  
    tis.m <- as.matrix(x$aux$mX[,tis.names])
    tis.index <- which(x$aux$mXnames %in% tis.names)
    colnames(tis.m) <- tis.names
    tis.date <- apply(tis.m,2, function(x) (which(x>0))[1])  
    tis.date.index <- tis.date
    tis.date <- x$aux$y.index[tis.date]
    
    tis.breaks$breaks <- tis.names
    tis.breaks$date <- tis.date
    tis.breaks$index <- tis.date.index
    tis.breaks$coef <- x$mean.results$coef[tis.index]
    tis.breaks$coef.se <- x$mean.results$std.error[tis.index]
    tis.breaks$coef.t <- x$mean.results$`t-stat`[tis.index]
    tis.breaks$coef.p <- x$mean.results$`p-value`[tis.index]
    
  } else {
    tis.breaks <- NULL 
  }
  
  out <- list(iis.breaks, sis.breaks, tis.breaks)
  names(out) <- c("iis", "sis", "tis")
  return(out)
  
}
