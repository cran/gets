regs.var <-
function(e, vc=TRUE, arch=NULL, asym=NULL,
  log.ewma=NULL, vxreg=NULL, p=2, zero.adj=0.1)
{
#logep:
logep <- glog.ep(e, zero.adj=zero.adj, p=p)
mX <- cbind(logep)
colnames(mX) <- if(p==2){"loge2"}else{"logep"}

#variance constant:
e.n <- length(e)
if(!is.null(vc)){vconst <- rep(1,e.n)}else{vconst <- NULL}
mX <- cbind(mX,vconst)

#create exponential arch terms:
archlags <- NULL
if(is.null(arch)){ archadj<-NULL ; s2<-1 }else{
  archadj <- setdiff(round(abs(arch)),0)
  for(i in 1:length(archadj))
    {
    archlags <- cbind(archlags,glag(logep,k=archadj[i]))
    }
  colnames(archlags) <- paste(c("arch"),archadj,sep="")
  s2 <- I(max(archadj)+1)
}
mX <- cbind(mX, archlags)

#create asymmetry terms:
asymlags <- NULL
if(is.null(asym)){asymadj <- NULL}else
  {
  asymadj <- setdiff(round(abs(asym)),0)
  for(i in 1:length(asymadj))
    {
    asymlags <- cbind(asymlags,glag(logep*(e < 0),k=asymadj[i]))
    }
  colnames(asymlags) <- paste(c("asym"),asymadj,sep="")
  s2 <- max(s2,I(max(asym)+1))
  }
mX <- cbind(mX, asymlags)

#create log-EWMA term:
if(is.null(log.ewma)){logEWMA <- NULL}else{
  logEWMA <- do.call(leqwma, c(list(e),log.ewma) )
}
mX <- cbind(mX, logEWMA)

#create matrix of variance regressors vxreg:
if(is.null(vxreg)){vX <- NULL}else{
 vxregnames <- colnames(vxreg)
 if(is.null(vxregnames))
   {
   vxregnames <- paste("vxreg",1:ncol(vxreg),sep="")
   }else
   {
   for(i in 1:length(vxregnames))
     {
     if(vxregnames[i]==""){vxregnames[i] <- paste("vxreg",i,sep="")}
     }
   }
 vX <- cbind(vxreg)
 colnames(vX) <- vxregnames
 }
mX <- cbind(mX, vX)

#out-matrix:
out <- mX
return(out)

}
