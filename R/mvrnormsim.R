mvrnormsim <-
function(n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE){
  
  p <- length(mu)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) {
    nm <- dn[[1L]]  
  }
  dimnames(X) <- list(nm, NULL)
  if (n == 1) {
    drop(X) } else {
      t(X)   
    } 
  #return(X)
}
