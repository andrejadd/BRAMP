##########################################
#mvrnorm, function de la librairie MASS
##########################################

mvrnorm<-function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) 
{
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE, EISPACK = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1]))) {
      cat("FIXME: in mvrnorm(), 'Sigma' is not positive definite:")
      print.table(Sigma)
      stop("stopping")
    }
    
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
        t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
        nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
        drop(X)
    else t(X)
}
