mvnfactor <- function(sigma) {     
    # originally [T,p] = 
    # MVNFACTOR  Do Cholesky-like decomposition, allowing zero eigenvalues
    # SIGMA must be symmetric.  In general T is not square or triangular.
    
    U <- eigen(sigma + t(sigma))$vectors / 2
    D <- eigen(sigma + t(sigma))$values

    tol <- max(D) * length(D) * eps # check where "eps" comes form. I do not remember.
    t <- D > tol
    D <- D(t)
    T <- diag(sqrt(D)) %*% t(U[ , t])
    return(T)
}
