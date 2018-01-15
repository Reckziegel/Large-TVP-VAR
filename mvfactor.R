mvnfactor <- function(sigma) {     
    # originally [T,p] = 
    # MVNFACTOR  Do Cholesky-like decomposition, allowing sigmaero eigenvalues
    # SIGMA must be symmetric.  In general T is not square or triangular.
    
    eps = 0.01
    
    U <- eigen(sigma + t(sigma))$vectors
    U <- U[ , order(ncol(U):1)] 
    D <- eigen(sigma + t(sigma))$values / 2
    D <- rev(D)

    tol <- max(D) * length(D) * eps # check where "eps" comes form. I do not remember.
    t <- D > tol
    D <- D[t]
    T <- -sqrt(D) * t(U[ , t])
    return(T)
}
