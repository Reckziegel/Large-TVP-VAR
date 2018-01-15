mvnrnd <- function(mu, sigma, cases) {
    # modified matlab code, lets almost singular matricies, that may be generated at some draws, go through

    # MVNRND Random matrices from the multivariate normal distribution.
    #   R = MVNRND(MU,SIGMA,CASES) returns a matrix of random numbers chosen from the multivariate 
    # normal distribution with mean vector, MU, and covariance matrix, SIGMA. CASES is the number 
    # of rows in R.
    #
    #  SIGMA is a symmetric positive semi-definite matrix with size equal to the length of MU.
    
    m1 <- nrow(mu)
    n1 <- ncol(mu)
    c  <- max(m1, n1)
    
    m <- nrow(sigma)
    n <- ncol(sigma)
    
    #safety check
    if (m1 * n1 != c) {
        stop('Mu must be a vector.')
    }
    
    if (m != n) {
        stop('Sigma must be square')
    }
    
    if (m != c) {
        stop('The length of mu must equal the number of rows in sigma.')
    }
    
    if (any(eigen(sigma, only.values = TRUE)$values %in% TRUE)) { # if it is positive definite
        TT <- chol(sigma)
    } else {
    # The input covariance has some type of perfect correlation.
    # If it is positive semi-definite, we can still proceed.
    # Find a factor T that has the property: sigma == T' * T.
    # Unlike a Cholesky factor, this T is not necessarily triangular or even square.
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
            TT <- -sqrt(D) * t(U[ , t])
            #return(TT)
        }
        TT <- mvnfactor(sigma)
    # if p > 0, error('Sigma must be a non-negative definite matrix.'); end
    }
    if (m1 == c) {
        mu <- t(mu)
    }
    
    mu <- mu[matrix(data = 1, nrow = cases, ncol = 1), ] 
    r <- matrix(data = rnorm(cases * nrow(TT)), nrow = cases, byrow = TRUE) %*% TT + mu 
    return(r)
    
}

