mlag2 <- function(X, p) { 
    
    Traw <- nrow(X)
    N <- ncol(X)
    
    Xlag <- matrix(data = 0, nrow = Traw, ncol = N * p)
    #colnames(Xlag) <- names(X)
    
    for (ii in 1:p) {
        Xlag[(p + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <- X[(p + 1 - ii):(Traw - ii), 1:N]
    }
    
    return(Xlag)
    
}
