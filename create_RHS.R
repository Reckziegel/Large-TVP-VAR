create_RHS <- function(YY, M, p, t) {
    
    K <- M + p * (M ^ 2) # K is the number of elements in the state vector
    x_t <- matrix(data = 0, nrow = (t - p) * M, ncol = K) # repeat (t-p)*M 'K' times
    
    for (i in 1:(t - p)) {
        ztemp <- diag(M)
        for (j in 1:p) {
            xtemp <- t(YY[i, ((j - 1) * M + 1):(j * M)])
            xtemp <- kronecker(diag(M), xtemp)
            ztemp <- c(ztemp, xtemp)
        }
        x_t[((i - 1) * M + 1):(i * M), ] <- ztemp
    }
    
    out <- list(x_t = x_t, K = K)
    return(out)
    
}
