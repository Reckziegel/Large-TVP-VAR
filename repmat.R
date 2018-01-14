repmat <- function(x, m, n) {
    
    if (class(data) != 'matrix') {
        x <- as.matrix(x)
    }
    
    mx <- dim(x)[1]
    nx <- dim(x)[2]
    matrix(data  = t(matrix(data = x, nrow = mx, ncol = nx * n)), 
           nrow  = mx * m, 
           ncol  = nx * n, 
           byrow = TRUE)
    
}
