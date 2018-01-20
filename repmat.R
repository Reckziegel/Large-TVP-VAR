repmat <- function(x, m, n) {
    
    mx <- dim(x)[1]
    nx <- dim(x)[2]
    matrix(data  = t(matrix(data = x, nrow = mx, ncol = nx * n)), 
           nrow  = mx * m, 
           ncol  = nx * n, 
           byrow = TRUE)
    
}
