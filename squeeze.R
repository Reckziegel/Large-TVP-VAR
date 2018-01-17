squeeze <- function (x)  {
    # The same as the squeeze function in MATLAB 
    d <- dim(x)
    dim(x) <- d[d > 1]
    x
}
