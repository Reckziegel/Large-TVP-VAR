# function [f] = kernel_PL(Z,X,bandw)
# %Kernel density estimator
# nobs = size(X,1);
# n = size(Z,1); %#ok<NASGU>
# Q1 = quantile(X,0.25);
# Q3 = quantile(X,0.75);
# 
# if nargin==2
# h=0.79*(Q3-Q1)*nobs^(-0.2); 
# else
#     h = bandw;
# end
# 
# v1 = normpdf((Z - X)/h);
# f = mean(v1)/h;

kernel_PL <- function(Z, X, bandw = 0) {
    
    #Kernel density estimator
    if (class(X) == 'numeric') {
        nobs <- length(X)
    } else if (class(X) == 'matrix') {
        nobs <- nrow(X)
    } else {
            stop("you should continue working with this function!")
    }
    
    #n <- nrow(Z)
    Q1 <- round(quantile(X, 0.25), 0)
    Q3 <- round(quantile(X, 0.75), 0)

    if (bandw == 0) {
        h <- 0.79 * (Q3 - Q1) * round(nobs ^ (-0.2), 0)
    } else {
        h <- bandw
    }
    
    v1 <- dnorm((Z - X) / h)
    f <- mean(v1) / h
    #return(f)
}


