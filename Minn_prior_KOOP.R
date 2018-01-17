Minn_prior_KOOP <- function(alpha_bar, gamma, M, p, K) {
    
    # This is the version of the Minnesota prior with no dependence on the
    # standard deviations of the univariate regressions. This prior allows
    # online estimation and forecasting of the large TVP-VAR.
    
    # 1. Minnesota Mean on VAR regression coefficients
    A_prior <- t(matrix(data = c(matrix(data = 0 , nrow = 1, ncol = M),  
                                 0 * diag(M), 
                                 matrix(data = 0, nrow = (p - 1) * M, ncol = M)), 
                        ncol = M, 
                        byrow = TRUE))
    a_prior <- as.vector(A_prior) 
    
    # 2. Minnesota Variance on VAR regression coefficients
    # Create an array of dimensions K x M, which will contain the K diagonal   
    # elements of the covariance matrix, in each of the M equations.
    V_i <- matrix(data = 0, nrow = K / M, ncol = M)
    for (i in 1:M) {                       # for each i-th equation
        for (j in 1:(K / M)) {             # for each j-th RHS variable
            if (j == 1) {                  # if there is constant, use this code
                    V_i[j, i] <- alpha_bar # variance on intercept
            } else {
                V_i[j, i] <- gamma / ((ceiling((j - 1) / M)) ^ 2) # variance on own lags           
                # Note: 
                # the "ceil((j-1)/M^2)" command finds the associated lag number for each parameter         
            }
        }
    }
    # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
    V_i_T <-  t(V_i)
    V_prior <- diag(as.vector(V_i_T))  # this is the prior variance of the vector alpha
    # This last command is wrong for sure. I will verify how to fix it latter.
    
    out <- list(a_prior = a_prior, V_prior = V_prior)
    return(out)
    
}

