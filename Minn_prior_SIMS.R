Minn_prior_SIMS <- function(Y, Ylag, M, p, K, t) {
    # 1. Minnesota Mean on VAR regression coefficients
    a <- 2 - sqrt(3)
    s <- matrix(data = 0, nrow = p + 1, ncol = 1)
    for (i in 0:p) {
        s[i + 1, ] <- (1 + a) * ((-a) ^ i)
    }
    A_prior <- t(matrix(data = c(s[1] * matrix(data = 1, nrow = 1, ncol = M), 
                                 kronecker(s[2:(p + 1)], diag(M))), 
                        ncol = M, 
                        byrow = TRUE))
    a_prior <- as.vector(A_prior)
    
    # 2. Minnesota Variance on VAR regression coefficients. First define the hyperparameters 'pi_i'
    p_1 <- 0.001
    p_2 <- 6
    p_3 <- 10
    # Now get residual variances of univariate p-lag autoregressions. Here we just run the AR(p) model on each equation, 
    # ignoring the constant and exogenous variables (if they have been specified for the original VAR model)
    p_MIN <- p
    sigma_sq <- matrix(data = 0, nrow = M, ncol = 1) # vector to store residual variances
    for (i in 1:M) {
        # Create lags of dependent variables
        Ylag_i <- Ylag[ , seq(from = i, to = M * p_MIN, by = M)]
        X_i <- cbind(matrix(data = 1, nrow = t - p_MIN + p, ncol = 1),  Ylag_i)
        Y_i = Y[ , i]
        # OLS estimates of i-th equation
        alpha_i <- solve(t(X_i) %*% X_i) %*% (t(X_i) %*% Y_i)
        sigma_sq[i, 1] <- (1 / (t - p_MIN + 1)) %*% t(Y_i - X_i %*% alpha_i) %*% (Y_i - X_i %*% alpha_i)
    }
    
    # Now define prior hyperparameters.
    # Create an array of dimensions K x M, which will contain the K diagonal   
    # elements of the covariance matrix, in each of the M equations.
    V_i <- matrix(data = 0, nrow = K / M, ncol = M)
    
    # index in each equation which are the own lags
    ind <- matrix(data = 0, nrow = M, ncol = p)
    for (i in 1:M) {
        ind[i, ] <- seq(from = i + 1, to = K / M, by = M)
    }
    
    for (i in 1:M) {                                           # for each i-th equation
        for (j in 1:(K / M)) {                                 # for each j-th RHS variable   
            if (j == 1) {                                      # if there is constant, use this code
                V_i[j, i] <- exp(-p_2 * log(1 / p_3))          # variance on intercept         
                #V_i(j,i) = sigma_sq(i,1)*exp(-p_2*log(1/p_3)) # variance on intercept
            } else if (which(j == ind[i, ]) > 0) {
                V_i[j, i] <- p_1 * exp(-p_2 * log(ceiling((j - 1) / M))) # variance on own lags           
                # Note: the "ceil((j-1)/M)" command finds the associated lag number for each parameter
            } else {
                for (kj in 1:M) {
                    if (which(j == ind[kj, ]) > 0) {
                        ll <- kj
                    }
                } # variance on other lags
                #V_i(j,i) = 0.1*exp(-p_2*log(ceil((j-1)/M)))
                V_i[j, i] <- (sigma_sq[i, 1] / sigma_sq[ll, 1]) * p_1 * exp(-p_2 * log(ceiling((j - 1) / M)))
            }
        }
    }
    
    # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
    V_i_T <- t(V_i)
    V_prior <- diag(as.vector(V_i_T))  # this is the prior variance of the vector alpha
    Sigma_0 <- diag(sigma_sq)
    
    out <- list(a_prior = a_prior, V_prior = V_prior, Sigma_0 = Sigma_0)
    return(out)
    
}
    
    