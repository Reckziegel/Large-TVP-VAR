Minn_prior_LITT <- function(Y, Ylag, alpha_bar, gamma, M, p, K, t) {
    
    # 1. Minnesota Mean on VAR regression coefficients
    A_prior <- t(matrix(data = c(matrix(data = 0 , nrow = 1, ncol = M),  
                                 0 * diag(M), 
                                 matrix(data = 0, nrow = (p - 1) * M, ncol = M)), 
                        ncol = M, 
                        byrow = TRUE))
    # A_prior = [mean(Y(1:40,:)); 0.9*eye(M); zeros((p-1)*M,M)]';
    # A_prior = [zeros(1,M) ; 0.9*repmat(eye(M),p,1)]';
    a_prior <- as.vector(A_prior)
    
    # 2. Minnesota Variance on VAR regression coefficients
    # Now get residual variances of univariate p-lag autoregressions. Here we just run the AR(p) model on each equation,
    # ignoring the constant and exogenous variables (if they have been specified for the original VAR model)
    p_MIN <- p
    sigma_sq <- matrix(data = 0, nrow = M, ncol = 1) # vector to store residual variances

    for (i in 1:M) {
        Ylag_i <- Ylag
        X_i <- cbind(matrix(data = 1, nrow = t - p_MIN + p, ncol = 1),  Ylag_i[ , seq(i, M * p_MIN, M)])
        Y_i <- Y[ , i]
        # OLS estimates of i-th equation
        alpha_i <- solve(t(X_i) %*% X_i) %*% (t(X_i) %*% Y_i)
        sigma_sq[i, 1] <-  (1 / (t - p_MIN + 1)) %*% t((Y_i - X_i %*% alpha_i)) %*% (Y_i - X_i %*% alpha_i)
    }
    
    # Now define prior hyperparameters. 
    # Create an array of dimensions K x M, which will contain the K diagonal   
    # elements of the covariance matrix, in each of the M equations.
    V_i <- matrix(data = 0, nrow = K / M, ncol = M)
    
    # index in each equation which are the own lags
    ind <- matrix(data = 0, nrow = M, ncol = p)
    for (i in 1:M) {
        ind[i, ] <- seq(i + 1, K / M, M)
    }
    for (i in 1:M) {            # for each i-th equation
        for (j in 1:(K / M)) {  # for each j-th RHS variable   
            if (j == 1) {       # if there is constant, use this code
                V_i[j, i] <- alpha_bar %*% sigma_sq[i, 1] # variance on intercept               
            } else if (which(j == ind[i, ]) > 0) {
                V_i[j, i] <- gamma / ((ceiling(j - 1) / M) ^ 2) # variance on own lags
                # Note: the "ceiling((j-1)/M)" command finds the associated lag number for each parameter
            } else {
                for (kj in 1:M) {
                    if (which(j == ind[kj, ]) > 0) {
                        ll <- kj
                    }
                }
                # variance on other lags
                V_i[j, i] <- (gamma %*% sigma_sq[i, 1]) / (((ceiling(j - 1) / M) ^ 2) %*% sigma_sq(ll, 1))
            }
        }
    }
    
    # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
    V_i_T <- t(V_i)
    # this is the prior variance of the vector alpha
    V_prior <-  diag(as.vector(V_i_T))  
    Sigma_0 <- diag(sigma_sq)
    
    out <- list(a_prior = a_prior, V_prior = V_prior, Sigma_0 = Sigma_0)
    return(out)
    
}
