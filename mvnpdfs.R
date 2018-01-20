mvnpdfs <- function(X, Mu = NULL, Sigma = NULL) {
    # MVNPDF Multivariate normal probability density function (pdf).
    # Y = MVNPDF(X) returns the probability density of the multivariate normal distribution with zero mean and identity covariance matrix, evaluated at
    # each row of X. Rows of the N-by-D matrix X correspond to observations or points, and columns correspond to variables or coordinates. Y is an
    # N-by-1 vector.
    #
    # Y = MVNPDF(X,MU) returns the density of the multivariate normal distribution with mean MU and identity covariance matrix, evaluated
    # at each row of X.  MU is a 1-by-D vector, or an N-by-D matrix, in which case the density is evaluated for each row of X with the corresponding
    # row of MU.  MU can also be a scalar value, which MVNPDF replicates to match the size of X.
    #
    # Y = MVNPDF(X,MU,SIGMA) returns the density of the multivariate normal distribution with mean MU and covariance SIGMA, evaluated at each row
    # of X.  SIGMA is a D-by-D matrix, or an D-by-D-by-N array, in which case the density is evaluated for each row of X with the corresponding page
    # of SIGMA, i.e., MVNPDF computes Y(I) using X(I,:) and SIGMA(:,:,I). If the covariance matrix is diagonal, containing variances along the 
    # diagonal and zero covariances off the diagonal, SIGMA may also be specified as a 1-by-D matrix or a 1-by-D-by-N array, containing 
    # just the diagonal. Pass in the empty matrix for MU to use its default value when you want to only specify SIGMA.
    #
    # If X is a 1-by-D vector, MVNPDF replicates it to match the leading dimension of MU or the trailing dimension of SIGMA.
    #
    # Example:
    #    mu = [1 -1]; Sigma = [.9 .4; .4 .3];
    #    [X1,X2] = meshgrid(linspace(-1,3,25)', linspace(-3,1,25)');
    #    X = [X1(:) X2(:)];
    #    p = mvnpdf(X, mu, Sigma);
    #    surf(X1,X2,reshape(p,25,25));
    
    if (nargs() < 1) {
        stop('stats:mvnpdf:TooFewInputs','Requires at least one input.')
    } else if (length(dim(X)) != 2) {
        stop('stats:mvnpdf:InvalidData','X must be a matrix.')
    }
    
    # Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
    n <- nrow(X)
    d <- ncol(X)
    if (d < 1) {
        stop('stats:mvnpdf:TooFewDimensions. X must have at least one column.')
    }
    
    # Assume zero mean, data are already centered
    if (nargs() < 2 || is.null(Mu) == TRUE) {
        X0 <- X
    # Get scalar mean, and use it to center data
    } else if (length(Mu) == 1) {
        X0 <- X - Mu
    # Get vector mean, and use it to center data
    } else if (length(dim(Mu)) == 2) {
        n2 <- nrow(Mu)
        d2 <- ncol(Mu)
        if (d2 != d) { # has to have same number of coords as X
            stop('stats:mvnpdf:InputSizeMismatch. X and MU must have the same number of columns.')
        } else if (n2 == n) {# lengths match
            X0 <- X - Mu
            X0 <- matrix(rep(X0, d2), nrow = n2, ncol = d2, byrow = TRUE)
        } else if (n2 == 1) {# mean is a single row, rep it out to match data
            X0 <- sweep(X, 2, Mu, '-') # bsxfun(@minus,X,Mu);
            X0 <- matrix(rep(X0, d2), nrow = n2, ncol = d2, byrow = TRUE)
        } else if (n == 1) { # data is a single row, rep it out to match mean
            n <- n2
            X0 <- sweep(X, 2, Mu, '-') # bsxfun(@minus,X,Mu);  
            X0 <- matrix(rep(X0, d2), nrow = n, ncol = d2, byrow = TRUE)
        } else {# sizes don't match
            stop('stats:mvnpdf:InputSizeMismatch. X or MU must be a row vector, or X and MU must have the same number of rows.')
        }
    } else {
        stop('stats:mvnpdf:BadMu','MU must be a matrix.')
    }
    
    # Assume identity covariance, data are already standardized
    if (nargs() < 3 || is.null(Sigma) == TRUE) {
        # Special case: if Sigma isn't supplied, then interpret X
        # and Mu as row vectors if they were both column vectors
        if (d == 1 & length(X) > 1) {
            X0 <- t(X0)
            d <- ncol(X0)
        }
        xRinv <- X0
        logSqrtDetSigma <- 0

    # Single covariance matrix
    } else if (length(dim(Sigma)) == 2) {
        sz <- dim(Sigma)
        if (sz[1] == 1 & sz[2] > 1) {
            # Just the diagonal of Sigma has been passed in.
            sz[1] <- sz[2]
            sigmaIsDiag <- TRUE
        } else {
            sigmaIsDiag <- FALSE
        }

        # Special case: if Sigma is supplied, then use it to try to interpret X and Mu as row vectors if they were both column vectors.
        if ((d == 1) && (length(X) > 1) && (sz[1] == n)) {
            X0 <- t(X0)
            d <- ncol(X0)
        }
    
        #Check that sigma is the right size
        if (sz[1] != sz[2]) {
            stop('stats:mvnpdf:BadCovariance. SIGMA must be a square matrix or a vector.')
        } else if (all(sz != d)) {
            stop('stats:mvnpdf:InputSizeMismatch. SIGMA must be a square matrix with size equal to the', '/n',
                  'number of columns in X, or a row vector with length equal to the number of columns in X.')
        } else {
            if (sigmaIsDiag == TRUE) {
                if (apply(Sigma <= 0, 2, any)) { # any(Sigma<=0)
                    stop('stats:mvnpdf:BadDiagSigma. All elements of diagonal SIGMA must be positive.')
                }
                R <- sqrt(Sigma)
                xRinv <- sweep(X0, MARGIN = 2, R, FUN = '/') # bsxfun(@rdivide, X0, R)
                logSqrtDetSigma <- sum(log(R))
            } else {
                # Make sure Sigma is a valid covariance matrix
                #             [R,err] = cholcov(Sigma,0);
                #             if err ~= 0
                #                 error('stats:mvnpdf:BadCovariance', ...
                #                     'SIGMA must be symmetric and positive definite.');
                #             end
                # Create array of standardized data, and compute log(sqrt(det(Sigma)))
                R <- chol(Sigma)
                xRinv <- X0 %*% solve(R)
                logSqrtDetSigma <- sum(log(diag(R)))
            }
        }
    } else if (length(dim(Sigma)) == 3) { # Multiple covariance matrices
        sz <- dim(Sigma)
        if (sz[1] == 1 & sz[2] > 1) {
            # Just the diagonal of Sigma has been passed in.
            sz[1] <- sz[2]
            Sigma <- t(matrix(data = Sigma, nrow = sz[2], ncol = sz[3]))
            sigmaIsDiag <- TRUE
        } else {
            sigmaIsDiag <- FALSE
        }
        
        # Special case: if Sigma is supplied, then use it to try to interpret
        # X and Mu as row vectors if they were both column vectors.
        if ((d == 1 & length(X) > 1) & (sz[1] == n)) {
            X0 <- t(X0)
            n <- nrow(X0)
            d <- ncol(X0)
        }
        
        # Data and mean are a single row, rep them out to match covariance
        if (n == 1) { # already know size(Sigma,3) > 1
            n <- sz[3]
            X0 <- repmat(X0, n, 1) # rep centered data out to match cov
        }
    
        # Make sure Sigma is the right size
        if (sz[1] != sz[2]) {
            stop('stats:mvnpdf:BadCovariance. Each page of SIGMA must be a square matrix or a row vector.')
        } else if ((sz(1) != d) | (sz(2) != d)) { # Sigma is a stack of dxd matrices
            stop('stats:mvnpdf:InputSizeMismatch. Each page of SIGMA must be a square matrix', '/n', 
                  'with size equal to the number of columns in X, or a row vector with length', '/n', 
                  'the same as the number of columns in X.')
        } else if (sz[3] != n) {
            stop('stats:mvnpdf:InputSizeMismatch. SIGMA must have one page for each row of X.')
        } else {
            if (sigmaIsDiag == TRUE) {
                if (any(apply(Sigma <= 0, 2, any))) {
                    stop('stats:mvnpdf:BadDiagSigma. All elements of diagonal SIGMA must be positive.')
                }
                R <- sqrt(Sigma)
                xRinv <- X0  / R
                logSqrtDetSigma <- apply(log(R), 1, sum)
            } else {
                # Create array of standardized data, and vector of log(sqrt(det(Sigma)))
                xRinv <- matrix(data = 0, nrow = n, ncol = d) #superiorfloat(X0,Sigma))
                logSqrtDetSigma <- matrix(data = 0, nrow = n, ncol = 1)  # class(Sigma))
                for (i in 1:n) {
                    #                 # Make sure Sigma is a valid covariance matrix
                    #                 [R,err] = cholcov(Sigma(:,:,i),0);
                    #                 if err ~= 0
                    #                     error('stats:mvnpdf:BadCovariance',...
                    #                         'Each page of SIGMA must be symmetric and positive definite.');
                    #                 end
                    xRinv[i, ] <- X0[i, ] / chol(Sigma)
                    logSqrtDetSigma[i, ] = sum(log(diag(R)))
                }
            }
        }
    } else if (length(dim(Sigma)) > 3) {
        stop('stats:mvnpdf:BadCovariance. SIGMA must be a vector, a matrix, or a 3 dimensional array.')
    }
    
    # The quadratic form is the inner products of the standardized data
    quadform <- apply(xRinv ^ 2, 1, sum) 
    
    y <- exp(-0.5 * quadform - logSqrtDetSigma - d * log(2 * pi) / 2)
    print(y)
    
}
    