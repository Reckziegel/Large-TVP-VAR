carter_kohn <-  function(y, Z, R, Qt, m, p, t, B0, V0) {
    # Carter and Kohn (1994), "On Gibbs sampling for state space models".
    
    library(mvtnorm) # install.packages("mvtnorm")
    
    # Kalman Filter
    bp <- B0
    Vp <- V0
    bt <- matrix(data = 0, nrow = t, ncol = m)
    Vt <- matrix(data = 0, nrow = m ^ 2, ncol = t)
    log_lik <- 0
    
    for (i in 1:t) {
        H <- Z[((i - 1) * p + 1):(i * p), ]
        # F = eye(m);
        cfe <-  y[ , i] - H %*% bp   # conditional forecast error
        f <- H %*% Vp %*% t(H) + R   # variance of the conditional forecast error
        K <- (Vp %*% t(H)) / f
        log_lik <- log_lik + log(det(f)) + t(cfe) %*% solve(f) %*% cfe
        btt <- bp + K %*% cfe
        Vtt <- Vp - K %*% (H %*% Vp)
        if (i < t) {
            bp <- btt
            Vp <- Vtt + Qt
        }
        bt[i, ] <- t(btt)
        Vt[ , i] <- matrix(data = Vtt, nrow = m ^ 2, ncol = 1)
    }

    # draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
    bdraw <- matrix(data = 0, nrow = t, ncol = m)
    bdraw[t, ] <- mvnrnd(btt, Vtt, 1) # This is one of Koops additional functions

    # Backward recurssions
    for (i in 1:(t - 1)) {
        bf <- t(bdraw[t - i + 1, ])
        btt <- t(bt[t - i, ])
        Vtt <- matrix(data = Vt[ , t - i], nrow = m, ncol = m)
        f <- Vtt + Qt
        K <- Vtt / f
        cfe <- bf - btt
        bmean <- btt + K %*% cfe
        bvar <- Vtt - K %*% Vtt
        bdraw[t - i, ] <- mvnrnd(bmean, bvar, 1)  #mvnrnd(bmean, bvar, 1)
    }
    bdraw <- t(bdraw)
    
    out <- list(bdraw = bdraw, log_lik = log_lik)
    return(out)
    
}
