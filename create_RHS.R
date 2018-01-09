# function [x_t,K] = create_RHS(YY,M,p,t)
# 
# K = M + p*(M^2); % K is the number of elements in the state vector
# 
# % Create x_t matrix.
# % first find the zeros in matrix x_t
# x_t = zeros((t-p)*M,K,'single');
# for i = 1:t-p
# ztemp = eye(M);
# for j = 1:p
# xtemp = YY(i,(j-1)*M+1:j*M);
# xtemp = kron(eye(M),xtemp);
# ztemp = [ztemp xtemp];
# end
# x_t((i-1)*M+1:i*M,:) = ztemp;
# end

create_RHS <- function(YY, M, p, t) {
    
    K <- M + p * (M ^ 2) # K is the number of elements in the state vector
    x_t <- matrix(data = 0, nrow = (t - p) * M, ncol = K) # repeat (t-p)*M 'K' times
    
    for (i in seq(from = 1, to = t - p, by = 1)) {
        ztemp <- diag(M)
        for (j in seq_along(p)) {
            xtemp <- YY[i, ((j - 1) * M + 1):(j * M)]
            xtemp <- t(kronecker(diag(M), xtemp))
            ztemp <- cbind(ztemp, xtemp)
        }
        x_t[ ((i - 1) * M + 1):(i * M), ] <- ztemp
    }
    return(x_t) 
}
