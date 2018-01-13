# function y = standardize1(x,T_thres)
# 
# % standardize using mean and variance from training sample
# % Data in x are Txp, i.e. T time series observations times p variables
# y = (x - repmat(mean(x(1:T_thres,:),1),size(x,1),1))./repmat(std(x(1:T_thres,:),1),size(x,1),1);

standardize1 <- function(x, T_thres) {
    # standardize using mean and variance from training sample
    # Data in x are Txp, i.e. T time series observations times p variables
    
    y <- (x - repmat(t(apply(x[1:T_thres, ], 2, mean)), nrow(x), 1)) / repmat(t(apply(x[1:T_thres, ], 2, sd)), nrow(x), 1)
    y
    
    #scale(x[1:T_thres, ], center = TRUE)
    
    # Still have problems with this one!
    # Answers are different then those provided by matlab.

}

