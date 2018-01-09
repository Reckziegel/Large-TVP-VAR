function y = standardize1(x,T_thres)

% standardize using mean and variance from training sample
% Data in x are Txp, i.e. T time series observations times p variables
y = (x - repmat(mean(x(1:T_thres,:),1),size(x,1),1))./repmat(std(x(1:T_thres,:),1),size(x,1),1);
