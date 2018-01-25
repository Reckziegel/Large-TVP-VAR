# DMA Results
msfe_focus_DMA <- apply(MSFE_DMA, 2, mean, na.rm = TRUE); msfe_focus_DMA
mafe_focus_DMA <- apply(MAFE_DMA, 2, mean, na.rm = TRUE); mafe_focus_DMA
lpls_focus_DMA <- apply(logpl_DMA, 2, sum, na.rm = TRUE); lpls_focus_DMA
lpl_focus_DMA  <- sum(LOG_PL_DMA, na.rm = TRUE); lpl_focus_DMA

# # DMS Results
msfe_focus_DMS <- apply(MSFE_DMS, 2, mean, na.rm = TRUE); msfe_focus_DMS
mafe_focus_DMS <- apply(MAFE_DMS, 2, mean, na.rm = TRUE); mafe_focus_DMS
lpls_focus_DMS <- apply(logpl_DMS, 2, sum, na.rm = TRUE); lpls_focus_DMS
lpl_focus_DMS  <- sum(LOG_PL_DMS, na.rm = TRUE); lpl_focus_DMS

# plot volatilities and lambda_t (but not theta_t which is huge)
library(zoo)

prob_pl <- matrix(data = 0, nrow = nrow(index_DMA), ncol = nos)
for (i in 1:nos) {
    prob_pl[ , i] <- omega_update[[i]][index_DMA[ , i]]
}

final_prob <- matrix(data = 0, nrow = t, ncol = nos)
for (sss in 1:nos) {
    final_prob[ , sss] <- prob_pl[ ,sss] / apply(prob_pl ,1 , sum)
}

lam_pl <- prob_pl <- matrix(data = 0, nrow = nrow(index_DMA), ncol = nos)
for (i in 1:nos) {
    lam_pl[ , i] <- lambda_t[[i]][index_DMA[ , i]]
}

autoplot.zoo(final_prob, facets = NULL)
autoplot.zoo(lam_pl, facets = NULL) # forgetting must be 2 to this graph vary.

# plot predictive likelihoods over time
autoplot.zoo(cumsum(LOG_PL_DMS)) # Cumulative sum of total predictive likelihood h=1

# print the Minnesota prior over time
vars_Min <- T_thres:t
for (ss in 1:nos) {
    vars_Min <- cbind(vars_Min, Minn_gamms[[ss]])
}
colnames(vars_Min) <- c('Period', 'Small VAR', 'Medium VAR', 'largeVAR')
autoplot.zoo(vars_Min[, -1])
