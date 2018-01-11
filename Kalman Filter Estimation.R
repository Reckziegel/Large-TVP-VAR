# ======================= BEGIN KALMAN FILTER ESTIMATION =======================

for (irep in 1:t) {
    for (ss in 1:nos) {  # LOOP FOR 1 TO NOS VAR MODELS OF DIFFERENT DIMENSIONS
        # Find sum of probabilities for DPS
        if (irep > 1) {
            sum_prob_omega[[ss]] <- sum((omega_update[[ss]][irep - 1, ]) ^ eta)  
            # this is the sum of the nom model probabilities 
            # (all in the power of the forgetting factor 'eta')
        }
        
        for (k in 1:nom) { # LOOP FOR 1 TO NOM VAR MODELS WITH DIFFERENT DEGREE OF SHRINKAGE
            
            # Predict
            if (irep == 1) {
                theta_pred[[ss]][ , irep, k] <- theta_0_prmean[[ss]][ , k]         
                R_t[[ss]][ , , k] <- theta_0_prvar[[ss]][ , , k]           
                omega_predict[[ss]][irep, k] <- 1 / nom
            } else {
                theta_pred[[ss]][ ,irep, k] <- theta_update[[ss]][ , irep - 1, k]
                R_t[[ss]][ , , k] <- (1 / lambda_t[[ss]][irep - 1, k]) %*% S_t[[ss]][ , , k]
                omega_predict[[ss]][irep, k] <- ((omega_update[[ss]][irep - 1, k]) ^ eta + offset) /
                    (sum_prob_omega[[ss]] + offset)
            }
            xx = x_t[[ss]][((irep - 1) %*% M(ss) + 1):(irep %*% M(ss)), ]
            y_t_pred[[ss]][ , irep, k] <- xx %*% theta_pred[[ss]][ , irep, k] 
            # this is one step ahead prediction

            # Prediction error
            e_t[[ss]][ , irep, k] <- t(y_t[[ss]][irep, ]) - y_t_pred[[ss]][ , irep, k]
            # this is one step ahead prediction error

            # Update forgetting factor
            if (forgetting == 2) {
                lambda_t[[ss]][irep, k] <- lambda_min + 
                    (1 - lambda_min) %*% (LL ^ (-round(alpha %*% t(e_t[[ss]][1:nfocus, irep, k]) %*% e_t[[ss]][1:nfocus, irep, k])))
            }
            
            # first update V[t], see the part below equation (10)
            A_t <- e_t[[ss]][ , irep, k] %*% t(e_t[[ss]][ , irep, k])
            if (irep == 1) {
                V_t[[ss]][ , , irep, k] <- kappa * Sigma_0[[ss]]
            } else {
                V_t[[ss]][ , , irep, k] <- kappa * squeeze(V_t[[ss]][ , , irep - 1, k]) + (1 - kappa) * A_t
            }
            
            #%         if all(eig(squeeze(V_t(:,:,irep,k))) < 0)
            #%             V_t(:,:,irep,k) = V_t(:,:,irep-1,k);       
            #%         end

            # update theta[t] and S[t]
            Rx <- R_t[[ss]][ , , k] %*% t(xx)
            KV <- squeeze(V_t[[ss]][ , , irep, k]) + xx %*% Rx 
            kG = Rx / KV
            theta_update[[ss]][ , irep, k] <- theta_pred[[ss]][ , irep, k] + (KG %*% e_t[[ss, 1]][ , irep, k])
            S_t[[ss]][ , , k] <- R_t[[ss]][ , ,k] - KG %*% (xx %*% R_t[[ss]][ , , k])

            # Find predictive likelihood based on Kalman filter and update DPS weights     
            if (irep == 1) {
                variance <- Sigma_0[[ss]][1:nfocus, 1:nfocus] + xx[1:nfocus, ] %*% Rx[ , 1:nfocus]
            } else {
                variance <- V_t[[ss]][1:nfocus, 1:nfocus, irep, k] + xx[1:nfocus, ] %*% Rx[ , 1:nfocus]
            }
            if (which((eig(variance) > 0) == 0) > 0) {
                variance <- abs(diag(diag(variance)))
            }
            ymean <- y_t_pred[[ss]][1:nfocus, irep, k]
            ytemp <- t(y_t[[ss]][irep, 1:nfocus])
            f_l[k, 1] <- mvnpdfs(ytemp, ymean, variance)
            w_t[[ss]][ , irep, k] <- omega_predict[[ss]][irep, k] %*% f_l[k, 1]

            # forecast simulation for h=1 only (instead of saving only the
            # mean and variance, I just generate 2000 samples from the
            # predictive density
            variance2 <- V_t[[ss]][ , , irep, k] + xx %*% Rx
            if (which((eig(variance2) > 0) == 0) > 0) { 
                variance2 <- abs(diag(diag(variance2)))
            }
            bbtemp <- theta_update[[ss]][(M[ss] + 1):K[ss], irep, k]  # get the draw of B(t) at time i=1,...,T  (exclude intercept)                   
            splace <- 0
            biga <- 0
            for (ii in 1:p) {
                for (iii in 1:M[ss]) {            
                    biga[iii, ((ii - 1) * M[ss] + 1):(ii * M[ss])] <- t(bbtemp(splace + 1:splace + M[ss], 1))
                    splace <- splace + M[ss]
                }
            }
            beta_fore <- rbind(t(theta_update[[ss]][1:M[ss], irep, k]), t(biga))
            ymean2 <- rbind(1 , y_t[[ss]][irep, ], x_f[[ss]][irep, 1:M[ss] * (p - 1)] %*% t(beta_fore))
            y_fore[[ss]][1, , k, ] <- repmat(ymean2, 1, 2000) + t(chol(variance2)) * randn(M[ss], 2000) # This line is wrong. Check latter.
        } # End cycling through nom models with different shrinkage factors

        # First calculate the denominator of Equation (19) (the sum of the w's)
        sum_w_t <- 0   
        for (k_2 in 1:nom) {       
            sum_w_t <- sum_w_t + w_t[[ss]][ , irep, k_2]
        }
              
        # Then calculate the DPS probabilities  
        for (k_3 in 1:nom) {
            omega_update[[ss]][irep, k_3] <- (w_t[[ss]][ , irep, k_3] + offset) / (sum_w_t + offset)  # this is Equation (19)
        }
        [max_prob[[ss]], k_max[[ss]]] <- max(omega_update[[ss]][irep, ])
        index_DMA[irep, ss] <- k_max[[ss]]
                                    
        # Use predictive likelihood of best (DPS) model, and fight the weight for DMA     
        w2_t[irep, ss] <- omega_predict[[ss]][irep, k_max[ss]] * f_l[k_max[ss], 1]
                                                
        # Find "observed" out-of-sample data for MSFE and MAFE calculations
        if (irep <= (t - nfore)) {
            Yraw_f[[ss]] <- y_t[[ss]][(irep + 1):(irep + nfore) , , ] #Pseudo out-of-sample observations                   
        } else {
            Yraw_f[[ss]] <- rbind(y_t[[ss]][(irep + 1):t, ], NaN(nfore - (t - irep), M[ss]))
        }
    }
                                    
    # First calculate the denominator of Equation (19) (the sum of the w's)
    sum_w2_t <- 0   
    for (k_2 in 1:nos) {       
        sum_w2_t <- sum_w2_t + w2_t[irep, k_2]
    }

    # Then calculate the DPS probabilities
    for (k_3 in 1:nos) {
        ksi_update[irep, k_3] <- (w2_t[irep, k_3] + offset) / (sum_w2_t + offset)  # this is Equation (19)
    }

    # Find best model for DMS
    [max_prob_DMS[irep, ],ss_max] <- max(ksi_update[irep, ])
    index_DMS[irep, 1] <- k_max[ss_max]

    # Now we cycled over NOM and NOS models, do DMA-over-DPS
    if (irep >= T_thres) {
        ii <- nfore
        weight_pred <- 0 * y_fore[[ss]][ii, 1:nfocus, k_max[ss], ]
        # DPS-DMA prediction
        for (ss in 1:nos) {
            temp_predict <- y_fore[[ss]][ii, 1:nfocus, k_max[ss], ] * ksi_update[irep, ss]
            weight_pred  <- weight_pred + temp_predict
            Minn_gamms[[ss]][irep - T_thres + 1, 1] <- gamma(1, k_max[[ss]])
        }

        y_t_DMA[ii, , irep - T_thres + 1] <- mean(weight_pred, 4)
        variance_DMA <- cov(t(squeeze(weight_pred)))
                   
        y_t_DMS[ii, , irep-T_thres + 1] <- mean(y_fore[[ss_max]][ii, 1:nfocus, k_max[[ss_max]], ], 4)
        variance_DMS <- cov(t(squeeze(y_fore[[ss_max]][ii, 1:nfocus, k_max[ss_max], ])))
            
            
        LOG_PL_DMS[irep-T_thres + 1, ii] <- log(mvnpdfs(t(Yraw_f[ss][ii, 1:nfocus]), t(y_t_DMS[ii, , irep - T_thres + 1]), variance_DMS) + offset)
        MAFE_DMS[irep-T_thres + 1, , ii] <- abs(Yraw_f[ss][ii, 1:nfocus] - squeeze(y_t_DMS[ii, , irep - T_thres + 1]))
        MSFE_DMS[irep-T_thres + 1, , ii] <- (Yraw_f[ss][ii, 1:nfocus] - squeeze(y_t_DMS[ii, , irep - T_thres + 1])) ^ 2
            
            
        LOG_PL_DMA[irep - T_thres + 1, ii] <- log(mvnpdfs(t(Yraw_f[ss][ii, 1:nfocus]), t(y_t_DMA(ii, , irep - T_thres + 1)), variance_DMA) + offset)
        MAFE_DMA[irep - T_thres + 1, , ii] <- abs(Yraw_f[ss][ii, 1:nfocus] - squeeze(y_t_DMA[ii, , irep - T_thres + 1]))
        MSFE_DMA[irep - T_thres + 1, , ii] <- (Yraw_f[ss][ii, 1:nfocus] - squeeze(y_t_DMA[ii, , irep - T_thres + 1])) ^ 2
        for (j in 1:nfocus) {
            logpl_DMA[irep - T_thres + 1, j, ii] <- log(mvnpdfs(t(Yraw_f[ss][ii, j]), t(y_t_DMA[ii, j, irep - T_thres + 1]), variance_DMA[j, j]) + offset)
            logpl_DMS[irep-T_thres   + 1, j, ii] <- log(mvnpdfs(t(Yraw_f[ss][ii, j]), t(y_t_DMS[ii, j, irep - T_thres + 1]), variance_DMS[j, j]) + offset)
        }
    }
}
#======================== END KALMAN FILTER ESTIMATION ========================