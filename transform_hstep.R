transform_hstep <-  function(ydata, tcode, nfore) {
    
    # RANSFORM Transform large dataset to stationarity
    # This code corrects the number of observations lost from transformations
    Yraw <- 0 * ydata
    
    for (i in 1:ncol(ydata)) {
        Yraw[ ,i] <- yfcst(ydata[ , i],tcode[i], nfore) 
    }

    Yraw <- Yraw[1:(end - nfore), ]

}
