# function [Yraw,yearlab] = transform(ydata,tcode,yearlab)
# %TRANSFORM Transform large dataset to stationarity
# %This code corrects the number of observations lost from transformations
# Yraw = 0*ydata;
# for i=1:size(ydata,2)
# Yraw(:,i) = transx(ydata(:,i),tcode(i)); 
# end
# end

transform <- function(ydata, tcode, yearlab) {
    # TRANSFORM Transform large dataset to stationarity
    # This code corrects the number of observations lost from transformations
    Yraw <- 0 * ydata
    for (i in 1:ncol(ydata)) {
        Yraw[ , i] <- transx(ydata[ , i], tcode[i])
    }
}

