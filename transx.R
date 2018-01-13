# function [y]=transx(x,tcode)
# %    Transform x
# %    Return Series with same dimension and corresponding dates
# %    Missing values where not calculated
# %    -- Tcodes:
# %             1 Level
# %             2 First Difference
# %             3 Second Difference
# %             4 Log-Level
# %             5 Log-First-Difference
# %             6 Log-Second-Difference
# %             7 Detrend Log Using 1-sided HP detrending for Monthly data
# %             8 Detrend Log Using 1-sided HP detrending for Quarterly data
# %            16 Log-Second-Difference
# %            17 (1-L)(1-L^12)
# 
# %  Translated from the Gauss procs of Stock&Watson(2005),'Implications of
# %  dynamic factor models for VAR analysis'
# %  Dimitris Korobilis, June 2007
# 
# small=1.0e-040;
# relvarm=.00000075;
# relvarq=.000625;     %HP parameter
# %.00000075 for monthly data
# %.000625 for quarterly data, see Harvey/Jeager (1993), page 234 @
#     n=size(x,1);
# y=zeros(n,1);        %storage space for y
# 
# if tcode == 1
# y=x;
# elseif tcode == 2
# y(2:n)=x(2:n)-x(1:n-1);
# elseif tcode == 3
# y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
# elseif tcode == 4
# if min(x) < small
# y=NaN; 
# end
# x=log(x);
# y=x;
# elseif tcode == 5
# if min(x) < small
# y=NaN; 
# end
# x=log(x);
# y(2:n)=x(2:n)-x(1:n-1);
# elseif tcode == 6
# if min(x) < small
# y=NaN; 
# end
# x=log(x);
# y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
# elseif tcode == 7
# if min(x) < small
# y=NaN; 
# end
# x=log(x);
# [y,t1]=detrend1(x,relvarm); %#ok<NASGU>
# elseif tcode == 8
# if min(x) < small
# y=NaN; 
# end
# x=log(x);
# [y,t1]=detrend1(x,relvarq); %#ok<NASGU>
# elseif tcode == 16
# if min(x) < small
# y=NaN; 
# end
# x=log(x);
# y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
# elseif tcode == 17
# if min(x) < small
# y=NaN; 
# end
# x=log(x);
# y(14:n)=x(14:n)-x(13:n-1)-x(2:n-12)+x(1:n-13);
# else
#     y=NaN;
# end

transx <- function(x, tcode) {
    #    Transform x
    #    Return Series with same dimension and corresponding dates
    #    Missing values where not calculated
    #    -- Tcodes:
    #             1 Level
    #             2 First Difference
    #             3 Second Difference
    #             4 Log-Level
    #             5 Log-First-Difference
    #             6 Log-Second-Difference
    #             7 Detrend Log Using 1-sided HP detrending for Monthly data
    #             8 Detrend Log Using 1-sided HP detrending for Quarterly data
    #            16 Log-Second-Difference
    #            17 (1-L)(1-L^12)

    #  Translated from the Gauss procs of Stock&Watson(2005),'Implications of
    #  dynamic factor models for VAR analysis'
    #  Dimitris Korobilis, June 2007

    small   <- 1.0e-040
    relvarm <- 0.00000075
    relvarq <- 0.000625     #HP parameter
    # 0.00000075 for monthly data
    # 0.000625   for quarterly data, see Harvey/Jeager (1993), page 234 @
    x <- as.matrix(x)
    n <- nrow(x)
    y <- matrix(data = 0, nrow = nrow(x), ncol = ncol(x))       # storage space for y
 
    if (tcode == 1) {
        y <- x
        y
    } else if (tcode == 2) {
        y[2:n] <- x[2:n] - x[1:(n - 1)]
        y
    } else if (tcode == 3) {
        y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
        y
    } else if (tcode == 4) {
        if (min(x) < small) {
            y <- NaN 
        }
        x <- log(x)
        y <- x
        y
    } else if (tcode == 5) {
        if (min(x) < small) {
            y <- NaN
        }
        x <- log(x)
        y[2:n] <- x[2:n] - x[1:(n - 1)]
        y
    } else if (tcode == 6) {
        if (min(x) < small) {
            y <- NaN 
        }
        x <- log(x)
        y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
        y
    } else if (tcode == 7) {
        if (min(x) < small) {
            y <- NaN 
        }
        x <- log(x)
        y <- matrix(data = 0, nrow = nrow(x), ncol = ncol(x))
        filtro_HP <- vector(mode = 'list', length = ncol(x))
        for (i in 1:ncol(x)) {
            filtro_HP[[i]] <- hpfilter(x, freq = 1600, type = 'lambda')
            y[i] <- filtro_HP[[i]]$cycle
        }
        y
        #tl <- filtro_HP$trend
        #[y,t1]=detrend1(x, relvarm) Verificar essa função detrendl
    } else if (tcode == 8) {
        if (min(x) < small) {
                y <- NaN 
        }
        x <- log(x)
        y <- matrix(data = 0, nrow = nrow(x), ncol = ncol(x))
        filtro_HP <- vector(mode = 'list', length = ncol(x))
        for (i in 1:ncol(x)) {
            filtro_HP[[i]] <- hpfilter(x, freq = 1600, type = 'lambda')
            y[i] <- filtro_HP[[i]]$cycle
        }
        y
        #tl <- filtro_HP$trend
        #[y,t1]=detrend1(x, relvarq)
    } else if (tcode == 16) {
        if (min(x) < small) {
            y <- NaN 
        }
        x <- log(x)
        y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
        y
    } else if (tcode == 17) {
        if (min(x) < small) {
            y = NaN 
        }
        x <- log(x)
        y[14:n] <- x[14:n] - x[13:(n - 1)] - x[2:(n - 12)] + x[1:(n - 13)]
        y
    } else {
        y <- NaN
        y
    }
}

