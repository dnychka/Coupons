##
##
## 
## Function to calculate test statistic for each 
## population of coupons
##
## INPUTS
## nlsCoupon - a matrix for an oriented, centered coupon where 
##             the columns are the x,y,z coords of the pores 
##
## OUTPUTS
## testStatistic - the numerical value for the test stat
##
## Uses chi-squared distribution and mvspec function
## to generate confidence intervals for each periodogram
## ordinate (different from the original moving average
## method used up until 09/27/2019)
## ------------------------------------------------------


getTestStatistic <- function(nlsCoupon){ 

    testStatistic <- NA

    
    #create series data (count values for hist bins)
    h <- hist(~nlsCoupon[,3], w=5, plot = FALSE) 
    
    counts <- h$counts
    coords <- h$mids
    
    
    # detrend the series using a spline 
    xGrid <- seq(min(coords), max(coords), length.out = length(coords))
    highpass <- splint(coords, counts, xGrid, lambda = 50) #lambda set to be very smooth
    
    detrendedData <- counts - highpass
    
    
    # calculate the dft and compute smoothed periodogram using mvspec(.)
    P <- mvspec(detrendedData, kernel('daniell',4), log='no', plot = FALSE)
    
    # calculate lower bounds of 95% confidence intervals for each
    # periodogram ordinate
    df = P$df
    L = qchisq(.975, df) 
    U = qchisq(.025, df)
    n = length(P$freq)
    
    lowerBound <- rep(NA, n)
    upperBound <- rep(NA, n)
    for(i in 1:n){
      lowerBound[i] <- df*P$spec[i]/L
      upperBound[i] <- df*P$spec[i]/U
    }
    
    
    # establish baseline for periodogram using linear model
    baseline <- lm(P$spec ~ P$freq)
    
    
    # find peaks in periodogram
    peaks <- which(lowerBound > baseline$fitted.values)
    
    
    # calculate test statistic
    testStatistic <- sum(P$spec[peaks])/sum(P$spec)

    return(testStatistic)

}
