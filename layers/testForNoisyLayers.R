##
##
##
## 
## Calculate TS for noisy synthetic coupons
##
## Uses chi-squared distribution and mvspec function
## to generate confidence intervals for each periodogram
## ordinate (different from the original moving average
## method used up until 09/27/2019)
## ------------------------------------------------------

library(astsa)
library(fields)
library(FSA)

# load required function files
setwd("C:/Users/barna/Documents/Coupons/layers/layerFunctions")
load("couponCov.rda")
source("calculateTestStatistic.R")
source("rotations.R")


## ---------------------------------------------------------------------------------
## ZERO DEGREE COUPONS
## ---------------------------------------------------------------------------------

## ------------------------------
## 0 degree coupons, noise 
## ------------------------------

noiseAmt <- matrix(nrow = 100, ncol = 10)

for(j in 1:10){
  
  setwd(paste0("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/noisySpacing40to60/noise0.", j))
  
  inputFiles <- list.files(full.names = TRUE)
  
  i = 1 #dummy index to increment TS
  
  for(coupon in inputFiles){
    
    load(coupon)
    
    # calculate the test statistic
    
    noiseAmt[i,j] <- getTestStatistic(nlsCoupon)
    i = i + 1
  }
}


boxplot(noiseAmt)



## ----------------------------------------
## save data
## ----------------------------------------
setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/noisySpacing40to60")
saveRDS(noiseAmt, "TSZnoise.rds")
