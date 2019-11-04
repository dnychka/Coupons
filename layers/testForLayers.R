##
##
##
## 
## Calculate test statistic for each population of 
##  coupons
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
## 0 degree coupons, layered
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/spacing40to60")

inputFiles <- list.files(full.names = TRUE)

N <- length(inputFiles) #number of coupons in the sample size

TSZlayers <- rep(NA, N) #test statistic (TS), zero degree, layers

i = 1 #dummy index to increment TS

for(coupon in inputFiles){
  
  load(coupon)
    
  # calculate the test statistic
    
  TSZlayers[i] <- getTestStatistic(nlsCoupon)
  i = i + 1
}


## ------------------------------
## 0 degree coupons,not layered
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/noLayers")

inputFiles <- list.files(full.names = TRUE)

N <- length(inputFiles) #number of coupons in the sample size

TSZnolayers <- rep(NA, N) #test statistic (TS), zero degree, not layered

i = 1 #dummy index to increment TS

for(coupon in inputFiles){
  
  load(coupon)
  
  # calculate the test statistic
  
  TSZnolayers[i] <- getTestStatistic(nlsCoupon)
  i = i + 1
}


## ------------------------------
## 0 degree coupons, real data
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

inputFiles <- list.files(full.names = TRUE)

zeros <- which(couponCov$polarAngle == 0)

N <- length(inputFiles[zeros]) #number of coupons in the sample size

TSZreal <- rep(NA, N) #test statistic (TS), zero degree, not layered

i = 1 #dummy index to increment TS

for(coupon in inputFiles[zeros]){
  
  load(coupon)
  
  # calculate the test statistic
  
  TSZreal[i] <- getTestStatistic(nlsCoupon)
  i = i + 1
}


## -------------------
## save all zero degree TS to data folder
## -------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData")

saveRDS(TSZlayers, "TSZlayers.rds")
saveRDS(TSZnolayers, "TSZnoLayers.rds")
saveRDS(TSZreal, "TSZreal.rds")






## ---------------------------------------------------------------------------------
## FORTY-FIVE DEGREE COUPONS
## ---------------------------------------------------------------------------------


## ------------------------------
## 45 degree coupons, layered
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/spacing40to60")

inputFiles <- list.files(full.names = TRUE)

N <- length(inputFiles) #number of coupons in the sample size

angleSeq <- seq(0, 2*pi, length.out = 100) #increment by which we change rotation

tempTS <- rep(NA, 100) #store TS for each rotation

TSFlayers <- rep(NA, N) #test statistic (TS), forty-five degree, layered

i = 1 #dummy index to increment TS

for(coupon in inputFiles){
  
  load(coupon)
  
  print(i)
  
  m = 1
  for(w in angleSeq){
    
    rotatedCoupon <- rotateFortyFive(w, nlsCoupon) #center and rotate the coupon
    
    tempTS[m] <- getTestStatistic(rotatedCoupon) #calc TS for rotated coupon
    
    m=m+1
  }
  
  TSFlayers[i] <- max(tempTS)
  i = i + 1
}



## ------------------------------
## 45 degree coupons, not layered
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/noLayers")

inputFiles <- list.files(full.names = TRUE)

N <- length(inputFiles) #number of coupons in the sample size

angleSeq <- seq(0, 2*pi, length.out = 100) #increment by which we change rotation

tempTS <- rep(NA, 100) #store TS for each rotation

TSFnolayers <- rep(NA, N) #test statistic (TS), forty-five degree, not layered

i = 1 #dummy index to increment TS

for(coupon in inputFiles){
  
  load(coupon)
  
  print(i)
  
  m = 1
  for(w in angleSeq){
    
    rotatedCoupon <- rotateFortyFive(w, nlsCoupon) #center and rotate the coupon
    
    tempTS[m] <- getTestStatistic(rotatedCoupon) #calc TS for rotated coupon
    
    m=m+1
  }
  
  TSFnolayers[i] <- max(tempTS)
  i = i + 1
}


## ------------------------------
## 45 degree coupons, real data
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

inputFiles <- list.files(full.names = TRUE)

fortyFive <- which(couponCov$polarAngle==45)

N <- length(inputFiles[fortyFive]) #number of coupons in the sample size

angleSeq <- seq(0, 2*pi, length.out = 100) #increment by which we change rotation

tempTS <- rep(NA, 100) #store TS for each rotation

TSFreal <- rep(NA, N) #test statistic (TS), forty-five degree, not layered

i = 1 #dummy index to increment TS

for(coupon in inputFiles[fortyFive]){
  
  load(coupon)
  
  print(i)
  
  m = 1
  for(w in angleSeq){
    
    rotatedCoupon <- rotateFortyFive(w, nlsCoupon) #center and rotate the coupon
    
    tempTS[m] <- getTestStatistic(rotatedCoupon) #calc TS for rotated coupon
    
    m=m+1
  }
  
  TSFreal[i] <- max(tempTS)
  i = i + 1
}



## -------------------
## save all forty-five degree TS to data folder
## -------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData")

saveRDS(TSFlayers, "TSFlayers.rds")
saveRDS(TSFnolayers, "TSFnoLayers.rds")
saveRDS(TSFreal, "TSFreal.rds")


