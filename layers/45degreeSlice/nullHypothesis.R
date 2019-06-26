## 6/25/2019
##
##
##
## working document for testing null hypothesis 
## distribution
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("threshold.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/nullHyp")

## ------------------------------------------------------



angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta

numSims <- 50
relPeakStrength <- matrix(NA, nrow = length(angleSeq), ncol = numSims)
relSignals <- NULL


layerSpace <- c(55,45,35,25,15,5)

for(i in layerSpace){

  for(n in 1:numSims){
    
    print(n)
  
    load(paste0(i,"spacingRep",n,".rda"))
    
    # center the coupon
    centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                          scale(newCoupon[,2], center = TRUE, scale = FALSE),
                          newCoupon[,3])
    
  
    histSave <- FFrotate(angleSeq, centerCoupon)
  
    relPeakStrength[,n] <-  makePeriodogram(angleSeq, histSave)
    
  }
  
  relPeakStrength <- rbind(rep(i, length(numSims)), relPeakStrength) #name the columns for fbplot
  
  relSignals <- cbind(relSignals, relPeakStrength)
  
  relPeakStrength <- relPeakStrength[-1,]
 
} 


bplot(t(relSignals[-1,]), by = relSignals[1,])

bplot(apply(relSignals[-1,], 2, max), by = relSignals[1,]) #showing its easier to pick up layer spacing in 25 than 10....

bplot(apply(relSignals[-1,], 2, mean), by = relSignals[1,])

saveRDS(relSignals, "relSignals.rda")





