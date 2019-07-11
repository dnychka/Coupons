## 6/25/2019
##
##
##
## test out our summary statistics in the synthetic layered
## coupons
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("findPeaks.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/spacingLessThan30")

## ------------------------------------------------------

angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta

layerSpacing <- seq(40,60,length.out = 150)

layerSpacingFine <- seq(10,30,length.out=50)

harmonicSignals <- matrix(NA, nrow = length(angleSeq), ncol = length(layerSpacingFine))

i=1
for(n in layerSpacingFine){

    print(n)
    
    load(paste0(round(n,3),"spacingSyn45.rda"))
    
    # center the coupon
    centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                          scale(newCoupon[,2], center = TRUE, scale = FALSE),
                          newCoupon[,3])
    
    
    histSave <- FFrotate(angleSeq, centerCoupon)
    
    harmonicSignals[,n] <-  makePeriodogram(angleSeq, histSave)
    i=i+1
}

boxplot(apply(harmonicSignals,2,max), main = "relative peak strength, synthetic layers")

saveRDS(harmonicSignals, "45sFineLayeredSignals.rds")
