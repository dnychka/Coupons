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
source("threshold.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData")

## ------------------------------------------------------

angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta

harmonicSignals <- matrix(NA, nrow = length(angleSeq), ncol = length(layerSpacing))


layerSpacing <- c(40,45,50,55,60)

i=1
for(n in layerSpacing){

    print(n)
    
    load(paste0(n,"spacingSyn45.rda"))
    
    # center the coupon
    centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                          scale(newCoupon[,2], center = TRUE, scale = FALSE),
                          newCoupon[,3])
    
    
    histSave <- FFrotate(angleSeq, centerCoupon)
    
    harmonicSignals[,i] <-  makePeriodogram(angleSeq, histSave)
    i=i+1
}

boxplot(apply(harmonicSignals,2,max), main = "relative peak strength, synthetic layers")


