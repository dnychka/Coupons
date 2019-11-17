##
##
##
##
##
##
##
## Bayseian uncertainty sampling
## -------------------------------------------------------

library(fields)
library(MASS)
library(FSA)
library(rgl)
library(astsa)


setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("initialParameters.R")
source("rotateCoupon.R")
source("getRadius.R")
source("newCoupon.R")
source("cropCoupon.R")
source("nlsAxisFit.R")


setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/spacing40to60")

inputFiles <- list.files(full.names = TRUE)

##--------------------------------------------------------------------
## look at RSE for the real coupons. Effect of interior crop on RSE
##--------------------------------------------------------------------

TSZnlsObj <- list()
rse <- rep(NA, length(inputFiles))


n=1
for(coupon in inputFiles){
  
  load(coupon)
  
  rse[n] <- summary(centerAxis)$sigma
  
  
  TSZnlsObj[[n]] <- centerAxis

n=n+1
  
}
