## 
## 28 Jan 2020
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## Find center axis of test coupons using non-linear least squares
##
##------------------------------------------------------------------------------------------------------------------------

library(fields)
library(plyr)

setwd("C:/Users/barna/Documents/Coupons/datasets")
poreData <- readRDS("porosityData.rds")


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty")
source("nlsFunctionsForSimulation.R") #lolol don't set wd in functions


for(n in 1:58){
  
print(n)
  
##--------------------------------------------------------------------
## crop coupon
##--------------------------------------------------------------------

poreCoordinates <- cropCoupon(n, poreData)

##--------------------------------------------------------------------

nlsObj <- nlsAxisFit(poreCoordinates)

nlsCoeff <- coef(nlsObj)


##--------------------------------------------------------------------
## store data
##--------------------------------------------------------------------

## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
## useful for generating surface plots and histograms for each coupon
nlsCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])

oldCoupon <- poreCoordinates

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")
save(oldCoupon, nlsCoupon, nlsCoeff, file = paste0("nlsCoupon", n, ".rda"))

} # end of for loop


