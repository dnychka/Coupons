## 
## 29 May 2019
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## Find center axis of test coupons using non-linear least squares
##
##------------------------------------------------------------------------------------------------------------------------

library(fields)
library(plyr)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("initialParameters.R")
source("rotateCoupon.R")
source("getRadius.R")
source("newCoupon.R")
source("nlsAxisFit.R")

setwd("C:/Users/barna/Documents/Coupons/datasets")
poreData <- readRDS("porosityData.rds")


j = 1

for(n in 1:58){
  

  
##--------------------------------------------------------------------
## crop coupon
##--------------------------------------------------------------------
  
ordered <- order(poreData[[n]]$comZ)

comX <- poreData[[n]]$comX[ordered]
comY <- poreData[[n]]$comY[ordered]
comZ <- poreData[[n]]$comZ[ordered]

oldCoupon <- cbind( comX,
                    comY,
                    comZ)

cropSections <- quantile(comZ, prob = seq(0, 1, length = 11), type = 5)

## subset the coupon to avoid the weld/support material remnants
good <- (comZ >= cropSections[3] & comZ <= cropSections[9])

poreCoordinates <- cbind( comX[good],
                          comY[good],
                          comZ[good])
##--------------------------------------------------------------------

nlsObj <- nlsAxisFit(poreCoordinates, n)

nlsCoeff <- coef(nlsObj[[1]])


##--------------------------------------------------------------------
## store data
##--------------------------------------------------------------------

## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
## useful for generating surface plots and histograms for each coupon
nlsCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")
save(oldCoupon, nlsCoupon, nlsCoeff, file = paste0("nlsCoupon", n, ".rda"))

j = j+1
} # end of for loop


