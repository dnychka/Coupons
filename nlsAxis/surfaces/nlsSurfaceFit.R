
## 
## 13 Feb 2020
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## Find center axis of test coupons (surface data) using non-linear least squares
##
##------------------------------------------------------------------------------------------------------------------------


library(fields)
library(plyr)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("nlsFunctions.R") #lolol don't set wd in functions

setwd("D:/server-files/SURFACES")

folderList <- list.files()

for(folder in folderList){
  
  setwd(paste0("D:/server-files/SURFACES/",folder))
  
  print(folder)
  
  poreCoordinates <- readRDS(paste0(folder,".rds"))
  
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
  
  setwd("C:/Users/barna/Documents/Coupons/nlsAxis/surfaces/nlsSurfaceData")
  save(oldCoupon, nlsCoupon, nlsCoeff, file = paste0("nlsSurfaceCoupon", folder, ".rda"))
  
} # end of for loop


