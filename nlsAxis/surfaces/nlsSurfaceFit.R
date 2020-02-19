
## 
## 13 Feb 2020
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## Find center axis of test coupons (surface data) using non-linear least squares
##
## NOTE: if nls chokes just run the for loop again; it probably just hit a local minimum
##
##------------------------------------------------------------------------------------------------------------------------


library(fields)
library(plyr)

setwd("C:/Users/barna/Documents/Coupons/datasets")
nPores <- readRDS("numberOfPoresPerCoupon.rds")
wantThese <- readRDS("0and45degreeNames.rds")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("nlsFunctions.R") #lolol don't set wd in functions

## ----------------------------------------------------------------------
## get the relevant coupon names
## --------------------------------------------------------------------

nPores <- nPores[names(nPores) %in% names(wantThese)]

missing <- which(names(nPores) == "H21" | names(nPores) == "O9" | names(nPores) == "T3")

nPores <- nPores[-missing] #remove the sadly missing coupons

nPores <- nPores[order(names(nPores))] #order to match alphabetic file structure

rm(wantThese, missing)

## --------------------------------------------------------------------


setwd("D:/server-files/SURFACES")

folderList <- list.files()

poreInd <- 1 #dummy index to access nPores vector

for(folder in folderList){
  
  setwd(paste0("D:/server-files/SURFACES/",folder))
  
  print(folder)
  
  poreCoordinates <- readRDS(paste0(folder,".rds"))
  
  ##--------------------------------------------------------------------
  
  nlsObj <- nlsAxisFit(poreCoordinates)
  
  nlsCoeffBig <- coef(nlsObj)
  
  
  ##--------------------------------------------------------------------
  ## store first round of data
  ##--------------------------------------------------------------------
  
  ## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
  ## useful for generating surface plots and histograms for each coupon
  bigCoupon <- newCoupon(poreCoordinates, nlsCoeffBig["centroidX"], nlsCoeffBig["centroidY"], 
                         nlsCoeffBig["axisVectorX"], nlsCoeffBig["axisVectorY"])
  
  oldCoupon <- poreCoordinates
  
  
  ##--------------------------------------------------------------------
  ## downsample surface proportional to num of pores and run nls
  ##--------------------------------------------------------------------
  
  poreCoordinates <- poreCoordinates[sample(poreCoordinates, nPores[poreInd], replace = F),]
  
  ##--------------------------------------------------------------------
  
  nlsObj <- nlsAxisFit(poreCoordinates)
  
  nlsCoeff <- coef(nlsObj)
  
  
  ##--------------------------------------------------------------------
  ## store second round of data
  ##--------------------------------------------------------------------
  
  nlsCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                         nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
  
  ##--------------------------------------------------------------------
  ## save it
  ##--------------------------------------------------------------------
  
  setwd("C:/Users/barna/Documents/Coupons/nlsAxis/surfaces/nlsSurfaceData")
  save(bigCoupon, oldCoupon, nlsCoupon, nlsCoeff, nlsCoeffBig, 
       file = paste0("nlsSurfaceCoupon", folder, ".rda"))
  
  poreInd <- poreInd + 1
  
} # end of for loop


