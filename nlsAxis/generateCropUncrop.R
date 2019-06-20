##
##
##
## create the cropped and uncropped versions of coupons and associated
## coeffificnets to feed into nls
## ---------------------------------------------------------------------------------


library(fields)
library(MASS)

setwd("C:/Users/barna/Documents/Coupons/datasets")
poreData <- readRDS("porosityData.rds")

Zero <- which(couponCov$polarAngle == 0)
Zero <- Zero[-9]

fortyFive <- which(couponCov$polarAngle==45)

Ninety <- which(couponCov$polarAngle==90)

for(n in c(1:22,24:40,42:58)){
## -----------------------------
## initial radius guess
r = 1000 # based on ideal coupon radius of 1000 micor-meters
## tolerance for radius diff
tol = 0.001
diffMedians = 100
## -----------------------------

ordered <- order(poreData[[n]]$comZ)

comX <- poreData[[n]]$comX[ordered]
comY <- poreData[[n]]$comY[ordered]
comZ <- poreData[[n]]$comZ[ordered]

deciles <- quantile(comZ, prob = seq(0, 1, length = 11), type = 5)

poreCoordinates <- cbind( comX,
                          comY,
                          comZ)


uncropPoreCoords <- poreCoordinates # save the uncropped coupon


## -------------------------------
## get coefficient values for uncropped

upperHalf <- (comZ >= deciles[6])

lowerHalf <- (comZ <= deciles[6])

axisVectorU <- colMeans(poreCoordinates[upperHalf,]) - colMeans(poreCoordinates[lowerHalf,])

axisVectorU <- axisVectorU / sqrt(sum(axisVectorU^2)) # make it unit length

centroid <- colMeans(uncropPoreCoords)

xyCentroidU <- c( (axisVectorU[1]*-centroid[3])/axisVectorU[3] + centroid[1],
                 (axisVectorU[2]*-centroid[3])/axisVectorU[3] + centroid[2],
                 0)


## -------------------------------
## get coefficient values for cropped


good <- (comZ >= deciles[3] & comZ <= deciles[9])

cropPoreCoords <-  cbind( comX[good],
                          comY[good],
                          comZ[good])

upperHalf <- (comZ >= deciles[6]) & good 

lowerHalf <- (comZ <= deciles[6]) & good

axisVectorC <- colMeans(poreCoordinates[upperHalf,]) - colMeans(poreCoordinates[lowerHalf,])

axisVectorC <- axisVectorC / sqrt(sum(axisVectorC^2)) # make it unit length

centroid <- colMeans(cropPoreCoords)

xyCentroidC <- c( (axisVectorC[1]*-centroid[3])/axisVectorC[3] + centroid[1],
                 (axisVectorC[2]*-centroid[3])/axisVectorC[3] + centroid[2],
                 0)



## --------------------------------------------
## run nlsAxisFunction on the crop uncrop

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("nlsAxisFunction.R")

uncrop <- nlsAxisFit(uncropPoreCoords, xyCentroidU, axisVectorU, diffMedians, tol, r)

crop <- nlsAxisFit(cropPoreCoords, xyCentroidC, axisVectorC, diffMedians, tol, r)

print(n)
print(uncrop)
print(crop)

# library(DescTools)
# 
# print(n)

# for(i in 1:4){
# print(uncrop[i,] %overlaps% crop[i,])
# }
}










