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

setwd("C:/Users/barna/Documents/Coupons/datasets")

poreData <- readRDS("porosityData.rds")


## select coupon number
#n = 4
nVec <- c(4, 12, 15, 17,22,30,31,37,41,43,44,45,51,54,56,57)

medianRadius <- rep(NA, length(nVec))
j = 1

for(n in nVec){

ordered <- order(poreData[[n]]$comZ)

comX <- poreData[[n]]$comX[ordered]
comY <- poreData[[n]]$comY[ordered]
comZ <- poreData[[n]]$comZ[ordered]

deciles <- quantile(comZ, prob = seq(0, 1, length = 11), type = 5)

## subset the coupon to avoid the weld/support material remnants
good <- (comZ >= deciles[3] & comZ <= deciles[9])

poreCoordinates <- cbind( comX,
                          comY,
                          comZ)


## get centers of mass for upper and lower half of the coupon to 
## compute axisVector, which is the direction vector for 
## the initial axis estimate
upperHalf <- (comZ >= deciles[6]) & good 

lowerHalf <- (comZ <= deciles[6]) & good

axisVector <- colMeans(poreCoordinates[upperHalf,]) - colMeans(poreCoordinates[lowerHalf,])

axisVector <- axisVector / sqrt(sum(axisVector^2))


poreCoordinates <- cbind( comX[good],
                          comY[good],
                          comZ[good])

## project the centroid into the x-y plane along axisVector.
## allows us the parameterize the centroid using two parameters
## (x,y,0) rather than three (x,y,z)
centroid <- colMeans(poreCoordinates)

xyCentroid <- c( (axisVector[1]*-centroid[3])/axisVector[3] + centroid[1],
                 (axisVector[2]*-centroid[3])/axisVector[3] + centroid[2],
                 0)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("alignCoupon.R")
source("radiusAligned.R")


N <- length(poreCoordinates[,1])
radiusTarget <- rep(1000, N) # based on assumption an ideal coupon radius is 1000 micro-m

centerAxis <- nls(radiusTarget~radiusAligned(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY),
              start = list(centroidX = xyCentroid[1], centroidY = xyCentroid[2],
                           axisVectorX = axisVector[1], axisVectorY = axisVector[2]))


nlsCoeff <- coef(centerAxis)

## compare radiusFinal with radiusInitial
radiusInitial <- radiusAligned(poreCoordinates, centroid[1], centroid[2], axisVector[1], axisVector[2])

## run 'radiusAligned' with the optimal parameter values found by nls
radiusFinal <- radiusAligned(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                             nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])


plot(radiusFinal, radiusInitial, main = "cropped")
xline(median(radiusFinal), col = "cornflowerblue")

medianRadius[j] <- (median(radiusFinal))
j = j+1
}


plot(nVec, medianRadius, type = "b", pch = 20)

#saveRDS(nlsCoeff, "nlsCoeff.rds") # for surface plotting

