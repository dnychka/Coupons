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

setwd("C:/Users/barna/Documents/Coupons/datasets")
poreData <- readRDS("porosityData.rds")

#Ninety <- c(2,  3,  7,  9, 13, 14, 16, 19, 21, 23, 24, 32, 33, 35, 38, 47, 48, 50)

# store radii for iteration comparison
radiusFirstIter <- list()
radiusLastIter <-  list()
SEfirst <- list()
SEfinal <- list()


j = 1

for(n in 1:58){
  
  ## initial radius guess
  r = 1000 # based on ideal coupon radius of 1000 micor-meters
  ## tolerance for radius diff
  tol = 0.001
  diffMedians = 100

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

oldCoupon <- poreCoordinates # to save the orginial coords for 3d plotting

## get centers of mass for upper and lower half of the coupon to 
## compute axisVector, which is the direction vector for 
## the initial axis estimate
upperHalf <- (comZ >= deciles[6]) & good 

lowerHalf <- (comZ <= deciles[6]) & good

axisVector <- colMeans(poreCoordinates[upperHalf,]) - colMeans(poreCoordinates[lowerHalf,])

axisVector <- axisVector / sqrt(sum(axisVector^2)) # make it unit length


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

while (diffMedians > tol){
radiusTarget <- rep(r, N) 

centerAxis <- nls(radiusTarget~radiusAligned(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY),
              start = list(centroidX = xyCentroid[1], centroidY = xyCentroid[2],
                           axisVectorX = axisVector[1], axisVectorY = axisVector[2]))


nlsCoeff <- coef(centerAxis)

## compare radiusFinal with radiusInitial
radiusInitial <- radiusAligned(poreCoordinates, centroid[1], centroid[2], axisVector[1], axisVector[2])

## run 'radiusAligned' with the optimal parameter values found by nls
radiusFinal <- radiusAligned(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                             nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])


# store the initial median for comaprison (first iteration)
ifelse(diffMedians == 100, radiusFirstIter[[j]] <- radiusFinal, NA) 
ifelse(diffMedians == 100, SEfirst[[j]] <- summary(centerAxis)$sigma, NA) 


diffMedians <- abs(r - median(radiusFinal))

r <- median(radiusFinal)

print(j)

if(j==23){break}
} # end of while loop

# store the final median for comparison (last iteration)
radiusLastIter[[j]] <- radiusFinal
SEfinal[[j]] <- summary(centerAxis)$sigma
j = j+1

## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
## useful for generating surface plots and histograms for each coupon
source("newCoupon.R")
setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData/cropped")
newCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])

save(oldCoupon, newCoupon, nlsCoeff, file = paste0("nlsCoupon", n, ".rda"))

} # end of for loop


radiusFirstIter <- data.frame(lapply(radiusFirstIter, "length<-", max(lengths(radiusFirstIter))))
names(radiusFirstIter) <- 1:58


radiusLastIter <- data.frame(lapply(radiusLastIter, "length<-", max(lengths(radiusLastIter))))
names(radiusLastIter) <- 1:58

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
save(radiusFirstIter, radiusLastIter, file = "radiusIterationsCrop.rda")


