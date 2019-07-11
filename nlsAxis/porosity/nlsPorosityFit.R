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

setwd("C:/Users/barna/Documents/Coupons/nlsAxisCopy")
source("initialParameters.R")
source("rotateCoupon.R")
source("getRadius.R")
source("newCoupon.R")

setwd("C:/Users/barna/Documents/Coupons/datasets")
poreData <- readRDS("porosityData.rds")

# store radii for iteration comparison
radiusFirstIter <- list()
radiusLastIter <-  list()


j = 1

for(n in 1:5){
  
  ## initial radius guess
  r = 1000 # based on ideal coupon radius of 1000 micor-meters
  ## tolerance for radius diff
  tol = 0.001
  diffMedians = 100

  
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



startValues <- getInitialParameters(poreCoordinates)

N <- length(poreCoordinates[,1])


while (diffMedians > tol){
      radiusTarget <- rep(r, N) 
      
      centerAxis <- nls(radiusTarget~radiusAligned(poreCoordinates, 
                                                   centroidX, centroidY, 
                                                   axisVectorX, axisVectorY),
                                      start = list(centroidX = startValues[2,1], 
                                                   centroidY = startValues[2,2],
                                                   axisVectorX = startValues[1,1], 
                                                   axisVectorY = startValues[1,2]))
      
      
      nlsCoeff <- coef(centerAxis)
      
      ## run 'getRadius' with the optimal parameter values found by nls
      radiusFinal <- radiusAligned(poreCoordinates, 
                                   nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                                   nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
      
      
      # store the initial median for comaprison (first iteration)
      ifelse(diffMedians == 100, radiusFirstIter[[j]] <- radiusFinal, NA) 
      
      # update radius and compare to last iteration
      diffMedians <- abs(r - median(radiusFinal))
      r <- median(radiusFinal)
      
      print(j)
      
      if(j==23){break}
} # end of while loop






##--------------------------------------------------------------------
## store data
##--------------------------------------------------------------------

# store the final median for comparison (last iteration)
radiusLastIter[[j]] <- radiusFinal


## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
## useful for generating surface plots and histograms for each coupon
upCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])

setwd("C:/Users/barna/Documents/Coupons/nlsAxisCopy/porosity/couponCaseStudies/caseStudyData/cropped")
save(oldCoupon, upCoupon, nlsCoeff, file = paste0("nlsCoupon", n, ".rda"))


j = j+1
} # end of for loop


radiusFirstIter <- data.frame(lapply(radiusFirstIter, "length<-", max(lengths(radiusFirstIter))))
names(radiusFirstIter) <- 1:58


radiusLastIter <- data.frame(lapply(radiusLastIter, "length<-", max(lengths(radiusLastIter))))
names(radiusLastIter) <- 1:58

setwd("C:/Users/barna/Documents/Coupons/nlsAxisCopy/couponCaseStudies/caseStudyData/cropped")
save(radiusFirstIter, radiusLastIter, file = "radiusIterations.rda")

##--------------------------------------------------------------------
