


library(fields)
library(plyr)

setwd("C:/Users/barna/Documents/Coupons/nlsAxisCopy")
source("initialParameters.R")
source("alignCoupon.R")
source("radiusAligned.R")
source("newCoupon.R")




surfaceCouponTest <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxisCopy/surfaceCouponTest.rds")



# store radii for iteration comparison
radiusFirstIter <- list()
radiusLastIter <-  list()

  
  ## initial radius guess
  r = 1000 # based on ideal coupon radius of 1000 micor-meters
  ## tolerance for radius diff
  tol = 0.001
  diffMedians = 100
  
  
  ##--------------------------------------------------------------------
  ## crop coupon
  ##--------------------------------------------------------------------
  
  ordered <- order(surfaceCouponTest[,3])
  
  comX <- surfaceCouponTest[,1][ordered]
  comY <- surfaceCouponTest[,2][ordered]
  comZ <- surfaceCouponTest[,3][ordered]
  
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
    
    ## run 'radiusAligned' with the optimal parameter values found by nls
    radiusFinal <- radiusAligned(poreCoordinates, 
                                 nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                                 nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
    
    
    # store the initial median for comaprison (first iteration)
    ifelse(diffMedians == 100, radiusFirstIter[[j]] <- radiusFinal, NA) 
    
    # evaluate change in radius
    diffMedians <- abs(r - median(radiusFinal))
    r <- median(radiusFinal)

  } # end of while loop
  
  # store the final median for comparison (last iteration)
  radiusLastIter[[j]] <- radiusFinal
  
  
  ## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
  ## useful for generating surface plots and histograms for each coupon
  setwd("C:/Users/barna/Documents/Coupons/nlsAxisCopy/surfaces")
  newCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                         nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
  
  save(oldCoupon, newCoupon, nlsCoeff, file = paste0("nlsSurfaceTest.rda"))
  
  
  j = j+1



radiusFirstIter <- data.frame(lapply(radiusFirstIter, "length<-", max(lengths(radiusFirstIter))))
names(radiusFirstIter) <- 1:58


radiusLastIter <- data.frame(lapply(radiusLastIter, "length<-", max(lengths(radiusLastIter))))
names(radiusLastIter) <- 1:58

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
save(radiusFirstIter, radiusLastIter, file = "radiusIterationsCrop.rda")


