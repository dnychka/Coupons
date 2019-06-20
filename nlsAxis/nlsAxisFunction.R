## 
## 17 June 2019
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## nlsCouponFit: cropped vs uncropped?
##
##------------------------------------------------------------------------------------------------------------------------

nlsAxisFit <- function(poreCoordinates, xyCentroid, axisVector,
                       diffMedians, tol, r){
  
  centroidX <- xyCentroid[1]
  centroidY <- xyCentroid[2]
  
  axisVectorX <- axisVector[1]
  axisVectorY <- axisVector[2]
  
setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("alignCoupon.R")
source("radiusAligned.R")
  
summaries <- list()
i = 1
N <- length(poreCoordinates[,1])
  
  while (diffMedians > tol){
    radiusTarget <- rep(r, N) 
    
    centerAxis <- nls(radiusTarget~radiusAligned(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY),
                      start = list(centroidX = xyCentroid[1], centroidY = xyCentroid[2],
                                   axisVectorX = axisVector[1], axisVectorY = axisVector[2]),
                      control = nls.control(maxiter=150, minFactor = 1/5000))
    nlsCoeff <- coef(centerAxis)
  
    radiusFinal <- radiusAligned(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                                 nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
    
    diffMedians <- abs(r - median(radiusFinal))
    
    r <- median(radiusFinal)
  
  }
  
summaries <- summary(centerAxis)$sigma

return(summaries)

}




