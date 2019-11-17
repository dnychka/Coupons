



getDistance <- function(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY, r){
  
  centroidGuess <- c(centroidX,centroidY,0)
  
  axisGuess <- c(axisVectorX,axisVectorY,1)
  
  newCoupon <- alignCoupon(poreCoordinates, centroidGuess, axisGuess)
  
  newRadius <- sqrt(newCoupon[,1]^2+newCoupon[,2]^2)
  
  dist <- (newRadius-rep(r,length(newRadius)))^2
  
  return(dist)
}


getInitialParameters <- function(poreCoordinates){
  
  ## get centers of mass for upper and lower half of the coupon to 
  ## compute axisVector, which is the direction vector for 
  ## the initial axis estimate
  deciles <- quantile(poreCoordinates[,3], prob = seq(0, 1, length = 11), type = 5)
  
  upperHalf <- poreCoordinates[,3] >= deciles[6]
  
  lowerHalf <- poreCoordinates[,3] <= deciles[6]
  
  axisVector <- colMeans(poreCoordinates[upperHalf,]) - colMeans(poreCoordinates[lowerHalf,])
  
  axisVector <- axisVector / sqrt(sum(axisVector^2)) # make it unit length
  
  
  ## project the centroid into the x-y plane along axisVector.
  ## allows us the parameterize the centroid using two parameters
  ## (x,y,0) rather than three (x,y,z)
  centroid <- colMeans(poreCoordinates)
  
  xyCentroid <- c( (axisVector[1]*-centroid[3])/axisVector[3] + centroid[1],
                   (axisVector[2]*-centroid[3])/axisVector[3] + centroid[2],
                   0)
  
  r = 1000 #hard coded for testing
  
  startingValues <- c(as.vector(axisVector), as.vector(xyCentroid), r)
  
  return(startingValues)
  
}

alignCoupon <- function(poreCoordinates, xyCentroid, axisVector){
  
  V1 <- axisVector*sign(axisVector[3])
  
  U1<- V1/ sqrt(sum( V1^2))
  
  I3<- diag( 1,3)
  
  V<- qr.qy(qr( cbind( U1, c(-1,0,0), c( 0,-1,0) )),I3) 
  
  V<- V[,c(2,3,1)]
  
  newXYZ <- t(V) %*% (t(poreCoordinates)-xyCentroid)
  
  newXYZ <- t(newXYZ)
  
  return(newXYZ)
  
}





nlsAxisFit <- function(poreCoordinates){  
  
  startValues <- getInitialParameters(poreCoordinates)
  
  N <- length(poreCoordinates[,1])
  
  
  distTarget <- rep(0,N)
  
  centerAxis <- nls(distTarget~getDistance(poreCoordinates, 
                                           centroidX, centroidY, 
                                           axisVectorX, axisVectorY, r),
                    start = list(centroidX = startValues[4], 
                                 centroidY = startValues[5],
                                 axisVectorX = startValues[1], 
                                 axisVectorY = startValues[2],
                                 r = startValues[7]),
                    control =  nls.control(minFactor = 1/10000000000))
  
  
  nlsCoeff <- coef(centerAxis)
  nlsCoeff
  
  
  return(centerAxis)
  
}



setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("newCoupon.R")
source("oldCoupon.R")
