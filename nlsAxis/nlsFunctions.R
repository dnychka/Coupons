



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




##
## function to crop the coupon to give to function nlsAxisFit
##
## n - index of coupon of interest
## poreData - a tall matrix containing x, y, and z coordinates of pores
##            (taken straight from unwrapped json files)  
##
## returns poreCoordinates, a tall matrix of x, y, and z coords of pores
## from the cropped coupon
##
## ---------------------------------------------------------------------------

cropCoupon <- function(n, poreData){
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
  
  return(poreCoordinates)
}


##
## This function takes nls-aligned coordinates (ex, nlsCoupon) and rotates 
## them back to their original tilt--important for generating relevant
## bootstrapped data sets
##
## poreCoordinates - a tall matrix containing x, y, and z coordinates of pores
## centroidX, centroidY - x and y coordinates of centroid, respectively
## axisVectorX, axisVectorY - x and y directions of axis direction vector
##
## returns coordinates of misaligned coupon
##
## ---------------------------------------------------------------------------


getOldCoupon <- function(newCoordinates, centroidX, centroidY, axisVectorX, axisVectorY){
  
  xyCentroid <- c(centroidX,centroidY,0)
  
  axisVector <- c(axisVectorX,axisVectorY,1)
  
  V1 <- axisVector*sign(axisVector[3])
  
  U1<- V1/ sqrt(sum( V1^2))
  
  I3<- diag( 1,3)
  
  V<- qr.qy(qr( cbind( U1, c(-1,0,0), c( 0,-1,0) )),I3) 
  
  V<- V[,c(2,3,1)]
  
  newXYZ <- t(newCoordinates)
  
  oldXYZ <- solve(t(V)) %*% newXYZ + xyCentroid
  
  oldXYZ <- t(oldXYZ)
  
  
  return(oldXYZ)
  
}


##
## This is just the radiusAligned function but altered sightly
## to return the coordinates of the "new" rotated coupon instead of the pore
## radius. For use in the surface plots
##
## poreCoordinates - a tall matrix containing x, y, and z coordinates of pores
## centroidX, centroidY - x and y coordinates of centroid, respectively
## axisVectorX, axisVectorY - x and y directions of axis direction vector
##
## returns coordinates of aligned coupon
##
## ---------------------------------------------------------------------------

newCoupon <- function(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY){
  
  centroidGuess <- c(centroidX,centroidY,0)
  
  axisGuess <- c(axisVectorX,axisVectorY,1)
  
  newCoupon <- alignCoupon(poreCoordinates, centroidGuess, axisGuess)
  
  return(newCoupon)
}
