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