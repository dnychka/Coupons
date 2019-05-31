##
## function for use in nls to optimize parameters 
## centroidX, centroidY, axisVectorX, axisVectorY to best
## find the center axis. 
##
## poreCoordinates - a tall matrix containing x, y, and z coordinates of pores
## centroidX, centroidY - x and y coordinates of centroid, respectively
## axisVectorX, axisVectorY - x and y directions of axis direction vector
##
## returns radius of aligned coupon from center axis
##
## ---------------------------------------------------------------------------

radiusAligned <- function(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY){
  
  centroidGuess <- c(centroidX,centroidY,0)
  
  axisGuess <- c(axisVectorX,axisVectorY,1)
  
  newCoupon <- alignCoupon(poreCoordinates, centroidGuess, axisGuess)
  
  saveRDS(newCoupon, "newCoupon.rds") # for surface example plotting
  
  radius <- sqrt(newCoupon[,1]^2+newCoupon[,2]^2)
  
  return(radius)
}