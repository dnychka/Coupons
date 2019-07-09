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