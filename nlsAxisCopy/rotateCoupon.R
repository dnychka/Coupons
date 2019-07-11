##
## rotate the coupon
##
## poreCoordinates - a tall matrix containing x, y, and z coordinates of pores
## xyCentroid - centroid coordinates in the xy plane: (x,y,0)
## axisVector - direction vector for initial axis guess
##
## returns new xyz pore coordinates 
##
## ----------------------------------------------------------------------------

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
