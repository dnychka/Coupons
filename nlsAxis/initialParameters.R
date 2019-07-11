##
## function to find starting values for parameters 
## centroidX, centroidY, axisVectorX, axisVectorY 
## for nls fit
##
## poreCoordinates - a tall matrix containing x, y, and z coordinates of pores
##                   (should already be cropped)   
##
## returns starting values for nls in matrix form
##
## ---------------------------------------------------------------------------


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


startingValues <- rbind(axisVector, xyCentroid)

return(startingValues)

}


