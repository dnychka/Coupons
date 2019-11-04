## 7/04/2019
##
##
##
## construct a test case for the 0 polar angle coupons
## ------------------------------------------------------


library(conicfit)
library(rgl)

# numLayer <- seq(10,30,length.out = 50)
# 
noiseAmt <- seq(0.1,1, by=0.1)

numLayer <- seq(40,60,length.out=100)

which(round(numLayer,3)==57.576)

for(l in numLayer[1:100]){
  
  z = seq(0,4000,by = l)
  
noise <- runif(length(z),-l*noiseAmt[10],l*noiseAmt[10])

  z = z+noise

y = rep(0, length(z))
x = y


a = 814  #radius

pts <- matrix(NA, nrow = 25, ncol = 3)
pts[,1:2] <- calculateEllipse(x[1],y[1],a,a, steps = 25, randomDist = TRUE)
pts[,3] <- rep(z[1], 25)

ellipseStack <- pts

for(i in 2:length(x)){
  
  pts <- matrix(NA, nrow = 25, ncol = 3)
  pts[,1:2] <- calculateEllipse(x[i],y[i],a,a, steps = 25, randomDist = TRUE)
  pts[,3] <- rep(z[i], 25)
  
  ellipseStack <- rbind(ellipseStack, pts)
  
}

#plot3d(ellipseStack[,1], ellipseStack[,2], ellipseStack[,3], type = "s", size = 0.45)

# open3d()
# plot3d(ellipseStack[,1], ellipseStack[,2], ellipseStack[,3], type = "s", size = 0.45,
#        zlab = "", xlab = "", ylab = "",
#        axes = FALSE)
# box3d(edges = "bbox", tick = FALSE, box = TRUE, col = "black")


#plot(ellipseStack[,1], ellipseStack[,2])

aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}

theta <- -runif(1,2,14)*pi/180 #rotate the coupon so it's almost upright to match the 
#usual orientation of the actual coupons. Works better with nls too
theta = -8*pi/180

Ry <- aboutY(theta)

turnedStack <- ellipseStack %*% Ry 

# open3d()
# plot3d(turnedStack[,1], turnedStack[,2], turnedStack[,3], type = "s", size = 0.45)

## ----------------------------------------------------------------------
## modified version of nlsCouponFit.R
## removed while loop--since the radius of the ellipse is set,
## nls finds the optimal orientation quite fast and there's zero
## distance between the ellipse points and the target radius,
## casuing the step size to be reduced past minimum threshold
## ----------------------------------------------------------------------

deciles <- quantile(turnedStack[,3], prob = seq(0, 1, length = 11), type = 5)

comX <- turnedStack[,1]
comY <- turnedStack[,2]
comZ <- turnedStack[,3]

## subset the coupon to avoid the weld/support material remnants

good <- (comZ >= deciles[3] & comZ <= deciles[9])

poreCoordinates <- cbind( comX,
                          comY,
                          comZ)

oldCoupon <- poreCoordinates # to save the orginial coords for 3d plotting

## get centers of mass for upper and lower half of the coupon to 
## compute axisVector, which is the direction vector for 
## the initial axis estimate
upperHalf <- (comZ >= deciles[6]) & good 

lowerHalf <- (comZ <= deciles[6]) & good

axisVector <- colMeans(poreCoordinates[upperHalf,]) - colMeans(poreCoordinates[lowerHalf,])

axisVector <- axisVector / sqrt(sum(axisVector^2)) # make it unit length


poreCoordinates <- cbind( comX[good],
                          comY[good],
                          comZ[good])


## project the centroid into the x-y plane along axisVector.
## allows us the parameterize the centroid using two parameters
## (x,y,0) rather than three (x,y,z)
centroid <- colMeans(poreCoordinates)

xyCentroid <- c( (axisVector[1]*-centroid[3])/axisVector[3] + centroid[1],
                 (axisVector[2]*-centroid[3])/axisVector[3] + centroid[2],
                 0)

## initial radius guess
r = 1000 # based on ideal coupon radius of 1000 micor-meters

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("rotateCoupon.R")
source("getRadius.R")


N <- length(poreCoordinates[,1])

radiusTarget <- rep(r, N) 



centerAxis <- NULL
attempt <- 1
while( is.null(centerAxis) && attempt <= 10 ) {
  print(l)
  attempt <- attempt + 1
  try(centerAxis <- nls(radiusTarget~radiusAligned(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY),
                  start = list(centroidX = xyCentroid[1], centroidY = xyCentroid[2],
                               axisVectorX = axisVector[1], axisVectorY = axisVector[2]),
                  control = nls.control(minFactor = 1/10000000000000)))
}
nlsCoeff <- coef(centerAxis)




## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
## useful for generating surface plots and histograms for each coupon
source("newCoupon.R")
setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/noisySpacing40to60/noise0.10")
nlsCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])

save(oldCoupon, nlsCoupon, nlsCoeff, centerAxis, file = paste0(round(l,3), "spacingNoise1.rda"))

}







