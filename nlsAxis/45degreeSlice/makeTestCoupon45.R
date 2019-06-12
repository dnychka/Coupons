##
##
## construct a test case for the 45 degree slice project
## ------------------------------------------------------


library(conicfit)
library(rgl)


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")
load("nlsCoupon4.rda")
setwd("C:/Users/barna/Documents/Coupons/nlsAxis/datasets")
load("radiusIterations.rda") #to get median radius from nls fit

# calc points on a 45 degree line
x = seq(0,4000,by = 200)
y = rep(0, length(x))
z = x


phi = 45*pi/180
a = 814*sec(phi) #major axis
b = 814 #minor axis

pts <- matrix(NA, nrow = 25, ncol = 3)
pts[,1:2] <- calculateEllipse(x[1],y[1],a,b, steps = 25, randomDist = TRUE)
pts[,3] <- rep(z[1], 25)

ellipseStack <- pts

for(i in 2:length(x)){
  
pts <- matrix(NA, nrow = 25, ncol = 3)
pts[,1:2] <- calculateEllipse(x[i],y[i],a,b, steps = 25, randomDist = TRUE)
pts[,3] <- rep(z[i], 25)

ellipseStack <- rbind(ellipseStack, pts)

}

plot(ellipseStack[1:25,1:2])
points(ellipseStack[26:50, 1:2], pch = 20)

plot(ellipseStack[,1:2])
plot(ellipseStack[,1], ellipseStack[,3])

axes3d(edges = "bbox")
points3d(ellipseStack[,1],ellipseStack[,2],ellipseStack[,3], size = 0.45)


hist(ellipseStack[,3], breaks = 90, col = "grey30")


aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

phi <- -41*pi/180

Ry <- aboutY(phi)

turnedStack <- ellipseStack %*% Ry

plot(turnedStack[,1], turnedStack[,3])

axes3d(edges = "bbox", labels = TRUE)
plot3d(turnedStack[,1],turnedStack[,2],turnedStack[,3], size = 2)

#saveRDS(ellipseStack, "testEllipse.rda")



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

plot3d(oldCoupon[,1], oldCoupon[,2], oldCoupon[,3], type = "s", size = 0.45)
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




circZ <- seq(0, max(poreCoordinates[,3]), length.out = 200)

plot3d(poreCoordinates[,1], poreCoordinates[,2], poreCoordinates[,3], type = "s", size = 0.45)
points3d(centroid[1], centroid[2], centroid[3], size = 3, col = "tomato")
points3d(xyCentroid[1], xyCentroid[2], xyCentroid[3], size = 3, col = "magenta")
lines3d(axisVector[1],
        axisVector[2],
        circZ, col = "darkorange",
        lwd = 3)


j = 1

## initial radius guess
r = 1000 # based on ideal coupon radius of 1000 micor-meters

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("alignCoupon.R")
source("radiusAligned.R")


N <- length(poreCoordinates[,1])

  radiusTarget <- rep(r, N) 
  
  centerAxis <- nls(radiusTarget~radiusAligned(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY),
                    start = list(centroidX = xyCentroid[1], centroidY = xyCentroid[2],
                                 axisVectorX = axisVector[1], axisVectorY = axisVector[2]))
  
  
  nlsCoeff <- coef(centerAxis)
  
  ## compare radiusFinal with radiusInitial
  radiusInitial <- radiusAligned(poreCoordinates, centroid[1], centroid[2], axisVector[1], axisVector[2])
  
  ## run 'radiusAligned' with the optimal parameter values found by nls
  radiusFinal <- radiusAligned(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                               nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
  
  
  # store the initial median for comaprison (first iteration)
  #ifelse(diffMedians == 100, radiusFirstIter[[j]] <- radiusFinal, NA) 
  
  
  plot(radiusFinal, radiusInitial, main = paste(r))
  xline(median(radiusFinal), col = "cornflowerblue")
  
  radiusLastIter <- radiusFinal


## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
## useful for generating surface plots and histograms for each coupon
source("newCoupon.R")
setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")
newCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
n="test"
save(oldCoupon, newCoupon, nlsCoeff, radiusLastIter, file = paste0("nlsCoupon", n, ".rda"))

