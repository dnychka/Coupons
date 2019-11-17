## 10/04/2019
##
##
##
## add internal random noise to synthetic 0 degree coupons
## -------------------------------------------------------

library(conicfit)
library(rgl)
library(dplyr)

l = 25 #number of layers 
z = seq(0,4000,by = l)
y = rep(0, length(z))
x = y

a = 814 #radius

rnoise <- 50

pts <- matrix(NA, nrow = 25, ncol = 3)
pts[,1:2] <- calculateEllipse(x[1],y[1],a,a, steps = 25, noiseFun=function(x) 
  (x+runif(1,-rnoise,rnoise)))
pts[,3] <- rep(z[1], 25)


ellipseStack <- pts

for(i in 2:length(x)){
  
  pts <- matrix(NA, nrow = 25, ncol = 3)
  pts[,1:2] <- calculateEllipse(x[i],y[i],a,a, steps = 25, noiseFun=function(x) 
    (x+runif(1,-rnoise,rnoise)))
  pts[,3] <- rep(z[i], 25)
  
  ellipseStack <- rbind(ellipseStack, pts)
  
}


plot3d(ellipseStack[,1], ellipseStack[,2], ellipseStack[,3],
       type = "s", size=0.45)
plot(ellipseStack[,1], ellipseStack[,2])

## store the OG radius values
r <- sqrt(ellipseStack[,1]^2 + ellipseStack[,2]^2)
theta <- atan2(ellipseStack[,2], ellipseStack[,1])
result <- data_frame(r = r, theta = theta, x = ellipseStack[,1], y = ellipseStack[,2])
result$theta <- result$theta + (result$y < 0) * 2 * pi

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty")
saveRDS(ellipseStack, "simulatedCouponInternalJitter.rds")


##
## check that our alignment method works 
## uses only one iteration of nls
## -------------------------------------------------------------------


deciles <- quantile(ellipseStack[,3], prob = seq(0, 1, length = 11), type = 5)

comX <- ellipseStack[,1]
comY <- ellipseStack[,2]
comZ <- ellipseStack[,3]


## get centers of mass for upper and lower half of the coupon to 
## compute axisVector, which is the direction vector for 
## the initial axis estimate
upperHalf <- (comZ >= deciles[6])

lowerHalf <- (comZ <= deciles[6])

axisVector <- colMeans(ellipseStack[upperHalf,]) - colMeans(ellipseStack[lowerHalf,])

axisVector <- axisVector / sqrt(sum(axisVector^2)) # make it unit length



## project the centroid into the x-y plane along axisVector.
## allows us the parameterize the centroid using two parameters
## (x,y,0) rather than three (x,y,z)
centroid <- colMeans(ellipseStack)

xyCentroid <- c( (axisVector[1]*-centroid[3])/axisVector[3] + centroid[1],
                 (axisVector[2]*-centroid[3])/axisVector[3] + centroid[2],
                 0)

## initial radius guess
r = 1000 # based on ideal coupon radius of 1000 micor-meters

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("rotateCoupon.R")
source("getRadius.R")


N <- length(ellipseStack[,1])

radiusTarget <- rep(r, N) 



centerAxis <- nls(radiusTarget~radiusAligned(ellipseStack, centroidX, centroidY, axisVectorX, axisVectorY),
                        start = list(centroidX = xyCentroid[1], centroidY = xyCentroid[2],
                                     axisVectorX = axisVector[1], axisVectorY = axisVector[2]),
                        control = nls.control(minFactor = 1/10000000000000))
  
nlsCoeff <- coef(centerAxis)

source("newCoupon.R")
nlsCoupon <- newCoupon(ellipseStack, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])


oldR <- abs(result$r - 1000)

newR <- residuals(centerAxis)

#new polar coords for nls aligned coupon
r <- sqrt(nlsCoupon[,1]^2 + nlsCoupon[,2]^2)
theta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
newResult <- data_frame(r = r, theta = theta, x = nlsCoupon[,1], y = nlsCoupon[,2])
newResult$theta <- newResult$theta + (newResult$y < 0)  * pi #NOTE: this is usualyl 2pi but
                                                             #nls rotated the coupon pi around z axis
                                                             #this time so to make pts match up

plot(result$theta[1:10], oldR[1:10], ylim = c(120, 250),
     ylab = "residual value", xlab = "theta",
     main = "'known' residuals vs nls residuals")
points(newResult$theta[1:10], newR[1:10], pch = 16, col = ifelse(newR > oldR, "violetred1", "tomato"))

