## 
## 3 June 2019
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## create 3d plots of select coupons, visulize both nls center axis and 
## ideal cylinder surface
## Also create unwrapped surface plots of coupon
##
##------------------------------------------------------------------------------------------------------------------------

library(rgl)
library(conicfit)
library(fields)
library(useful)
library(plyr)
library(dplyr)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData/cropped")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/45degreeSlice/45degreeData")


load("45spacingSyn45.rda")

load("nlsCoupon56.rda")

load("C:/Users/barna/Documents/Coupons/nlsAxisCopy/surfaces/nlsSurfaceTest.rda")

 
## calculate points on ideal coupon surface
circXY <- calculateCircle(nlsCoeff["axisVectorX"],
                nlsCoeff["axisVectorY"], 1000, steps = 100)

circZ <- seq(min(newCoupon[,3]), max(newCoupon[,3]), length.out = 200)

xyzcoord <- cbind(rep(circXY[,1], 200), rep(circXY[,2],200), rep(circZ, each = 100))

idealCyl <- cylinder3d(center = xyzcoord, radius = 2, closed = TRUE)

open3d()
par3d(cex=0.7)
plot3d(oldCoupon[,1],
       oldCoupon[,2],
       oldCoupon[,3],
       type = "s", size = 0.45,
       xlab = " ", ylab = " ", zlab = " ")

open3d()
par3d(cex=0.7)
plot3d(newCoupon[,1],
       newCoupon[,2],
       newCoupon[,3],
       type = "s", size = 0.45,
       zlab = "z", xlab = "x", ylab = "y")
lines3d(nlsCoeff["axisVectorX"],
        nlsCoeff["axisVectorY"],
        circZ, col = "darkorange",
        lwd = 3)
rgl.material(alpha = 0.5, lit = FALSE)
shade3d(idealCyl, col = "cornflowerblue")



## --------------------------------------------------
## unzipped coupon surface

dev.off()

polar <- useful::cart2pol(newCoupon[,1], newCoupon[,2])

polar <- polar[-order(polar$r, decreasing = TRUE)[1:5],]

zAx <- newCoupon[,3][-order(polar$r, decreasing = TRUE)[1:5]]

plot(polar$theta, newCoupon[,3], 
     col = color.scale(polar$r),  pch = 16,
     main = "surface plot, test coupon",
     xlab = "angle (radians)", ylab = "z axis",
     cex=2)


 plot(polar$theta, polar$r, pch = 20)
 
 
 
 
 ## center the data
 X <- newCoupon[,1] - nlsCoeff["axisVectorX"]
 Y <- newCoupon[,2] - nlsCoeff["axisVectorY"]
 
 r <- sqrt(X^2 + Y^2)
 theta <- atan2(Y, X)
 result <- data_frame(r = r, theta = theta, x = X, y = Y)
 result$theta <- result$theta + (result$y < 0) * 2 * pi
 
 plot(result$theta, newCoupon[,3], 
      col = color.scale(result$r),  pch = 16,
      main = "surface plot, test coupon",
      xlab = "angle (radians)", ylab = "z axis",
      cex=2)
 
 

#max coupon radius is 1328.943, but setting this as zlim makes
#surface plots mostly yellow and hard to read

