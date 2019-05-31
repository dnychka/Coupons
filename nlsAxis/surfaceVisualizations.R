

library(rgl)
library(conicfit)


newCoupon <- readRDS("newCoupon.rds")
nlsCoeff <- readRDS("nlsCoeff.rds")


 
## calculate points on ideal coupon surface
circXY <- calculateCircle(nlsCoeff["axisVectorX"],
                nlsCoeff["axisVectorY"], 1000, steps = 100)

circZ <- seq(min(newCoupon[,3]), max(newCoupon[,3]), length.out = 200)

xyzcoord <- cbind(rep(circXY[,1], 200), rep(circXY[,2],200), rep(circZ, each = 100))

idealCyl <- cylinder3d(center = xyzcoord, radius = 2, closed = TRUE)

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



# plot3d(poreCoordinatesCrop[,1], poreCoordinatesCrop[,2], poreCoordinatesCrop[,3], type = "s", size = 0.45)
# points3d(poreCoordinates[,1], poreCoordinates[,2], poreCoordinates[,3], col = "magenta")
# 
# points3d(centroid[1], centroid[2], centroid[3], col = "violetred1", size = 4)
# points3d(xyCentroid[1], xyCentroid[2], xyCentroid[3], col = "cornflowerblue", size = 4)
# lines3d(rbind(centroid, xyCentroid))



# centroidX <- nlsCoeff["centroidX"]
# centroidY <- nlsCoeff["centroidY"]
# axisVectorX <- nlsCoeff["axisVectorX"]
# axisVectorY <- nlsCoeff["axisVectorY"]





# A trefoil knot
open3d()
theta <- seq(0, 2*pi, len = 25)
knot <- cylinder3d(
  center = cbind(
    sin(theta) + 2*sin(2*theta), 
    2*sin(3*theta), 
    cos(theta) - 2*cos(2*theta)),
  e1 = cbind(
    cos(theta) + 4*cos(2*theta), 
    6*cos(3*theta), 
    sin(theta) + 4*sin(2*theta)),
  radius = 0.8, 
  closed = TRUE)

shade3d(addNormals(subdivision3d(knot, depth = 2)), col = "green") 