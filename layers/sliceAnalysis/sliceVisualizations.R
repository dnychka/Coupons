##
##
##
##
## make plots for section "layers" in thesis write up
##
## ---------------------------------------------------------------------



##
## very obvious layers in pore structure to highlight 
## differences in the 45 and 0 degree coupons
## --------------------------------------------------

zeroEx <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/layerExample0.rds")
fEx <- readRDS("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/layerExample45.rds")


plot3d(fEx[,1], fEx[,2], fEx[,3], type = "s", size = 0.45)

## zero coupon
circXY <- calculateCircle(0,
                          0, 1000, steps = 100)

circZ <- seq(1650, 1850, length.out = 200)

xyzcoord <- cbind(rep(circXY[,1], 200), rep(circXY[,2],200), rep(circZ, each = 100))

idealCylz <- cylinder3d(center = xyzcoord, radius = 2, closed = TRUE)

open3d()
par3d(cex=0.7)
#axes3d(edges="bbox")
points3d(zeroEx[1:750,1],
       zeroEx[1:750,2],
       zeroEx[1:750,3],
        size = 2, axes=FALSE,
       xlab = "", ylab= "", zlab= "")
lines3d(0,
        0,
        seq((min(zeroEx[,3])-400), (max(zeroEx[,3])+200), length.out = 300), col = "darkorange",
        lwd = 3)
rgl.material(alpha = 0.1, lit = FALSE)
shade3d(idealCylz, col = "cornflowerblue")
rgl.viewpoint(theta=0, phi=90)





##forty-five coupon (upright)
phi = 45*pi/180
a = 1000*sec(phi) #major axis
b = 1000 #minor axis


ellpXY <- calculateEllipse(0,0,a,b, steps = 25, randomDist = FALSE)

ellpZ <- seq(-2807, -3025, length.out = 200)

xyzEllp <- cbind(rep(ellpXY[,1], 200), rep(ellpXY[,2],200), rep(ellpZ, each = 100))

Rxyzcoord <- xyzEllp %*% aboutY(-42*(pi/180))

Rxyzcoord[,1] <- Rxyzcoord[,1] - mean(Rxyzcoord[,1])

idealCyle <- cylinder3d(center = Rxyzcoord, radius = 2, closed = TRUE)

open3d()
par3d(cex=0.7)
#axes3d(edges="bbox")
points3d(fEx[,1],
       fEx[,2],
       fEx[,3],
       size = 2, axes=FALSE,
       xlab = "", ylab= "", zlab= "")
lines3d(0,
        0,
        seq((min(fEx[,3])-600), (max(fEx[,3])+400), length.out = 300), col = "darkorange",
        lwd = 3)
rgl.material(alpha = 0.1, lit = FALSE)
shade3d(idealCyle, col = "cornflowerblue")
rgl.viewpoint(theta=0, phi=90)


##forty-five coupon (at a tilt)
phi = 45*pi/180
a = 1000*sec(phi) #major axis
b = 1000 #minor axis


theta = 0
Rz <- aboutZ(theta)

RzCoupon <- fEx %*% Rz


phi <- 45*pi/180 #45 degree tilt

Ry <- aboutY(phi)

RzyCoupon <- RzCoupon %*% Ry

axisMat <- cbind(rep(0,length(fEx[,1])), rep(2,length(fEx[,1])), 
                 seq((min(fEx)-600),(max(fEx[,3])+400), length.out = length(fEx[,1])))

Taxis <- axisMat %*% Rz %*% Ry

ellpXY <- calculateEllipse(Taxis[350,1],Taxis[350,2],a,b, steps = 25, randomDist = FALSE)

ellpZ <- seq(Taxis[328,3], Taxis[360,3], length.out = 200)

xyzEllp <- cbind(rep(ellpXY[,1], 200), rep(ellpXY[,2],200), rep(ellpZ, each = 100))


idealCyle <- cylinder3d(center = xyzEllp, radius = 2, closed = TRUE)


plot3d(RzyCoupon[,1], RzyCoupon[,2], RzyCoupon[,3], type = "s", size = 0.45)

open3d()
par3d(cex=0.7)
#axes3d(edges="bbox")
points3d(RzyCoupon[,1],
         RzyCoupon[,2],
         RzyCoupon[,3],
         size = 2, axes=FALSE,
         xlab = "", ylab= "", zlab= "")
points3d(Taxis[,1],
         Taxis[,2],
         Taxis[,3], col = "darkorange",
         size = 2)
rgl.material(alpha = 0.1, lit = FALSE)
shade3d(idealCyle, col = "cornflowerblue")
rgl.viewpoint(theta=0, phi=90)



