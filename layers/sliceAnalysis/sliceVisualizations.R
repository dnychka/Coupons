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

circZ <- seq(2650, 2850, length.out = 200)

xyzcoord <- cbind(rep(circXY[,1], 200), rep(circXY[,2],200), rep(circZ, each = 100))

idealCyl <- cylinder3d(center = xyzcoord, radius = 2, closed = TRUE)

open3d()
par3d(cex=0.7)
axes3d(edges="bbox")
plot3d(zeroEx[,1],
       zeroEx[,2],
       zeroEx[,3],
       type="s", size = 0.45, axes=FALSE,
       xlab = "", ylab= "", zlab= "")
lines3d(0,
        0,
        seq((min(zeroEx[,3])-600), (max(zeroEx[,3])+400), length.out = 300), col = "darkorange",
        lwd = 3)
rgl.material(alpha = 0.1, lit = FALSE)
shade3d(idealCyl, col = "cornflowerblue")





##forty-five coupon
circXY <- calculateCircle(0,
                          0, 814, steps = 100)

circZ <- seq(-2650, -2850, length.out = 200)

xyzcoord <- cbind(rep(circXY[,1], 200), rep(circXY[,2],200), rep(circZ, each = 100))

Rxyzcoord <- xyzcoord %*% aboutY(-45)

Rxyzcoord[,1] <- Rxyzcoord[,1] - mean(Rxyzcoord[,1])

idealCyl <- cylinder3d(center = Rxyzcoord, radius = 2, closed = TRUE)

open3d()
par3d(cex=0.7)
axes3d(edges="bbox")
points3d(fEx[,1],
       fEx[,2],
       fEx[,3],
       size = 0.2, axes=FALSE,
       xlab = "", ylab= "", zlab= "")
lines3d(0,
        0,
        seq((min(fEx[,3])-600), (max(fEx[,3])+400), length.out = 300), col = "darkorange",
        lwd = 3)
rgl.material(alpha = 0.1, lit = FALSE)
shade3d(idealCyl, col = "cornflowerblue")




