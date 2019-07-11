



library(rgl)
library(conicfit)
library(fields)
library(useful)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")


load("nlsCoupon41.rda")


## calculate points on ideal coupon surface
circXY <- calculateCircle(nlsCoeff["axisVectorX"],
                          nlsCoeff["axisVectorY"], 1000, steps = 100)

circZ <- seq(min(newCoupon[,3]), max(newCoupon[,3]), length.out = 80)

xyzcoord <- cbind(rep(circXY[,1], 80), rep(circXY[,2],80), rep(circZ, each = 100))

xyzcoord <- scale(xyzcoord, scale = TRUE, center = TRUE) + 2

xyzcoord[,3] <- xyzcoord[,3] + 3.5

newCoupon <- scale(newCoupon, scale = TRUE, center = TRUE) + 2

newCoupon[,3] <- newCoupon[,3] + 4

aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutX <- function(phi) {matrix( c(1, 0, 0, 0, cos(phi), sin(phi), 0, -sin(phi), cos(phi)), 3, 3)}

phi = 20*pi/180

Ry <- aboutY(phi)
Rx <- aboutX(phi)

TnewCoupon <- newCoupon %*% Ry %*% Rx

axisMat <- cbind(rep(2,length(xyzcoord[,1])), rep(2,length(xyzcoord[,1])), 
                 seq(0,max(xyzcoord[,3]), length.out = length(xyzcoord[,1])))

psi <- 10.7*pi/180

Ry <- aboutY(psi)
Rx <- aboutX(psi)

Txyzcoord <- xyzcoord %*% Ry %*% Rx 

Txyzcoord[,1:2] <- Txyzcoord[,1:2] + 0.5

Taxis <- axisMat %*% Ry %*% Rx

Taxis <- Taxis[-Taxis[,3]<0,]

idealCyl <- cylinder3d(center = Txyzcoord, radius = 0.02, closed = TRUE)

open3d()
par3d(cex=0.7)
points3d(TnewCoupon[,1],
       TnewCoupon[,2],
       TnewCoupon[,3],
       size = 2)
points3d(Taxis[,1]+0.5,
        Taxis[,2]+0.5,
        Taxis[,3], col = "darkorange",
        size = 2)
points3d(min(Taxis[,1]+0.5),
         min(Taxis[,2]+0.5),
         min(Taxis[,3]), col = "tomato",
         size = 7)
rgl.light(90, 90)
rgl.material(alpha = 0.2, lit = FALSE)
shade3d(idealCyl, col = "cornflowerblue")
rgl.lines(c(0, 8), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(0,10), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(0, 0), c(0,8), color = "black")
axes <- rbind(c(8, 0, 0), c(0, 10, 0), 
              c(0, 0, 8))
rgl.points(axes, color = "black", size = 3)
rgl.texts(axes, text = c("x", "y","z"), color = "black",
          adj = c(0.5, -0.8), size = 2)




## --------------------------------------------------------
## radial dispersion example
## --------------------------------------------------------

load("nlsCoupon4.rda")

par(pty = "s")
par(mfrow = c(1,2))
plot(oldCoupon[,1], oldCoupon[,2], pch = 20, 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(ylab="y axis", line=1, cex.lab=1.2)
title(xlab="x axis", line=1, cex.lab=1.2)
lmFit <- LMcircleFit(cbind(oldCoupon[,1], oldCoupon[,2]))
draw.circle(lmFit[1], lmFit[2], 1000, lwd = 2, border = "cornflowerblue")

plot(newCoupon[,1], newCoupon[,2], pch = 20, 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
title(ylab="y axis", line=1, cex.lab=1.2)
title(xlab="x axis", line=1, cex.lab=1.2)
lmFit <- LMcircleFit(cbind(newCoupon[,1], newCoupon[,2]))
draw.circle(lmFit[1], lmFit[2], 814, lwd = 2, border = "cornflowerblue")



plot(newCoupon[,1], newCoupon[,2])

plot3d(oldCoupon[,1], oldCoupon[,2], oldCoupon[,3], type = "s", size = 0.45)



