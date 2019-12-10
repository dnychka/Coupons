##
##
##
##
##
##
## create cool graphics with pore uncertainty
## -------------------------------------------------

library(scales)
library(conicfit)

draw <- sample(1:1000, 50, replace=F)

#sampDraw <- draw

#draw <- 1:1029

par(pty="s")


Xlo <- vector()
Ylo <- vector()

Xup <- vector()
Yup <- vector()

j = 1

for (i in draw){
 
  Xlo[j] <- pivotCIpore[1,i] * cos(nlsTheta[i])
  Xup[j] <- pivotCIpore[2,i] * cos(nlsTheta[i])
  
  Ylo[j] <- pivotCIpore[1,i] * sin(nlsTheta[i])
  Yup[j] <- pivotCIpore[2,i] * sin(nlsTheta[i])
  
  j = j + 1
   
}



circXY <- calculateCircle(nlsCoeff["axisVectorX"],
                          nlsCoeff["axisVectorY"], 1000, steps = 100)



plot(nlsCoupon[draw,1], nlsCoupon[draw,2], 
     pch=20,
     ylim = c(min(circXY[,2]), max(circXY[,2])),
     xlim = c(min(circXY[,1]), max(circXY[,1])),
     main = "subset of pores and confidence intervals, coupon E22",
     xlab = "x",
     ylab = "y")
lines(circXY, col = "cornflowerblue", lwd=2)
points(nlsCoeff[3], nlsCoeff[4], pch=16, col = "darkorange", cex=1.4)
arrows(x0=Xlo, y0=Ylo, 
       x1=Xup, y1=Yup, 
       code=3, angle=90, length=0.03, lty =1, lwd = 2, col = "red")
points(nlsCoupon[draw,1], nlsCoupon[draw,2], pch=16)



plot(order(nlsCoupon[draw,1]), nlsCoupon[,3])



## 
## make the r vs theta plots
## -----------------------------------------------
draw <- sample(1:1927,20, replace=F)

plot(nlsTheta[draw], sqrt(nlsCoupon[,1]^2 + nlsCoupon[,2]^2)[draw], 
     pch=16, col="darkorange", cex=1,
     main = "coupon 4, 95% confidence intervals",
     xlab = "theta",
     ylab = "radius")

arrows(x0=nlsTheta[draw], y0=pivotCIpore[1,draw], 
       x1=nlsTheta[draw], y1=pivotCIpore[2,draw], 
       code=3, angle=90, length=0.1, lty = 2, col = "steelblue")



widthCI <- apply(pivotCIpore, 2, diff)







