##
##
##
##
##
##
## exploratory anlaysis: uncertainty around unit vector
## ------------------------------------------------------------------

library(rgl)
library(scales)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/uncertaintyData")

load("coupon4poreCI.rda")

uDiff <- vector()

for(i in 1:1000){
  uDiff[i] <- sqrt(sum(  (nlsCoeff[3:4] - bootNlsCoef[i,3:4])^2 ))
}

hist(uDiff)


quantile(uDiff,c(.025,.975))


## plot out the possible center axes
## ---------------------------------
axisVecs <- cbind(bootNlsCoef[,3:4], rep(1, 1000))
m <- 10

points3d(oldCoupon[,1], oldCoupon[,2], oldCoupon[,3], 
         size = 5, alpha = 0.1)
abclines3d( bootNlsCoef[1:m,1], bootNlsCoef[1:m,2], rep(0,m), 
            axisVecs[1:m,1], axisVecs[1:m,2], axisVecs[1:m,3], alpha=0.7)
abclines3d(nlsCoeff[1], nlsCoeff[2], 0,
           nlsCoeff[3], nlsCoeff[4], 1, 
           col = "maroon3", lwd=4)
box3d(edges = "bbox", tick = T, box = TRUE, col = "black")




