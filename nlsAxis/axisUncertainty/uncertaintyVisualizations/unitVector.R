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
library(FSA)
library(fields)

nlsCoeff1ring <- nlsCoeffS
bootNlsCoefS1ring <- bootNlsCoefS
bootNlsCoefS1ringCI <- bootCoefCIS


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/surfaces/nlsSurfaceData/F3")
load("F3surfaceData.rda")

nlsCoeffS <- nlsCoeff
bootNlsCoefS <- bootNlsCoef
bootCoefCIS <- bootCoefCI

rm(bootCoefCI, bootNlsCoef, nlsCoeff, bootRadius, nlsCoupon)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData")

load("coupon11poreCI.rda")

nlsCoeffP <- nlsCoeff
bootNlsCoefP <- bootNlsCoef
bootCoefCIP <- bootCoefCI

rm(nlsCoeff, bootNlsCoef, bootCoefCI, bootRadius, bootRadiusCI, nlsCoupon)



uDiffP <- vector()
uDiffS <- vector()


for(i in 1:1000){
  uDiffP[i] <- sqrt(sum(  (nlsCoeffP[3:4] - bootNlsCoefP[i,3:4])^2 ))
  uDiffS[i] <- sqrt(sum(  (nlsCoeffS[3:4] - bootNlsCoefS[i,3:4])^2 ))
}

hist(~uDiffP , col="grey20", w=0.0025,
     xlim = c(0,0.065), ylim = c(0,350),
     main = "uncertainty around center axis for both surface and pore fit, coupon F3",
     xlab = "distance from center axis")

legend("topright", c("pore fit", "surface fit"), pch = c(15,15), col = c("grey20", alpha("grey40", 0.3)))

hist(~uDiffS,add=T, w=0.0025,col=alpha("grey40",0.4) )

xline(quantile(uDiffS, c(0.025, 0.975)), lty=2, col=alpha("grey40",0.4))
xline(quantile(uDiffP, c(0.025, 0.975)), lty=2, col="grey20")


## plot out the possible center axes
## ---------------------------------

library(rgl)
library(fields)

load("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped/nlsCoupon4.rda")

rm(nlsCoupon, nlsCoeff)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData")

load("coupon4poreCI.rda")

uDiff <- vector()
for(i in 1:1000){
  uDiff[i] <- sqrt(sum(  (nlsCoeff[3:4] - bootNlsCoef[i,3:4])^2 ))
}

bigDiff <- order(uDiff, decreasing = T)
ctab <- two.colors(n=5, start="maroon3", end="cornflowerblue", middle="white",
                   alpha=1.0)

#ctab <- c("#CD2990", "#E369AF" ,"#F6A5CD", "#FFD5E6", "#FFF6F9", "#F7FFFF" ,
 #         "#DEF9FF" ,"#BBE1FF", "#91BDF8" ,"#6495ED")

axisVecs <- cbind(bootNlsCoef[,3:4], rep(1, 1000))
m <- bigDiff[1:5]

points3d(oldCoupon[,1], oldCoupon[,2], oldCoupon[,3], 
         size = 5, alpha = 0.1, zlim = c(-200, 4000))
for(i in 1:5){
abclines3d( bootNlsCoef[m[i],1], bootNlsCoef[m[i],2], 0, 
            axisVecs[m[i],1], axisVecs[m[i],2], axisVecs[m[i],3], 
            col = "cornflowerblue", alpha = 0.7, lwd = 2)
}
abclines3d(nlsCoeff[1], nlsCoeff[2], 0,
           nlsCoeff[3], nlsCoeff[4], 1, 
           col = "maroon3", lwd=4)
box3d(edges = "bbox", tick = T, box = TRUE, col = "black")




