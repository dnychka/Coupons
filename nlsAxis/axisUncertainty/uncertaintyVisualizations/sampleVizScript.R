##
##
##
##
##
##
##
## script to run interactive 3D visual of pore network
## in two coupons: 'realCoupon', the data collected from
## XCT scans, and 'simulatedCoupon', data simulated using
## a bivariate kernel density estimate (KDE) on z and theta
## --------------------------------------------------------


install.packages(rgl)
library(rgl)

## put path here
myDir <- #"C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData"

## load data file 'couponsUsingKDESampling.rda'
load(file.path(myDir, "/couponsUsingKDESampling.rda"))


open3d()
plot3d(simulatedCoupon[,1], simulatedCoupon[,2], simulatedCoupon[,3], type = "s", size = 0.45,
       xlab = "simulated coupon", ylab= "", zlab="")


open3d()
plot3d(realCoupon[,1], realCoupon[,2], realCoupon[,3], type = "s", size=0.45,
       xlab = "real coupon", ylab=" ", zlab=" ")









