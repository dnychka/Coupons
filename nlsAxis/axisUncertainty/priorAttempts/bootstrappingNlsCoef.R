## 11/04/2019
##
##
##
## statistical inference for the nls coeff
## create CI for each parameter using non-parametric bootstrap
##
## uses simulated coupons with internal jitter 
## -------------------------------------------------------

library(conicfit)
library(rgl)
library(dplyr)

myDir <- "C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty"

# load the sample simulated coupon
trueCoupon <- readRDS(file.path(myDir, "/simulatedCouponInternalJitter.rds")) #has NOT been nls-ed

# knock it off its axis
aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}
theta = -40*pi/180

Ry <- aboutY(theta)

sampleCoupon <- trueCoupon %*% Ry 


## ---------------------------------------
## fit the model to the data using nls WITH 
## iterative median approach to optimize
## the parameter values
## ----------------------------------------

# path <- "C:/Users/barna/Documents/Coupons/nlsAxis/"
# setwd(path)
# filesToSource = list.files(pattern="*.R")
# sapply(filesToSource, source)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("initialParameters.R")
source("rotateCoupon.R")
source("getRadius.R")
source("newCoupon.R")
source("nlsAxisFit.R")


nlsObj <- nlsAxisFit(sampleCoupon)

nlsCoeff <- coef(nlsObj[[1]])


## scratch work, are the pores in the 
## same positions?

globalMedian <- median(nlsObj[[2]])

r <- sqrt(trueCoupon[,1]^2 + trueCoupon[,2]^2)
theta <- atan2(trueCoupon[,2], trueCoupon[,1])

median(r) #so original median radius value v close to nls found one...

nlsCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])


plot(nlsObj[[2]][1:100], pch=20, col="violetred1")
points(r[1:100])



# ahhhh the nlsCoupon is rotated by a factor of pi
# see this:
plot3d(trueCoupon[,1], trueCoupon[,2], trueCoupon[,3], type = "s", size =0.45,
       col = ifelse(trueCoupon[,1] == trueCoupon[1,1], "violetred1", "black"))
open3d()
plot3d(nlsCoupon[,1], nlsCoupon[,2], nlsCoupon[,3], type = "s", size = 0.45,
       col = ifelse(nlsCoupon[,1] == nlsCoupon[1,1], "violetred1", "black"))


aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}

nlsCouponZ <- nlsCoupon %*% aboutZ(pi)



## -------------------------
## create parameter values with residual based bootstrapping
## -------------------------

#err <- sampleCoupon - nlsCouponZ #errors for optimized param.
err <- trueCoupon - nlsCouponZ

## begin the bootstrap for loop ------------

nsamples = 1:1000

bootNlsCoef <- matrix(ncol = 5, nrow = length(nsamples), NA)

for(i in nsamples){
  print(i)
  
  j = 1
  while(j < 10){ #try nls up to 10 times

  # resample from nlsRes using non-parametric sample to get
  # new error values for every pore location
  bootErr <- err[sample(nrow(err),size=length(err[,1]),replace=TRUE),]
  
  # fabricate new data by adding the sampled errors to the 
  # model output
  
  bootCoupon <- nlsCouponZ + bootErr
  
  # estimate the parameter values for each new fabricated data set
  N <- length(bootCoupon[,1])
  
  startValues <- getInitialParameters(bootCoupon)
  
  
  distTarget <- rep(0,N)
  
  nlsObj <- try(nls(distTarget~getDistance(bootCoupon, 
                                           centroidX, centroidY, 
                                           axisVectorX, axisVectorY, r),
                    start = list(centroidX = startValues[4], 
                                 centroidY = startValues[5],
                                 axisVectorX = startValues[1], 
                                 axisVectorY = startValues[2],
                                 r = startValues[7]),
                    control =  nls.control(minFactor = 1/10000000000)))
  
  
 
  
  if(class(nlsObj) != "try-error"){
    j <- j + 1
    bootNlsCoef[i,] <- coef(nlsObj)
  }
  
  
  
  } #end of while loop
  
} #end of for loop

#saveRDS(bootNlsCoef, "bootNlsCoef.Rds")

origBootNlsCoef <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/bootNlsCoef.rds")


## create a bootstrap percentile CI

bootCentroidX <- quantile(bootNlsCoef[,1], prob = seq(0, 1, by= 0.1))
bCXlow <- bootCentroidX[1]
bCXup <- bootCentroidX[11]


bootCentroidY <- quantile(bootNlsCoef[,2], prob = seq(0, 1, by= 0.1))
bCYlow <- bootCentroidY[1]
bCYup <- bootCentroidY[11]


bootAxisVectorX <- quantile(bootNlsCoef[,3], prob = seq(0, 1, by= 0.1))
bAVXlow <- bootAxisVectorX[1]
bAVXup <- bootAxisVectorX[11]


bootAxisVectorY <- quantile(bootNlsCoef[,4], prob = seq(0, 1, by= 0.1))
bAVYlow <- bootAxisVectorY[1]
bAVYup <- bootAxisVectorY[11]


# 10th percentile coupon
lowBootCoupon <- newCoupon(sampleCoupon, bCXlow, bCYlow, bAVXlow, bAVYlow)

# 90th percentile coupon
upBootCoupon <- newCoupon(sampleCoupon, bCXup, bCYup, bAVXup, bAVYup)


plot(lowBootCoupon[,1], lowBootCoupon[,2], pch =20, col = "black",
     main = "XY projection of a=10% CI",
     xlab = "x", ylab = "y")
points(upBootCoupon[,1], upBootCoupon[,2], pch=1, col="cornflowerblue")
points(bCXlow, bCYlow, pch=20)
points(bCXup, bCYup, col="cornflowerblue")
