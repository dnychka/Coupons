## 6/27/2019
##
##
##
## working document for getting relative peak intenstiy
## from real 45 degree coupons
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("findPeaks.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/datasets")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData")
nullSignals <- readRDS("45sNullSignals.rds")
layeredCoupons <- readRDS("45sLayeredSignals.rds")


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData/cropped")

## ------------------------------------------------------

fortyFive <- which(couponCov$polarAngle == 45)

angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta

numCoupons <- length(fortyFive)
realSignals <- matrix(NA, nrow = length(angleSeq), ncol = numCoupons)

i=1
for(n in fortyFive){
    
    print(n)
    
    load(paste0("nlsCoupon", n, ".rda"))
    
    # center the coupon
    centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                          scale(newCoupon[,2], center = TRUE, scale = FALSE),
                          newCoupon[,3])
    
    
    histSave <- FFrotate(angleSeq, centerCoupon)
    
    realSignals[,i] <-  makePeriodogram(angleSeq, histSave)
    
  i=i+1
} 


layeredCoupons <- apply(layeredCoupons,2,max)
nullCoupons <- c(apply(nullSignals[-1,], 2, max), rep(NA,
                                                      (length(layeredCoupons)-length(apply(nullSignals[-1,],2,max)))))
realCoupons <- c(apply(realSignals,2,max), 
                 rep(NA, (length(nullCoupons)-length(apply(realSignals,2,max)))))

dat <- cbind(nullCoupons, realCoupons)

dev.off()
boxplot(dat, main = "relative peak strength, 45 degree coupons",
        pch = 16, boxwex = 0.7)


dat <- cbind(dat, layeredCoupons)

boxplot(dat, main = "relative peak strength, 45 degree coupons",
        pch = 16, boxwex = 0.7)

t.test(nullCoupons, realCoupons)
t.test(layeredCoupons, realCoupons)
