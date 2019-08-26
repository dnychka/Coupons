##
##
##
##
##
##
##
## Bayseian uncertainty sampling
## -------------------------------------------------------

library(fields)
library(MASS)
library(FSA)
library(rgl)

setwd("C:/Users/barna/Documents/Coupons/datasets")
poreData <- readRDS("porosityData.rds")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("initialParameters.R")
source("rotateCoupon.R")
source("getRadius.R")
source("newCoupon.R")
source("cropCoupon.R")
source("nlsAxisFit.R")

layerSpacing <- seq(40,60,length.out = 150)


##--------------------------------------------------------------------
## get idea of RSE for synthetic coupons (optimizing to actual radius)
##--------------------------------------------------------------------
rse <- vector()
j=1

files <- list.files("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/test")

for(i in files){
load(paste0("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/test/",
            i))
  
  #load("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/test/40spacingSyn0.rda")
  
  rse[j] <- summary(centerAxis)$sigma
  j=j+1
}

plot(rse)

load("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/spacing40to60/testCenter.rda")
summary(centerAxis)$sigma

rseReal <- vector()
niceCoupon <- c(which(couponCov$polarAngle==0), which(couponCov$polarAngle==45))
for(n in niceCoupon){
  poreCoordinates <- cropCoupon(n, poreData)
  nlsObj <- nlsAxisFit(poreCoordinates)
  
  rseReal[n] <- summary(nlsObj)$sigma
}

plot(rseReal)

##--------------------------------------------------------------------

bxplots <- list()
w = 1
for (n in niceCoupon){

##--------------------------------------------------------------------
## fit the cylinder by LS and return a nls object
##--------------------------------------------------------------------

poreCoordinates <- cropCoupon(n, poreData)
firstFit <- nlsAxisFit(poreCoordinates)

## fit again by nls once interior pores are removed
interior <- median(firstFit[[2]])-sd(firstFit[[2]])
discard <- which(firstFit[[2]] < interior)

outerPores <- poreCoordinates[-discard,]
nlsObj <- nlsAxisFit(outerPores)

nlsObj <- nlsObj[[1]]

##--------------------------------------------------------------------

##--------------------------------------------------------------------
## OR pull the nls object from synthetic coupon
##--------------------------------------------------------------------
# nlsObj <- centerAxis
# poreCoordinates <- nlsCoupon
##--------------------------------------------------------------------



##--------------------------------------------------------------------
## compute beta_hat, covariance matrix, and estimated standard deviation
##--------------------------------------------------------------------
nlsSum <- summary(nlsObj)

betaHat <- coef(nlsObj)
Vbeta <- vcov(nlsObj)
sigmaHat <- nlsSum$sigma

n <- length(outerPores[,1]) # number of data points
k <- 4 # number of parameters

nSims <- 100 # number of simulations

sigma <- rep(NA, nSims)
beta <- matrix(NA, nrow = nSims, ncol = k)

##
## create nSims number of random simulations of beta and sigma
##

for (s in 1:nSims){
  sigma[s] <- sigmaHat*sqrt((n-k)/rchisq(1,n-k))
  beta[s,] <- mvrnorm (1, betaHat, Vbeta*sigma[s]^2)
  
  nlsCoupon <- newCoupon(poreCoordinates, beta[s,1], beta[s,2], 
                        beta[s,3], beta[s,4])
  
  saveRDS(nlsCoupon, 
          paste0("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/bayes/bayesCoupon",s,".rds"))
}
##--------------------------------------------------------------------



##--------------------------------------------------------------------
## compute signal strength for each orientation 
##--------------------------------------------------------------------
setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice")
source("findPeaks.R")
source("0degreeFunctions.R")

sampledSignals <- vector()

for(i in 1:nSims){
  coupon <- readRDS(paste0("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/bayes/bayesCoupon",i,".rds"))
  
  histSeries <- getHistogram(coupon)
  sampledSignals[i] <- getSignal(histSeries)
}


bxplots[[w]] <- sampledSignals
w=w+1
#boxplot(sampledSignals, ylim = c(0,1))
#pointEst <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosLayeredSignals.rds")
#yline(harmonicSignals[2], col = "tomato", lty=3, lwd=2)
##--------------------------------------------------------------------
}






