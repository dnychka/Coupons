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

load(paste0("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/spacing40to60/",
            round(layerSpacing[3],3),
            "spacingSyn0.rda"))

##--------------------------------------------------------------------
## fit the cylinder by LS and return a nls object
##--------------------------------------------------------------------
n = 4
poreCoordinates <- cropCoupon(n, poreData)
nlsObj <- nlsAxisFit(poreCoordinates)
##--------------------------------------------------------------------

##--------------------------------------------------------------------
## OR pull the nls object from synthetic coupon
##--------------------------------------------------------------------
nlsObj <- centerAxis
poreCoordinates <- newCoupon
##--------------------------------------------------------------------



##--------------------------------------------------------------------
## compute beta_hat, covariance matrix, and estimated standard deviation
##--------------------------------------------------------------------
nlsSum <- summary(nlsObj)

betaHat <- coef(nlsObj)
Vbeta <- vcov(nlsObj)
sigmaHat <- nlsSum$sigma

n <- length(poreCoordinates[,1]) # number of data points
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

boxplot(sampledSignals, ylim = c(0,1))
#pointEst <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosLayeredSignals.rds")
yline(harmonicSignals[2], col = "tomato", lty=3, lwd=2)
##--------------------------------------------------------------------







