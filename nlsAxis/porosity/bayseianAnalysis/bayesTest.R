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
library(astsa)


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


##--------------------------------------------------------------------
## look at RSE for the real coupons. Effect of interior crop on RSE
##--------------------------------------------------------------------

rseReal <- vector()
rseOuter <- vector()
niceCoupon <- c(which(couponCov$polarAngle==0), which(couponCov$polarAngle==45))

realCouponNlsObj <- list()
numPores <- vector()

firstFitNlsObj <- list()

for(n in niceCoupon){
  
  ## fit the cylinder by LS and return a nls object
  ##------------------------------------------------
  poreCoordinates <- cropCoupon(n, poreData)
  firstFit <- nlsAxisFit(poreCoordinates)
  
  rseReal[n] <- summary(firstFit[[1]])$sigma
  
  firstFitNlsObj[[n]] <- firstFit[[1]]
  
  ## fit again by nls once interior pores are removed
  ##------------------------------------------------
  interior <- median(firstFit[[2]])-sd(firstFit[[2]])
  discard <- which(firstFit[[2]] < interior)
  
  outerPores <- poreCoordinates[-discard,]
  nlsObj <- nlsAxisFit(outerPores)
  numPores[n] <- length(outerPores[,1])
  
  rseOuter[n] <- summary(nlsObj[[1]])$sigma
  
  ## save the no-interior-pores nls object
  ##------------------------------------------------
  realCouponNlsObj[[n]] <- nlsObj[[1]] 

}

##--------------------------------------------------------------------
## a few quick histograms to show no interior pores
##--------------------------------------------------------------------

hist(~firstFit[[2]], w=10, col = "grey30",
     freq = FALSE, xlab = "radius (microns)",
     main = "pore radii distribution with interior pores",
     xlim = c(100,900))

hist(~nlsObj[[2]], w = 10, col = "grey30",
     freq = FALSE, xlab = "radius (microns)",
     main = "pore radii distribution without interior pores",
     xlim = c(100,900))


##--------------------------------------------------------------------
## check to make sure interior pores aren't influencing coeff fit
##--------------------------------------------------------------------

recordTTest <- matrix(nrow=4,ncol=length(niceCoupon))
firstTemp <- matrix(NA,4,4)
realTemp <- matrix(NA,4,4)

k=1

for(i in niceCoupon){

  firstTemp <-  summary(firstFitNlsObj[[i]])$coefficient
  
  realTemp <- summary(realCouponNlsObj[[i]])[[10]]
  
  for(j in 1:4){
    tObj <- t.test(firstTemp[j,], realTemp[j,])
    
    recordTTest[j,k] <- tObj[[3]] 
  }
  
  k=k+1
}


##--------------------------------------------------------------------
## RSE plots
##--------------------------------------------------------------------
# with color for polar angle
plot(rseReal, ylim = c(0,150),
     main = "RSE for nls fit, with and without interior pores",
     ylab = "RSE", xlab = "coupon number",
     col = ifelse(couponCov$polarAngle == 0, "cornflowerblue", "violetred1"))
points(rseOuter, 
       col = ifelse(couponCov$polarAngle == 0, "cornflowerblue", "violetred1"), pch = 16)
legend("topright", c("0 degree", "45 degree"), col = c("cornflowerblue", "violetred1"), pch = c(16,16))

# without color for polar angles
plot(rseReal, ylim = c(0,150),
     main = "RSE for nls fit, with and without interior pores",
     ylab = "RSE", xlab = "coupon number")
points(rseOuter, 
        pch = 16)
legend("topright", c("interior pores", "no interior pores"),  pch = c(1,16))



##--------------------------------------------------------------------

##--------------------------------------------------------------------
## UNCERTAINTY TESTING
##--------------------------------------------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerFunctions")
source("calculateTestStatistic.R")
source("rotations.R")


uncertaintyInSignals <- matrix(nrow = length(niceCoupon), ncol = 100)
zeros <- which(couponCov$polarAngle==0)
fortyFives <- which(couponCov$polarAngle==45)
j=1

for (n in fortyFives){
  
  print(n)

nlsObj <- realCouponNlsObj[[n]]


##--------------------------------------------------------------------
## compute beta_hat, covariance matrix, and estimated standard deviation
##--------------------------------------------------------------------
nlsSum <- summary(nlsObj)

betaHat <- coef(nlsObj)
Vbeta <- vcov(nlsObj)
sigmaHat <- nlsSum$sigma

m <- numPores[n] # number of data points
k <- 4 # number of parameters

nSims <- 100 # number of simulations

sigma <- rep(NA, nSims)
beta <- matrix(NA, nrow = nSims, ncol = k)

##
## create nSims number of random simulations of beta and sigma
##

bayesCoupon <- list()

for (s in 1:nSims){
  sigma[s] <- sigmaHat*sqrt((m-k)/rchisq(1,m-k))
  beta[s,] <- mvrnorm (1, betaHat, Vbeta*sigma[s]^2)
  
  bayesCoupon[[s]] <- newCoupon(poreCoordinates, beta[s,1], beta[s,2], 
                        beta[s,3], beta[s,4])

}
##--------------------------------------------------------------------



##--------------------------------------------------------------------
## compute signal strength for each orientation 
##--------------------------------------------------------------------

sampledSignals <- rep(NA, nSims)
angleSeq <- seq(0, 2*pi, length.out = 100) #increment by which we change rotation
tempTS <- rep(NA, 100) #store TS for each rotation

for(s in 1:nSims){
  coupon <- bayesCoupon[[s]]
  
  if(n %in% zeros){
  sampledSignals[s] <- getTestStatistic(coupon)
  }
  else {
    m = 1
    for(w in angleSeq){
      
      rotatedCoupon <- rotateFortyFive(w, coupon) #center and rotate the coupon
      
      tempTS[m] <- getTestStatistic(rotatedCoupon) #calc TS for rotated coupon
      
      print(m)
      m=m+1
    }
    sampledSignals[s] <- max(tempTS)
  }
}


uncertaintyInSignals[j,] <- sampledSignals
j = j+1

##--------------------------------------------------------------------
}


# saveRDS(uncertaintyInSignals, "porosityUncertainty.rds")

TSZlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZlayers.rds")
TSZreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZreal.rds")

TSZuncertain <- uncertaintyInSignals[1:16,]

TSZ <- rbind(TSZuncertain, TSZreal)

boxplot(t(TSZ))

boxplot(as.vector(TSZuncertain))


TSFun <- as.vector(uncertaintyInSignals[1:24,])

TSFun <- TSFun[order(TSFun, decreasing=T)]

boxplot(TSFun[1:500])
boxplot(TSFreal, add=T, border = alpha("black", 0.6))


