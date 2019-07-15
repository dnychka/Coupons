
##
##
## working document for experimenting with layers in
## 90 degree coupons
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)
library(plyr)
library(dplyr)

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice")
source("findPeaks.R")

setwd("C:/Users/barna/Documents/Coupons/datasets")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

## ------------------------------------------------------

ninety <- which(couponCov$polarAngle == 90)

numCoupons <- length(ninety)
realSignals <- rep(NA, ncol = length(ninety))

m=1
for(n in ninety){
  print(n)
  
  load(paste0("nlsCoupon", n, ".rda"))
  
  
  ## center the data
  X <- newCoupon[,1] - nlsCoeff["axisVectorX"]
  Y <- newCoupon[,2] - nlsCoeff["axisVectorY"]
  
  r <- sqrt(X^2 + Y^2)
  theta <- atan2(Y, X)
  result <- data_frame(r = r, theta = theta, x = X, y = Y)
  result$theta <- result$theta + (result$y < 0) * 2 * pi
  
  plot(result$theta, newCoupon[,3], 
       col = color.scale(result$r),  pch = 16,
       main = paste0(n),
       xlab = "angle (radians)", ylab = "z axis")
  
  boundBox <- which(result$theta >=3.7 & result$theta <=4)
  
  h<-hist(~result$theta[boundBox], w=.01, plot = FALSE)
  
  plot(h, col = "grey40")
  
  histSave <- rbind(h$counts, h$mids)
  
  
  threshold <- 10
  influence <- 0
  
  xGrid <- seq(min(histSave[2,]), max(histSave[2,]), length.out = length(histSave[1,]))
  
  lowpass <- splint(histSave[2,], histSave[1,], xGrid, lambda = 0.001)
  
  lines(xGrid, lowpass, col="tomato", lwd=2)

  frDt <- fft(histSave[1,]-lowpass)
  
  Dtlen <- length(histSave[1,])
  
  Fr <- (1:Dtlen/Dtlen)[1:(Dtlen/2)]
  P <- (Mod(2*frDt/Dtlen)^2)[1:(Dtlen/2)]
  
  plot(Fr,P,type="l", main = paste0("layer spacing at ", round(n,3)))
  
  lag <- length(Fr)/10
  
  peak <- findFreq(P,lag,threshold,influence)
  
  ind <- which(peak$signals==1)
  
  
  
  
  findHarmonic <- vector()
  findBand <- vector()
  
  if(length(ind) != 0){# check to make sure there's actually peaks
    
    Ppeaks <- P[ind]
    Frpeaks <- Fr[ind]
    
    origFreq <- Frpeaks[which.max(Ppeaks)]
    
    for(j in c(1,1/2,1/3,1/4,1/5)){
      fundFreq <- origFreq*j
      if(fundFreq <= .2){break}
    }
    ##find the harmonics
    for(j in 1:5){
      for(i in 1:length(Fr)){
        
        if(isTRUE(all.equal(fundFreq*j, Fr[i], tol = 0.05))){
          if(Fr[i] %in% findBand){
            NA #don't want to repeat values
          } else {findHarmonic <- c(findHarmonic, P[i])
          findBand <- c(findBand, Fr[i])}
        } else {NA}
        
      }
    }
    xline(findBand, col = "violetred1", lty=3)
    
  } else {
    fundFreq <- Fr[which.max(P[which(Fr<=0.2)])]
    ##find the harmonics
    for(j in 1:5){
      for(i in 1:length(Fr)){
        
        if(isTRUE(all.equal(fundFreq*j, Fr[i], tol = 0.05))){
          if(Fr[i] %in% findBand){
            NA #don't want to repeat values
          } else {findHarmonic <- c(findHarmonic, P[i])
                  findBand <- c(findBand, Fr[i])}
        } else {NA}
             
       }
    }
    xline(findBand, col = "violetred1", lty=3)
  }

 realSignals[m] <- sum(findHarmonic)/sum(P)
 
 m=m+1
 }
realSignals

boxplot(realSignals, main = "relative peak intensity, 90 degree coupons")



## when n=23...and boundBox is on

plot(Fr,P,type="l", main = "periodogram coupon 23 (top only), signal = 0.6835686")
xline(findBand[1], col = "violetred1", lty=3)


plot(result$theta[boundBox], newCoupon[,3][boundBox], 
     col = color.scale(result$r),  pch = 16,
     main = paste0(n),
     xlab = "angle (radians)", ylab = "z axis")
xline(seq(3.7,4,by=0.0625), col = "grey", lty=3)


plot(result$theta, newCoupon[,3], 
     col = color.scale(result$r),  pch = 16,
     main = paste0(n),
     xlab = "angle (radians)", ylab = "z axis")

