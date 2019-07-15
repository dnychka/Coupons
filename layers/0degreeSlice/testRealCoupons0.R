##
##
## working document for getting relative peak intenstiy
## from real 0 degree coupons
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice")
source("findPeaks.R")

setwd("C:/Users/barna/Documents/Coupons/datasets")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

## ------------------------------------------------------

Zero <- which(couponCov$polarAngle == 0)

numCoupons <- length(Zero)
realSignals <- rep(NA, ncol = length(Zero))

k=1

for(n in Zero){
  
  print(n)
  
  load(paste0("nlsCoupon", n, ".rda"))
  
  # center the coupon
  centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                        scale(newCoupon[,2], center = TRUE, scale = FALSE),
                        newCoupon[,3])
  
  
  h<-hist(~newCoupon[,3], w=5, plot = FALSE)
  
  histSave <- rbind(h$counts, h$mids)
  
  
  threshold <- 10
  influence <- 0
  
  xGrid <- seq(min(histSave[2,]), max(histSave[2,]), length.out = length(histSave[1,]))
  
  lowpass <- splint(histSave[2,], histSave[1,], xGrid, lambda = 50)
  
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
      if(fundFreq <= .15){break}
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
    fundFreq <- Fr[which.max(P[which(Fr<=0.15)])]
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
  
  
  realSignals[k] <- sum(findHarmonic)/sum(P)

 k=k+1

}

saveRDS(realSignals, "zerosRealSignals.rds")
