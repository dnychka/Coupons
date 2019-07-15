## 7/04/2019
##
##
##
## working document for testing null hypothesis 
## distribution for 0 coupons
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice")
source("findPeaks.R")


setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/noLayering")

## ------------------------------------------------------


layerSpacing <- c(40,45,50,55,60)

#layerSpacing <- c(40,50,60)
numSims <- 10
harmonicSignalsNull <- rep(NA, (length(layerSpacing)*numSims))

m=1
for(k in layerSpacing){
  
  for(n in 1:numSims){
    
    print(n)
    
    load(paste0(k,"spacingRep",n,".rda"))
    
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
          ifelse(isTRUE(all.equal(fundFreq*j, Fr[i], tol = 0.05)), findHarmonic <- c(findHarmonic, P[i]),  NA)
          
        }
      }
      
    } else {
      fundFreq <- Fr[which.max(P[which(Fr<=0.15)])]
      ##find the harmonics
      for(j in 1:5){
        for(i in 1:length(Fr)){
          ifelse(isTRUE(all.equal(fundFreq*j, Fr[i], tol = 0.05)), findHarmonic <- c(findHarmonic, P[i]),  NA)
          
        }
      }
    }
    harmonicSignalsNull[m] <- sum(findHarmonic)/sum(P)
    
    m=m+1
  }
  
}

boxplot(harmonicSignalsNull)

saveRDS(harmonicSignalsNull,"zerosNullSignals.rds")


