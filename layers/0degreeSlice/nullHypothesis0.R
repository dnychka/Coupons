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

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("threshold.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/nullHyp")

## ------------------------------------------------------


layerSpacing <- c(40,45,50,55,60)

layerSpacing <- c(40,50,60)
numSims <- 5
harmonicSignalsNull <- vector()


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
    
    plot(Fr,P,type="l", main = paste0(n))
    
    lag <- length(Fr)/10
    
    peak <- findFreq(P,lag,threshold,influence)
    
    ind <- which(peak$signals==1)
    
    Ppeaks <- P[ind]
    Frpeaks <- Fr[ind]
    
    
    
    findHarmonic <- vector()
  
    origFreq <- Frpeaks[which.max(Ppeaks)]
    
    for(j in c(1,1/2,1/3,1/4,1/5)){
      fundFreq <- origFreq*j
      if(fundFreq <= .15){break}
    }
    
    xline(fundFreq*c(1,2,3,4,5), col = "tomato", lty=3)
    
    ##find the harmonics
    for(j in 1:5){
      for(i in 1:length(Fr)){
        ifelse(all.equal(fundFreq*j, Fr[i], tol = 0.05)==TRUE, findHarmonic <- c(findHarmonic, P[i]),  NA)
        
      }
    }
    
    harmonicSignalsNull <- c(harmonicSignalsNull, sum(findHarmonic)/sum(P))
   
  }

} 

#saveRDS(harmonicSignalsNull, "harmonicSignalsNull.rds")

boxplot(harmonicSignalsNull)



