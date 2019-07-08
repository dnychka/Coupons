## 7/04/2019
##
##
##
## test out our summary statistics in the synthetic layered
## coupons for layered 0-0
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("threshold.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData")

## ------------------------------------------------------

layerSpacing <- seq(40,60,length.out = 150)

harmonicSignals <- rep(NA, ncol = length(layerSpacing))

n=59.060402684538
n=59.8657718
n=48.32215
n=51.27517
n=46.1745
n=59.329
n=59.195
n=43.087

k=1
for(n in layerSpacing){
  
  print(round(n,3))
  
  load(paste0(round(n,3),"spacingSyn0.rda"))
  
  #load("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/nullHyp/50spacingRep5.rda")
  
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
  
  
  
  Ppeaks <- P[ind]
  Frpeaks <- Fr[ind]
  
  
   #fundFreq <- Frpeaks[which.max(Ppeaks)]
  
   
   
   #fundFreq <- Fr[which.max(P[which(Fr<=0.15)])]
     
   findHarmonic <- vector()


   origFreq <- Frpeaks[which.max(Ppeaks)]
   
   for(j in c(1,1/2,1/3,1/4,1/5)){
     fundFreq <- origFreq*j
     if(fundFreq <= .15){break}
   }
   
   xline(fundFreq*c(1,2,3,4,5), col = "tomato", lty=3)
   
  ##check for lowest freq (i.e. fundamental frequency)

  # for(j in c(1/5,1/4,1/3,1/2,2,3,4,5)){
  #   for(i in 1:length(Frpeaks)){
  #     ifelse(all.equal(fundFreq*j, Frpeaks[i], tol = 0.05)==TRUE,  fundFreq <- c(fundFreq,Frpeaks[i]),  NA)
  #   }
  # }

  # xline(Frpeaks, col="steelblue", lty=3)
  
 # fundFreq <- min(fundFreq)
  
 # xline(fundFreq*c(1,2,3,4,5), col = "grey", lty=3)
  
  # xline(min(Frpeaks), col = "cornflowerblue", lty=3)
  # 
  # fundFreq <- min(fundFreq)
  # 
  # xline(fundFreq*c(1,2,3,4), col = "violetred1", lty=3)
  # 
  
  
  

  ##find the harmonics
  for(j in 1:5){
    for(i in 1:length(Fr)){
      ifelse(isTRUE(all.equal(fundFreq*j, Fr[i], tol = 0.05)), findHarmonic <- c(findHarmonic, P[i]),  NA)
      
    }
  }
  
  #xline(Frpeaks[c(4,5,6,7,9,10,11,12)], col="grey", lty=3)
  
  #findHarmonic <- c(findHarmonic, Ppeaks[i])
  
  harmonicSignals[k] <- sum(findHarmonic)/sum(P)

  k=k+1
}

harmonicSignals

boxplot(harmonicSignals, main = "relative peak strength, synthetic layers")

#saveRDS(harmonicSignals, "synthetic0layers.rds")
#