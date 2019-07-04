

library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/datasets")
load("couponCov.rda")

fortyFive <- which(couponCov$polarAngle==45)



## ------------------------------------------------------------------------
## exploring pore density
## ------------------------------------------------------------------------

s <- 20

dense <- NULL

for(s in 45:55){
numPore <- length(seq(0,4000,by = s))
dense[s] <- numPore*25 #find density of pores in layered version of coupon
}

plot( dense, ylim = c(0,3000))


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData/cropped")

len <- rep(length(fortyFive))
i=1
for(n in fortyFive){
  load(paste0("nlsCoupon", n, ".rda"))
  len[i]<-length(newCoupon[,1])
  i=i+1
}


plot(len)

  
  
## ------------------------------------------------------------------------
## show binning / lowpass filter issues
## ------------------------------------------------------------------------
  
setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("threshold.R")
source("45degreeFunctions.R")
setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData")

load("55spacingSyn45.rda") #Branden said layer spacing is most likely 50 microns

centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                      scale(newCoupon[,2], center = TRUE, scale = FALSE),
                      newCoupon[,3])


  theta = 0
  
  Rz <- aboutZ(theta)
  
  RzCoupon <- centerCoupon %*% Rz
  
  phi <- 45*pi/180 #45 degree tilt
  
  Ry <- aboutY(phi)
  
  RzyCoupon <- RzCoupon %*% Ry
  
## test out histogram bins  
  par(mfrow=c(1,2))
  
  hCoarse<-hist(~RzyCoupon[,3], w=5,
  col = "grey30", main = "coarse binning")
  xGridC <- seq(min(hCoarse$mids), max(hCoarse$mids), length.out = length(hCoarse$counts))
  lowpassC <- splint(hCoarse$mids, hCoarse$counts, xGridC, lambda = 100)
  lines(xGridC, lowpassC, col = "firebrick1", lwd=2)
  
  
  # hFine<-hist(~RzyCoupon[,3], w=0.5,
  #               col = "grey30", main = "fine binning", ylim = c(0,max(hCoarse$counts)))
  # xGridF <- seq(min(hFine$mids), max(hFine$mids), length.out = length(hFine$counts))
  # lowpassF <- splint(hFine$mids, hFine$counts, xGridF, lambda = 100)
  # lines(xGridF, lowpassF, col = "firebrick1", lwd=2)

  
  
  ffCoarse <- fft(hCoarse$counts-lowpassC)
  lenCoarse <- length(hCoarse$mids)
  Fr <- (1:lenCoarse/lenCoarse)[1:(lenCoarse/2)]
  P <- (Mod(2*ffCoarse/lenCoarse)^2)[1:(lenCoarse/2)]
  plot(Fr, P, type = "l", main="coarse binning: periodogram")
  
  

  
  lag = length(Fr)/10
  threshold <- 10
  influence <- 0
  
  peak <- findFreq(P,lag,threshold,influence)
  
  ind <- which(peak$signals==1)
 
  Ppeaks <- P[ind]
  Frpeaks <- Fr[ind]
  
  fundFreq <- Frpeaks[which.max(Ppeaks)]
  findHarmonic <- vector()
  
  for(j in 1:5){
    for(i in 1:length(Fr)){
    ifelse(all.equal(fundFreq*j, Fr[i], tol = 0.01)==TRUE,  findHarmonic <- c(findHarmonic, P[i]),  NA)
    }
  }
  
  sum(findHarmonic)/sum(P)

  #xline(c(0.5,1,2,3,4,5)*fundFreq, col = "grey", lty = 3)
 
  
  
  
  # #get an idea of how many frequencies we're looking for
  # t <-  hist(Frpeaks, plot = FALSE)
  # numH <- length(which(t$density != 0))
  # 
  # tempH <- rep(NA, length(Frpeaks))
  # findHarmonic <- rep(NA, numH-1)
  # 
  # for(j in 2:numH){
  #   for(i in 1:length(Frpeaks)){
  #     tempH[i] <- all.equal(fundFreq*j, Frpeaks[i], tol = 0.01)
  #   }
  #   
  #   #check for double peak
  #   
  #   dbl <- which(tempH==TRUE)
  #   print(dbl)
  #   
  #   
  #   if( length(dbl) > 1){
  #     closePeak <- which.min(fundFreq*j - Frpeaks[dbl])
  #     findHarmonic[j-1] <- dbl[closePeak]
  #   } else if (length(dbl)==1){
  #     findHarmonic[j-1] <- dbl
  #   }
  # }
  # 
  # findHarmonic
  
  