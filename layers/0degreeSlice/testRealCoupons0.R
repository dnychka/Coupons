##
##
## working document for getting relative peak intenstiy
## from real 0 degree coupons
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("threshold.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/datasets")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/nullHyp")
harmonicSignalsNull <- readRDS("harmonicSignalsNull.rds")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData/cropped")

## ------------------------------------------------------

Zero <- which(couponCov$polarAngle == 0)

numCoupons <- length(Zero)
realSignals <- vector()


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
  
  plot(Fr,P,type="l", main = paste0(n))
  
  lag <- length(Fr)/10
  
  peak <- findFreq(P,lag,threshold,influence)
  
  ind <- which(peak$signals==1)
  
  Ppeaks <- P[ind]
  Frpeaks <- Fr[ind]
  
  # fundFreq <- Frpeaks[which.max(Ppeaks)]
  # findHarmonic <- vector()
  # 
  # 
  # ##check for lowest freq (i.e. fundamental frequency)
  # 
  # for(j in c(1/4,1/3,1/2,2,3,4,5)){
  #   for(i in 1:length(Frpeaks)){
  #     ifelse(all.equal(fundFreq*j, Frpeaks[i], tol = 0.01)==TRUE,  fundFreq <- c(fundFreq,Frpeaks[i]),  NA)
  #   }
  # }
  # 
  # xline(fundFreq, col = "grey", lty=3)
  # 
  # fundFreq <- min(fundFreq)
  # 
  # xline(fundFreq, col = "violetred1", lty=3)
  # 
  # ##find the harmonics
  # for(j in 1:5){
  #   for(i in 1:length(Fr)){
  #     ifelse(all.equal(fundFreq*j, Fr[i], tol = 0.01)==TRUE,  findHarmonic <- c(findHarmonic, P[i]),  NA)
  #   }
  # }
  
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
  

realSignals <- c(realSignals, sum(findHarmonic)/sum(P))

}
layeredCoupons <- harmonicSignals
nullCoupons <- c(harmonicSignalsNull, rep(NA, length(layeredCoupons)-length(harmonicSignalsNull)))
realCoupons <- c(realSignals, rep(NA, (length(nullCoupons)-length(realSignals))))

dat <- cbind(nullCoupons, realCoupons)

#dev.off()
boxplot(dat, main = "relative peak strength, 0 degree coupons",
        pch = 16, boxwex = 0.7)

dat <- cbind(dat, layeredCoupons)

boxplot(dat, main = "relative peak strength, 0 degree coupons",
        pch = 16, boxwex = 0.7)

t.test(nullCoupons, realCoupons)
t.test(harmonicSignals, realCoupons)
