## 6/25/2019
##
##
##
## working document for testing null hypothesis 
## distribution
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("threshold.R")
source("45degreeFunctions.R")

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/nullHyp")

## ------------------------------------------------------



angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta

numSims <- 10
tempSignals <- matrix(NA, nrow = length(angleSeq), ncol = numSims)
relSignals <- NULL


layerSpacing <- seq(40,60,by=2.5)

for(i in layerSpacing){

  for(n in 1:numSims){
    
    print(n)
  
    load(paste0(i,"spacingRep",n,".rda"))
    
    # center the coupon
    centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                          scale(newCoupon[,2], center = TRUE, scale = FALSE),
                          newCoupon[,3])
    
  
    histSave <- FFrotate(angleSeq, centerCoupon)
  
    tempSignals[,n] <-  makePeriodogram(angleSeq, histSave)
    
  }
  
  tempSignals <- rbind(rep(i, length(numSims)), tempSignals) #name the columns for fbplot
  
  relSignals <- cbind(relSignals, tempSignals)
  
  tempSignals <- tempSignals[-1,]
 
} 

saveRDS(relSignals, "relSignalsHarmn.rda")

## ------------------------------------------
## examining distribution by pore density grouping
## ------------------------------------------

bplot((relSignals[-1,]), by = relSignals[1,])

b<-bplot(apply(relSignals[-1,], 2, max), by = relSignals[1,]) 

# some messy t-tests
t.test(b$stats[,1], b$stats[,2])
t.test(b$stats[,1], b$stats[,3])
t.test(b$stats[,1], b$stats[,4])
t.test(b$stats[,1], b$stats[,5])
t.test(b$stats[,1], b$stats[,6])

