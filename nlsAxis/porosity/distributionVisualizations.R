## 
## 3 June 2019
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## Create summary plots(boxplots, histograms) nls fit
##
##------------------------------------------------------------------------------------------------------------------------

library(fields)
library(plyr)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/datasets")

load("radiusIterations.rda")
load("couponCov.rda")

medianWithoutNA<-function(x) {
  median(x[which(!is.na(x))])
}


##---------------------------------------------------------
## all group (58 coupons) diagnostic plots
##---------------------------------------------------------

bplot(t(radiusLastIter), by = couponCov$polarAngle,
      boxwex = 0.5,
      pch = 20,
      main = "radial distance from center axis, nls fit",
      ylab = "radial distance",
      xlab = "polar angle")


vertMedians <- apply(radiusLastIter[,c(Zero,fortyFive)],2, medianWithoutNA)

hist(vertMedians, breaks = 20, freq = FALSE)

hist(t(radiusLastIter[,c(Zero,fortyFive)]), breaks = 50, freq=FALSE)
xline(median(vertMedians), col = "tomato", lty=3)

bplot(vertMedians, 
      by = couponCov$polarAngle[c(Zero,fortyFive)],
      boxwex = 0.5, pch = 20,
      frame = F)


##--------------
## polar angle 0 
##--------------
Zero <- which(couponCov$polarAngle==0)

## histogram, final radial distance
hist(as.vector(t(radiusLastIter[,Zero])), breaks = seq(0,1400, by = 25), freq = FALSE,
     col = "grey40", main = "", xlab = "radial distance")
mtext("radial pore distance to nls fit axis", side=3, adj=0, line=1.6, cex=1.4, font=1)
mtext("polar angle 0 ", side=3, adj=0, line=0.3, cex=1, font=1)

xline(median(apply(radiusLastIter[,Zero],2, medianWithoutNA)), col = "tomato")
xline(median(apply(radiusLastIter[,fortyFive],2, medianWithoutNA)), col = "cornflowerblue")

## boxplots before/after median iteration
par(mfrow=c(1,2))
bplot(radiusFirstIter[,Zero], pch = 20, 
      main = "first iteration of targetRadius, polar angle 0",
      ylim = c(0,1320))

yline(median(apply(radiusFirstIter[,Zero],2, medianWithoutNA)), col = "tomato")

bplot(radiusLastIter[,Zero], pch = 20, 
      main = "last iteration of targetRadius, polar angle 0",
      ylim = c(0,1320))

yline(median(apply(radiusLastIter[,Zero],2, medianWithoutNA)), col = "tomato")




##--------------
## polar angle 45 
##--------------
fortyFive <- which(couponCov$polarAngle==45)

par(mfrow=c(1,1))

## histogram, final radial distance
hist(as.vector(t(radiusLastIter[,fortyFive])), breaks = seq(0,1400, by = 25), freq = FALSE,
     col = "grey40", main = "", xlab = "radial distance")
mtext("radial pore distance to nls fit axis", side=3, adj=0, line=1.6, cex=1.4, font=1)
mtext("polar angle 45 ", side=3, adj=0, line=0.3, cex=1, font=1)


## boxplots before/after median iteration
par(mfrow=c(1,2))
bplot(radiusFirstIter[,fortyFive], pch = 20, 
      main = "first iteration of targetRadius, polar angle 45",
      ylim = c(0,1320))

yline(median(apply(radiusFirstIter[,fortyFive],2, medianWithoutNA)), col = "tomato")

bplot(radiusLastIter[,fortyFive], pch = 20, 
      main = "last iteration of targetRadius, polar angle 45",
      ylim = c(0,1320))

yline(median(apply(radiusLastIter[,fortyFive],2, medianWithoutNA)), col = "tomato")


##--------------
## polar angle 90 
##--------------
Ninety <- which(couponCov$polarAngle==90)

par(mfrow=c(1,1))

##histogram, final radial distance
hist(as.vector(t(radiusLastIter[,Ninety])), breaks = seq(0,1400, by = 25), freq = FALSE,
     col = "grey40", main = "", xlab = "radial distance")
mtext("radial pore distance to nls fit axis", side=3, adj=0, line=1.6, cex=1.4, font=1)
mtext("polar angle 90 ", side=3, adj=0, line=0.3, cex=1, font=1)

## boxplots before/after median iteration
par(mfrow=c(1,2))
bplot(radiusFirstIter[,Ninety], pch = 20, 
      main = "first iteration of targetRadius, polar angle 90",
      ylim = c(0,1320))

yline(median(apply(radiusFirstIter[,Ninety],2, medianWithoutNA)), col = "tomato")

bplot(radiusLastIter[,Ninety], pch = 20, 
      main = "last iteration of targetRadius, polar angle 90",
      ylim = c(0,1320))

yline(median(apply(radiusLastIter[,Ninety],2, medianWithoutNA)), col = "tomato")






##---------------------------------------------------------
## individual coupon histograms
##---------------------------------------------------------

n = 4

hist(as.vector(t(radiusLastIter[,4])), breaks = seq(0,1400, by = 25), freq = FALSE,
     col = "grey40", main = "", xlab = "radial distance")
mtext("radial pore distance to nls fit axis", side=3, adj=0, line=1.6, cex=1.4, font=1)
mtext(paste("coupon",  n, "polar angle 0 "), side=3, adj=0, line=0.3, cex=1, font=1)




for(i in Ninety){
  plot(poreData[[i]]$comZ, poreData[[i]]$comX,
       main = paste0(i))
}
  




