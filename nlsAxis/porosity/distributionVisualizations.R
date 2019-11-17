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

setwd("C:/Users/barna/Documents/Coupons/datasets")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData")
load("radiusIterations.rda")

medianWithoutNA<-function(x) {
  median(x[which(!is.na(x))])
}

Zero <- which(couponCov$polarAngle==0)
fortyFive <- which(couponCov$polarAngle==45)




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

h<-hist(t(radiusLastIter[,c(Zero,fortyFive)]), breaks = 50, freq=FALSE,
        col = "grey40", main = "",
        xlab = "radius (micrometers)")
mtext("pore radii, center axis found by nonlinear least squares", side=3, adj=0, line=1.6, cex=1.4, font=1)
mtext("polar angles 0 and 45", side=3, adj=0, line=0.3, cex=1, font=1)
xline(median(vertMedians), col = "tomato", lty=2, lwd=2)
legend("bottomleft", "median pore radius", col = "tomato", bty= "n", lty=2, lwd = 2, inset = 0.07)

bplot(vertMedians, 
      by = couponCov$polarAngle[c(Zero,fortyFive)]*-1,
      boxwex = 0.5, pch = 16,
      frame = F,horizontal = TRUE,
      axes=FALSE, xlab = "micrometers",cex.lab = 1.2,
      main = "")
mtext("median pore radius by polar angle", side=3, adj=0, line=1.6, cex=1.4, font=1)
axis(1, at=seq(830,890, by=10), labels=seq(830,890, by=10))
axis(2, at=c(1,2), labels = c("forty-five", "zero"), tick = FALSE, cex.axis=1.2)

xline(median(apply(radiusLastIter[,c(Zero)],2, medianWithoutNA)), col = "cornflowerblue", lty=3)
xline(median(apply(radiusLastIter[,c(fortyFive)],2, medianWithoutNA)), col = "tomato", lty=3)

kruskal.test(apply(radiusLastIter[,c(Zero,fortyFive)],2,medianWithoutNA), c(rep(1,16), rep(2,24)))

unzip <- as.vector(cbind(radiusLastIter[,Zero]))

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
  




