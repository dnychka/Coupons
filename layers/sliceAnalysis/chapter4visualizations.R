##
##
##
##
##
## make pretty histograms / periodograms for ch 4
##
## --------------------------------------------------------

library(pracma)
library(rgl)
library(fields)
library(FSA)
library(scales)


## --------------------------------------------------
## make histogram for 45 degree coupon with spline
## --------------------------------------------------
setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("findPeaks.R")
source("45degreeFunctions.R")
setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/spacing40to60")

load("59.866spacingSyn45.rda") #Branden said layer spacing is most likely 50 microns

centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                      scale(newCoupon[,2], center = TRUE, scale = FALSE),
                      newCoupon[,3]*-1)


theta = 0

Rz <- aboutZ(theta)

RzCoupon <- centerCoupon %*% Rz

phi <- 45*pi/180 #45 degree tilt

Ry <- aboutY(phi)

RzyCoupon <- RzCoupon %*% Ry

windowsFonts(A = windowsFont("Times New Roman"))
hCoarse<-hist(~RzyCoupon[,3], w=5,
              col = "grey30", main = "", xlab = "z coordinates", cex.lab = 1.2, family = "A")
mtext("Histogram of pore positions for layered 45-degree coupon", side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")
xGridC <- seq(min(hCoarse$mids), max(hCoarse$mids), length.out = length(hCoarse$counts))
lowpassC <- splint(hCoarse$mids, hCoarse$counts, xGridC, lambda = 100)
lines(xGridC, lowpassC, col = "firebrick1", lwd=3)


## --------------------------------------------------
## make histogram for 0 degree coupon with spline
## --------------------------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/spacing40to60")

load("59.866spacingSyn0.rda")

centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                      scale(newCoupon[,2], center = TRUE, scale = FALSE),
                      newCoupon[,3])

centerCoupon <- cbind(scale(nlsCoupon[,1], center = TRUE, scale = FALSE),
                      scale(nlsCoupon[,2], center = TRUE, scale = FALSE),
                      nlsCoupon[,3])


h<-hist(~newCoupon[,3], w=5, plot = FALSE)
h<-hist(~nlsCoupon[,3], w=5, plot = FALSE)


histSave <- rbind(h$counts, h$mids)

xGrid <- seq(min(histSave[2,]), max(histSave[2,]), length.out = length(histSave[1,]))

lowpass <- splint(histSave[2,], histSave[1,], xGrid, lambda = 50)

windowsFonts(A = windowsFont("Times New Roman"))
plot(h, col = "grey30", main = "", xlab = "z coordinates", cex.lab = 1.2, family = "A", ylim = c(0, max(hCoarse$counts)))
mtext("Histogram of pore positions for layered 0-degree coupon", side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")
lines(xGrid, lowpass, col = "firebrick1", lwd=3)




## --------------------------------------------------
## plot a signal periodogram and non signaled periodogram
## --------------------------------------------------

##
## non signaled
## ------------
load("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/noLayering/60spacingRep10.rda")

centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                      scale(newCoupon[,2], center = TRUE, scale = FALSE),
                      newCoupon[,3])

h<-hist(~newCoupon[,3], w=5, plot = FALSE)

histSave <- rbind(h$counts, h$mids)

xGrid <- seq(min(histSave[2,]), max(histSave[2,]), length.out = length(histSave[1,]))

lowpass <- splint(histSave[2,], histSave[1,], xGrid, lambda = 50)

ff <- fft(h$counts-lowpass)
len <- length(h$mids)
Fr <- (1:len/len)[1:(len/2)]
P <- (Mod(2*ff/len)^2)[1:(len/2)]
plot(Fr, P, type = "l", main="", 
     ylim = c(0,1.05), xlab = "frequency", ylab = "strength", family = "A", cex.lab = 1.2)
mtext("Periodogram for Non-layered Coupon", side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")


##
## signaled
## --------
load("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/spacing40to60/59.597spacingSyn0.rda")

centerCoupon <- cbind(scale(nlsCoupon[,1], center = TRUE, scale = FALSE),
                      scale(nlsCoupon[,2], center = TRUE, scale = FALSE),
                      nlsCoupon[,3])

h<-hist(~nlsCoupon[,3], w=5, plot = FALSE)

histSave <- rbind(h$counts, h$mids)

xGrid <- seq(min(histSave[2,]), max(histSave[2,]), length.out = length(histSave[1,]))

lowpass <- splint(histSave[2,], histSave[1,], xGrid, lambda = 50)

ff <- fft(h$counts-lowpass)
len <- length(h$mids)
Fr <- (1:len/len)[1:(len/2)]
P <- (Mod(2*ff/len)^2)[1:(len/2)]
plot(Fr, P, type = "l", main="", 
      xlab = "frequency", ylab = "strength", family = "A", cex.lab = 1.2)
mtext("Periodogram for Layered Coupon", side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")


##
## add in harmonic finding bands
## -----------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/0degreeSlice")
source("findPeaks.R")

plot(Fr, P, type = "l", main="", 
     xlab = "frequency", ylab = "strength", family = "A", cex.lab = 1.2)
mtext("Periodogram with Harmonics Highlighted", side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")


lag = length(Fr)/10
threshold <- 10
influence <- 0

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
  
  xline(findBand, col = "cornflowerblue", lty=3)
  
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

xline(c(0.5,1,2,3,4,5)*fundFreq, lty = 2, col = "cornflowerblue")


## --------------------------------
## show plot with higher harmonic than ff
## --------------------------------
load("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/spacing40to60/58.658spacingSyn0.rda")

centerCoupon <- cbind(scale(nlsCoupon[,1], center = TRUE, scale = FALSE),
                      scale(nlsCoupon[,2], center = TRUE, scale = FALSE),
                      nlsCoupon[,3])

h<-hist(~nlsCoupon[,3], w=5, plot = FALSE)

histSave <- rbind(h$counts, h$mids)

xGrid <- seq(min(histSave[2,]), max(histSave[2,]), length.out = length(histSave[1,]))

lowpass <- splint(histSave[2,], histSave[1,], xGrid, lambda = 50)

ff <- fft(h$counts-lowpass)
len <- length(h$mids)
Fr <- (1:len/len)[1:(len/2)]
P <- (Mod(2*ff/len)^2)[1:(len/2)]
plot(Fr, P, type = "l", main="", 
     xlab = "frequency", ylab = "strength", family = "A", cex.lab = 1.2)
mtext("Simulated Layered Coupon, 60 micron Spacing", side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")
xline(0.087, lty = 2, lwd = 3, col = alpha("cornflowerblue",0.4))
lines(Fr, P, lwd = 1.5)
