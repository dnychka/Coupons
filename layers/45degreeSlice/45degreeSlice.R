## 6/13/2019
##
##
##
## test for 45 degree slices by rotating coupon along a 
## sequence of z-axis rotations and slicing each one. 
## maximum number of zero bins in the transformed z-coord
## histogram is used  to indicate optimal rotation for finding layers
##
## in synthetic coupons, rotational angle is always pi or 2pi
## 
## updated to add spectral decomposition analysis
## ------------------------------------------------------


library(pracma)
library(rgl)
library(fields)
library(FSA)

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice")
source("threshold.R")

aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData")

numLayer <- seq(100,5,by=-5)

angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta
angleDouble <- c(pi,2*pi) #for testing synthetics


j=1
relSignal <- matrix(NA, nrow = length(angleSeq), ncol = length(numLayer))



for(n in numLayer){
  
  print(n)
  
  #load(paste0(n,"spacingSyn45.rda"))
  load(paste0(n,"spacingNoLayer.rda"))
  
  # center the coupon
  centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                        scale(newCoupon[,2], center = TRUE, scale = FALSE),
                        newCoupon[,3])
  
  
## ---------------------------------------------------------------
## start of for loop for angle Grid

  
  histSave <- list()
  
  i = 1
  
  for(theta in angleSeq){
    
    Rz <- aboutZ(theta)
    
    RzCoupon <- centerCoupon %*% Rz
    
    
    phi <- 45*pi/180 #45 degree tilt
    
    Ry <- aboutY(phi)
    
    RzyCoupon <- RzCoupon %*% Ry
    

    # open3d()
    # par3d(cex=0.7)
    # points3d(RzyCoupon[,1],
    #        RzyCoupon[,2],
    #        RzyCoupon[,3],
    #        size = 2,
    #        xlab = "x", ylab = "y", zlab = "z")

    h<-hist(~RzyCoupon[,3], w=n/10, plot = FALSE)
    #col = "darkgrey", main = paste0("z axis rotation (radians) = ",round(theta,3)))
     
    # h<-hist(~RzyCoupon[,3], w=n/10,
    # col = "darkgrey", main = paste0("z axis rotation (radians) = ",round(theta,3), " ", i), ylim = c(0,25))

    histSave[[i]] <- rbind(h$counts, h$mids)
    
    i = i + 1
    
  }

  
##
## find the periodogram that has meaningful frequency spikes
##
  
findTheta <- rep(NA,length(angleSeq))
frDt <- list()
Dtlen <- rep(NA, length(angleSeq))

high <- rep(NA, length(angleSeq))

threshold <- 10
influence <- 0

  for(m in 1:length(angleSeq)){
    
  xGrid <- seq(min(histSave[[m]][2,]), max(histSave[[m]][2,]), length.out = length(histSave[[m]][1,]))
  
  lowpass <- splint(histSave[[m]][2,], histSave[[m]][1,], xGrid, lambda = 50)
  
  frDt[[m]] <- fft(histSave[[m]][1,]-lowpass)
  
  Dtlen[m] <- length(histSave[[1]][1,])
  
  Fr <- (1:Dtlen[m]/Dtlen[m])[1:(Dtlen[m]/2)]
  P <- (Mod(2*frDt[[m]]/Dtlen[m])^2)[1:(Dtlen[m]/2)]
  
  # plot(Fr,  P,  type = "l", main = paste0("layers spaced ",n," apart", " ", m),
  #      ylab = "strength", xlab = "layer frequency", ylim = c(0,1.5))
   
  lag <- length(Fr)/10
  
  peak <- findFreq(P,lag,threshold,influence)
  
  high[m] <- (max(P[which(peak$signals==1)]))
  
  #relSignal[m,j] <- sum(P[which(peak$signals==1)])/sum(P)
  
  }

# plot(high, pch = 20)

ind <- which.max(high)

Fr <- (1:Dtlen[ind]/Dtlen[ind])[1:(Dtlen[ind]/2)]
P <- (Mod(2*frDt[[ind]]/Dtlen[ind])^2)[1:(Dtlen[ind]/2)]

# plot(Fr,  P,  type = "l", main = paste0("layers spaced ",n," apart"),
#      ylab = "strength", xlab = "layer frequency")



j = j+1

}


matplot(t(relSignal), pch=20, col = "black", xaxt = "n",
        xlab = "layer spacing", ylab = "sum(peaks)/sum(all points)", main = "relative peak strength")
points(1:20, apply(relSignal, 2, max), pch = 16, col = "mediumturquoise")
axis(1,at=c(1,5,10,15,20), labels = c("100", "80", "55", "30", "5"))



