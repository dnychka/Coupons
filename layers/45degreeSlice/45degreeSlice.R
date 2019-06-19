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
library(useful)
library(rgl)
library(fields)



aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/45degreeSlice/45degreeData")

 
numLayer <- seq(100,5,by=-5)

for(n in numLayer){
  
  load(paste0(n,"spacingSyn45.rda"))
  
  # center the coupon
  centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                        scale(newCoupon[,2], center = TRUE, scale = FALSE),
                        newCoupon[,3])
  
  
## ---------------------------------------------------------------
## start of for loop for angle Grid
  
  #angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta
  angleDouble <- c(pi,2*pi) #for testing synthetics
  
  histSave <- list()
  
  i = 1
  
  for(theta in angleDouble){
    
    Rz <- aboutZ(theta)
    
    RzCoupon <- centerCoupon %*% Rz
    
    
    phi <- 45*pi/180 #45 degree tilt
    
    Ry <- aboutY(phi)
    
    RzyCoupon <- RzCoupon %*% Ry
    

    # open3d()
    # par3d(cex=0.7)
    # plot3d(RzyCoupon[,1],
    #        RzyCoupon[,2],
    #        RzyCoupon[,3],
    #        type = "s", size = 0.25,
    #        xlab = "x", ylab = "y", zlab = "z")

    h<-hist(~RzyCoupon[,3], w=n/10, plot = FALSE)
    #col = "darkgrey", main = paste0("z axis rotation (radians) = ",round(theta,3)))

    histSave[[i]] <- h$counts
    
    i = i + 1
    
  }


  for(m in 1:2){
  
  frDt <- fft(histSave[[m]]-mean(histSave[[m]]))
  
  Dtlen <- length(histSave[[1]])
  
  Fr <- 1:Dtlen/Dtlen
  P <- Mod(2*frDt/Dtlen)^2
  
  plot(Fr[1:(Dtlen/2)],P[1:(Dtlen/2)],  type = "l", main = paste0("layers",n,"_",m))
  firstHarmonic <- which.max(P[1:(Dtlen/2)])
  xline(Fr[firstHarmonic], col = "grey", lty = 3)
  
  #h$breaks
  
  Fr[firstHarmonic] 
  #interp: every 10 bins we have a full cycle. Each binwidth is 10,
  #and the coupon has a layer spacing of 100 = checks out
  }
}

