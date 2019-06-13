## 6/13/2019
##
##
##
## test for 45 degree slices by rotating coupon along a 
## sequence of z-axis rotations and slicing each one. 
## maximum number of zero bins in the transformed z-coord
## histogram is used  to indicate optimal rotation for finding layers
## ------------------------------------------------------


library(pracma)
library(useful)
library(rgl)
library(scales)


aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/45degreeSlice/45degreeData")

# setwd("C:/Users/barna/Documents/Coupons/nlsAxis/datasets")
# load("couponCov.rda")
# fortyFive <- which(couponCov$polarAngle==45)
# 
# setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")

numLayer <- seq(100,20,by=-5)

for(n in numLayer){
   #load(paste0("nlsCoupon",n,".rda"))
   #load("nlsCoupon5.rda")
  
  load(paste0(n,"spacingSyn45.rda"))
  # load("75stacks25points.rda")
  
  # plot old coupon
  # open3d()
  # par3d(cex=0.7)
  # plot3d(newCoupon[,1],
  #        newCoupon[,2],
  #        newCoupon[,3],
  #        type = "s", size = 0.45,
  #        xlab = " ", ylab = " ", zlab = " ")
  
  # center the coupon
  centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                        scale(newCoupon[,2], center = TRUE, scale = FALSE),
                        newCoupon[,3])
  

  
## ---------------------------------------------------------------
## start of for loop for angle Grid
  
  angleSeq <- seq(0, 2*pi, by = pi/20) #test out a sequence of theta
  
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
    # plot3d(RzyCoupon[,1],
    #        RzyCoupon[,2],
    #        RzyCoupon[,3],
    #        type = "s", size = 0.25,
    #        xlab = "x", ylab = "y", zlab = "z")

  
    h<-hist(RzyCoupon[,3], breaks = 300, plot = FALSE)
            #col = "darkgrey", main = paste0(theta))
    
    # h<-hist(RzyCoupon[,3], breaks = 300,
    # col = "darkgrey", main = paste0("z axis rotation (radians) = ",round(theta,3)))


    histSave[[i]] <- h$counts
    
    
    i = i + 1
    
  }
  

#seq(min(trunc(RzyCoupon[,3])), max(trunc(RzyCoupon[,3])), by = 1)
  
histSave <- data.frame(lapply(histSave, "length<-", max(lengths(histSave))))

names(histSave) <- round(angleSeq, digits=3)

plot(round(angleSeq, digits=3),apply(histSave, 2, function(x) table(x)[1]), type = "b", pch = 20, main = paste0(n))

}
which.max(apply(histSave, 2, function(x) table(x)[1]))


