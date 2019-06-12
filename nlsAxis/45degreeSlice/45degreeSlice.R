
library(pracma)
library(useful)
library(rgl)

fortyFive <- which(couponCov$polarAngle==45)

aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}

#look at a few 45 degree coupons
for(n in fortyFive[1]){
  setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")
  
  load(paste0("nlsCoupon",n,".rda"))
  
  load("nlsCoupontest.rda")
  
  # plot old coupon
  # open3d()
  # par3d(cex=0.7)
  # plot3d(oldCoupon[,1],
  #        oldCoupon[,2],
  #        oldCoupon[,3],
  #        type = "s", size = 0.45,
  #        xlab = " ", ylab = " ", zlab = " ")
  
  # center the coupon
  centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                        scale(newCoupon[,2], center = TRUE, scale = FALSE),
                        newCoupon[,3])
  
  
  angleSeq <- seq(0,pi, by = pi/30) #test out a sequence of theta
  
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
    #        type = "s", size = 0.45,
    #        xlab = "x", ylab = "y", zlab = "z")

  
    h<-hist(trunc(RzyCoupon[,3]), breaks = 130, col = "darkgrey")
    
  }
  
}




