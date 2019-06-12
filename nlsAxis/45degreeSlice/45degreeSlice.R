
library(pracma)
library(useful)
library(rgl)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/datasets")
load("radiusIterations.rda") #to get median radius from nls fit
load("couponCov.rda") #to get covariates


fortyFive <- which(couponCov$polarAngle==45)

#look at a few 45 degree coupons
for(n in fortyFive[1]){
  setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")
  
  load(paste0("nlsCoupon",n,".rda"))
  
  load("nlsCoupontest.rda")
  
  # plot old coupon
  open3d()
  par3d(cex=0.7)
  plot3d(oldCoupon[,1],
         oldCoupon[,2],
         oldCoupon[,3],
         type = "s", size = 0.45,
         xlab = " ", ylab = " ", zlab = " ")
  
  # center the coupon
  centerCoupon <- cbind(scale(newCoupon[,1], center = TRUE, scale = FALSE),
                        scale(newCoupon[,2], center = TRUE, scale = FALSE),
                        newCoupon[,3])
  
  open3d()
  par3d(cex=0.7)
  plot3d(newCoupon[,1],
         newCoupon[,2],
         newCoupon[,3],
         type = "s", size = 0.45,
         xlab = " ", ylab = " ", zlab = " ")
  
  open3d()
  par3d(cex=0.7)
  plot3d(centerCoupon[,1],
         centerCoupon[,2],
         centerCoupon[,3],
         type = "s", size = 0.45,
         xlab = " ", ylab = " ", zlab = " ")
  
  theta <- pi #the angle that we iterate over
  
  angleSeq <- seq(0,pi, by = pi/30) #test out a sequence of theta
  
  for(theta in angleSeq){
    
    Rxy <- matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )
    
    RxyCoupon <- centerCoupon %*% Rxy
    
    
    phi <- 45*pi/180 #45 degree tilt
    
    Ry <- aboutY(phi)
    
    RzCoupon <- RxyCoupon %*% Ry
    
    
    #get new center axis for tilted coupon
    axisMat <- cbind(rep(0,length(centerCoupon[,1])), rep(0,length(centerCoupon[,1])), centerCoupon[,3])
    
    newAxis <- axisMat %*% Rz
    
    # open3d()
    # par3d(cex=0.7)
    # plot3d(RzCoupon[,1],
    #        RzCoupon[,2],
    #        RzCoupon[,3],
    #        type = "s", size = 0.45,
    #        xlab = "x", ylab = "y", zlab = "z")
    # points3d(newAxis[,1],
    #          newAxis[,2],
    #          newAxis[,3],
    #          size = 3, col = "tomato")


    
    # find angle to origin (newAxis) and horizontal distance
    
    centerForPolar <- RzCoupon[,1:2]-newAxis[,1:2]
    
    polar <- useful::cart2pol(centerForPolar[,1], centerForPolar[,2])
    
    
    #plot(polar$theta, polar$r, ylim = c(0,1200), main = paste0(n))
    
    medianWithoutNA<-function(x) {
      median(x[which(!is.na(x))])
    }
    
    b = medianWithoutNA(radiusLastIter) #instead of 1000 use median radius from algined coupon
    a = medianWithoutNA(radiusLastIter)*sec(phi)
    
    ellipseR <- function(theta){a*b/sqrt( a^2*sin(theta)^2 + b^2*cos(theta)^2)}
    
    #points(polar$theta, ellipseR(polar$theta), col = "red")
    
    
    h<-hist(trunc(RzCoupon[,3]), breaks = 130, col = "darkgrey")
    
  }
  
}


open3d()
par3d(cex=0.7)
plot3d(RzCoupon[,1],
       RzCoupon[,2],
       RzCoupon[,3],
       type = "s", size = 0.45,
       xlab = "x", ylab = "y", zlab = "z")


