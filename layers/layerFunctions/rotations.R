##
##
##
##
##
##
## functions to rotate the coupons
## ----------------------------------------

## --------------rotation matrices------------------------
aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}
## ------------------------------------------------------


## ------------- 45 rotation function------------------------
rotateFortyFive <- function(angle, nlsCoupon){
  
  
  centerCoupon <- cbind(scale(nlsCoupon[,1], center = TRUE, scale = FALSE),
                        scale(nlsCoupon[,2], center = TRUE, scale = FALSE),
                        nlsCoupon[,3])
  
  theta <- angle
    
  Rz <- aboutZ(theta)
    
  RzCoupon <- centerCoupon %*% Rz
    
  phi <- 45*pi/180 #45 degree tilt
    
  Ry <- aboutY(phi)
    
  RzyCoupon <- RzCoupon %*% Ry
    
    
  
  return(RzyCoupon)
  
}
## ------------------------------------------------------


## ------------- 0 rotation function------------------------
rotateZero <- function(angle, nlsCoupon){
  
  
  centerCoupon <- cbind(scale(nlsCoupon[,1], center = TRUE, scale = FALSE),
                        scale(nlsCoupon[,2], center = TRUE, scale = FALSE),
                        nlsCoupon[,3])
  
  theta <- angle
  
  Rz <- aboutZ(theta)
  
  RzCoupon <- centerCoupon %*% Rz
  
  return(RzCoupon)
  
}
## ------------------------------------------------------
