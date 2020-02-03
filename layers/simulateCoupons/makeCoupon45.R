## 6/13/2019
##
##
##
## construct a test case for the 45 degree slice project
## ------------------------------------------------------


library(conicfit)
library(rgl)


setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("nlsFunctions.R")


numLayer <- seq(40,60,length.out = 10000)

rnoise <- 50


for(l in numLayer){
# calc points on a 45 degree line
x = seq(0,4000,by = l)
y = rep(0, length(x))

z=x
#noise <- runif(length(x),-l*noiseAmt,l*noiseAmt)
#z = x + noise


phi = 45*pi/180
a = 814*sec(phi) #major axis
b = 814 #minor axis

pts <- matrix(NA, nrow = 25, ncol = 3)
pts[,1:2] <- calculateEllipse(x[1],y[1],a,a, steps = 25, noiseFun=function(x) 
  (x+runif(1,-rnoise,rnoise)))
pts[,3] <- rep(z[1], 25)

ellipseStack <- pts

for(i in 2:length(x)){
  
pts <- matrix(NA, nrow = 25, ncol = 3)
pts[,1:2] <- calculateEllipse(x[i],y[i],a,b, steps = 25, randomDist = TRUE)
pts[,3] <- rep(z[i], 25)

ellipseStack <- rbind(ellipseStack, pts)

}

aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}

theta <- -runif(1,38,42)*pi/180 #rotate the coupon so it's almost upright to match the 
                    #usual orientation of the actual coupons. Works better with nls too

Ry <- aboutY(theta)

turnedStack <- ellipseStack %*% Ry 


while(class(nlsObj) == "try-error"){ 
  print(l)
  nlsObj <- try(nlsAxisFit(turnedStack))
}

nlsCoeff <- coef(nlsObj)

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/spacing40to60")

## store the old coupon coordinates, the "new" rotated coupon coords, and the nls coeff
## useful for generating surface plots and histograms for each coupon
nlsCoupon <- newCoupon(turnedStack, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])

oldCoupon <- turnedStack

print(round(l,3))

save(oldCoupon, nlsCoupon, nlsCoeff, file = paste0(round(l,3),"spacingSyn45.rda"))

}
  

