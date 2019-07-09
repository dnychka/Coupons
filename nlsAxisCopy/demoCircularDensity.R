




library(rgl)
library(circular)
library(fields)
library(useful)


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")

load("nlsCoupon5.rda")

for(i in Zero){
  load(paste0("nlsCoupon",i,".rda"))

polar <- cart2pol(newCoupon[,1], newCoupon[,2], degrees = TRUE)

circ <- circular(trunc(polar$theta), units = "degrees")

plot(circ, stack = FALSE, shrink = 1.3, cex = 1.03, bins=24, main = paste0(i))
points(circ, rotation = "clock", stack = TRUE)
}


for(i in fortyFive){
  load(paste0("nlsCoupon",i,".rda"))
  
  polar <- cart2pol(newCoupon[,1], newCoupon[,2], degrees = TRUE)
  
  circ <- circular(trunc(polar$theta), units = "degrees")
  
  plot(circ, stack = FALSE, shrink = 1.3, cex = 1.03, bins=24, main = paste0(i))
  points(circ, rotation = "clock", stack = TRUE)
  
  hist(trunc(polar$theta),breaks = 70, col = "grey40", main = "41")
}












library(LHSpline)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/datasets")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")

fortyFive <- which(couponCov$polarAngle==45)
Zero <- which(couponCov$polarAngle==0)


splineSave <- matrix(NA, nrow = 54, ncol = length(fortyFive))
midsSave <- matrix(NA, nrow = 54, ncol = length(fortyFive))

i = 1

for(n in fortyFive){
load(paste0("nlsCoupon",n,".rda"))


# center the coupon
centerCoupon <- cbind(newCoupon[,1]-nlsCoeff["axisVectorX"],
                      newCoupon[,2]-nlsCoeff["axisVectorY"],
                      newCoupon[,3])




polar <- useful::cart2pol(centerCoupon[,1], centerCoupon[,2], degrees = TRUE)

periodicHist <- c(trunc(polar$theta)-360, trunc(polar$theta), trunc(polar$theta)+360)

h <- hist(periodicHist,
     breaks = 40, plot = FALSE)



LHFit <- splineDensity(periodicHist, N = length(h$breaks))


s <- h$mids
fHat1 <- dsplineDensity(LHFit, s)

hist(periodicHist,
     freq=FALSE, breaks = 80, col = "grey40", main = paste0(n))
lines(s, fHat1, col = "tomato", lwd = 2)
xline(c(0,360), col = "grey", lty = 2)

midsSave[,i] <- s
splineSave[,i] <- fHat1

i = i+1
}

load("nlsCoupon56.rda")
open3d()
par3d(cex=0.7)
plot3d(newCoupon[,1],
       newCoupon[,2],
       newCoupon[,3],
       type = "s", size = 0.45,
       xlab = " ", ylab = " ", zlab = " ")
polar <- useful::cart2pol(newCoupon[,1], newCoupon[,2])

plot(polar$theta, newCoupon[,3], 
     col = color.scale(polar$r, zlim = c(0,1170)),  pch = 16,
     main = "",
     xlab = "angle (radians)", ylab = "z axis",
     cex=1.3)

real <- which(midsSave[,1] >= 0 & midsSave[,1] <= 360)

nice <- which(fortyFive %in% c(55,53,34,5,1))

matplot(midsSave[,1][real], splineSave[real,nice], type = "l", lwd = 2)


## testing registration

J <- splineSave[real,nice]

dataGrid <- seq(0,365, length.out = 18)

matplot(dataGrid, J, type = "l", lwd = 2)
## want ind to reflect the max distance between bumps

ind <- rep(NA,length(nice))
for(i in 1:length(nice)){
  ind[i] <- findpeaks(J[,i])[,2][1]
}

ind <- c(findpeaks(J[,1])[,2][1], findpeaks(J[,2])[,2][1])

xline(midsSave[,1][ind]+360, col = "grey", lty = 3)

plot(dataGrid, J[,1], type = "l", lwd = 2, col = "violetred1", ylim = c(0.0004, 0.0015))
lines(dataGrid - diff(transPar), J[,2], lwd = 2, col = "violetred4")




transPar <- midsSave[,1][ind]+360



uGrid<-seq(-200,200,,18)
bigZ<-matrix(0,length(nice),length(uGrid))

matplot(dataGrid, J, type = "l", lwd = 2)



for( k in 1:length(nice)){
  # translated x values so maximum is at zero
  xTemp<- dataGrid-transPar[k]
  bigZ[k,]<-splint( xTemp,J[,k],uGrid)
  # set to missing any shifts outside data range
  indNA<- (uGrid< min( xTemp))|(uGrid> max( xTemp))
  bigZ[k,indNA]<- NA
  }


matplot(uGrid, t(bigZ), type = "l", lwd = 2)





##
## trying to register just the one large bumps now

#these we will have to wrap the curves. mayeb from -200 to 560

real <- which(midsSave[,1] >= -360 & midsSave[,1] <= 720)

nice <- which(fortyFive %in% c(58,49,42,40,28,26,25,11))

ctab = color.scale(nice)

matplot(midsSave[,1][real], splineSave[real,nice], type = "l", lwd = 2, col = ctab)
xline(c(0,365), col = "grey")

J <- splineSave[real,nice]
dataGrid <- midsSave[,1][real]

middle <- which(dataGrid >= 0 & dataGrid <= 360)

ind <- apply(J[middle,], 2, which.max)

xline(dataGrid[middle][ind], col = ctab)


transPar <- dataGrid[middle][ind]

uGrid<-seq(-200,200,,18)
bigZ<-matrix(0,length(nice),length(uGrid))




for( k in 1:length(nice)){
  # translated x values so maximum is at zero
  xTemp<- dataGrid-transPar[k]
  bigZ[k,]<-splint( xTemp,J[,k],uGrid)
  # set to missing any shifts outside data range
  indNA<- (uGrid< min( xTemp))|(uGrid> max( xTemp))
  bigZ[k,indNA]<- NA
}


matplot(uGrid, t(bigZ), type = "l", lwd = 2, col = ctab)
bumpShape <- colMeans( bigZ, na.rm=TRUE )

lines(uGrid, bumpShape, col = "darkgrey", lwd = 3)



three <- c(55,53,34,46,39,36,29,27,20,18,10,8,6,5,1)
one <- c(58,49,42,40,28,26,25,11)



