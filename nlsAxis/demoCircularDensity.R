




library(rgl)
library(circular)
library(fields)


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/couponCaseStudies/caseStudyData")

load("nlsCoupon55.rda")

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
  
  hist(trunc(polar$theta),breaks = 70, col = "grey40", main = paste0(i))
}




kde<- density(circ, bw = 100,control.circular=list(units="degrees"))

lines(kde, col = "grey70", lwd = 2)

rose.diag(circ,bins=24,shrink=0.23,xlim=c(-2,2),ylim=c(-2,2),
          axes=FALSE,prop=1.5)


