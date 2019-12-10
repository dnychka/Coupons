##
##
##
##
##
## bootstrap (really MC) by sampling theta from a uniform (0,2pi)
### -------------------------------------------------------------------


myDir <- "C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/uncertaintyData"

setwd(myDir)

# load the sample simulated coupon
trueCoupon <- readRDS(file.path(myDir, "/simulatedCouponInternalJitter.rds")) #has NOT been nls-ed

# knock it off its axis
aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}
theta = -13*pi/180

Ry <- aboutY(theta)

poreCoordinates <- trueCoupon %*% Ry 

rm(myDir, aboutY, theta, Ry)


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty")
source("functionsForBootstrapping.R")




##
## fit it
## -----------------------------------------


nlsObj <- nlsAxisFit(poreCoordinates)

nlsCoeff <- coef(nlsObj)

nlsCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])



##
## residual based bootstrapping for nlsObj
## -----------------------------------------


## model output from first round of nls
nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
nlsZ <- nlsCoupon[,3]
nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)


## begin the bootstrap for loop ------------

nsamples = 1:10

bootNlsCoef <- matrix(ncol = 5, nrow = length(nsamples), NA)
bootRadius <- matrix(ncol = length(nlsCoupon[,1]), nrow = length(nsamples), NA)


for(i in nsamples){
  print(i)
  
  j = 1
  while(j < 10){ #try nls up to 10 times
    
    # generate stochastic data set by sampling possible
    # thetas from uniform(0, 2pi)
    bootTheta <- runif( length(nlsCoupon[,1]),0,2*pi)
    
    
    
    # fabricate new data by adding the sampled errors to the 
    # model output
    bootX <- nlsR*cos(bootTheta)
    bootY <- nlsR*sin(bootTheta)
    
    bootCoupon <- cbind(bootX, bootY, nlsZ)    
    
    # orient the bootCoupon so it matches the tilt in the 
    # original data; that way, nlsCoeff will match taking 
    # coupons from this original tilt -> alignment
    bootCoupon <- getOldCoupon(bootCoupon, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                               nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
    
    
    
    # estimate the parameter values for each new fabricated data set
    N <- length(bootCoupon[,1])
    
    startValues <- getInitialParameters(bootCoupon)
    
    
    distTarget <- rep(0,N)
    
    nlsObj <- try(nls(distTarget~getDistance(bootCoupon, 
                                             centroidX, centroidY, 
                                             axisVectorX, axisVectorY, r),
                      start = list(centroidX = startValues[4], 
                                   centroidY = startValues[5],
                                   axisVectorX = startValues[1], 
                                   axisVectorY = startValues[2],
                                   r = startValues[7]),
                      control =  nls.control(minFactor = 1/10000000000)))
    
    
    
    
    if(class(nlsObj) != "try-error"){
      j <- j + 1
      bootNlsCoef[i,] <- coef(nlsObj)
      
      tempNlsCoeff <- coef(nlsObj)
      
      tempCoupon <- newCoupon(bootCoupon, tempNlsCoeff["centroidX"], tempNlsCoeff["centroidY"], 
                              tempNlsCoeff["axisVectorX"], tempNlsCoeff["axisVectorY"])
      tempRadius <- sqrt(tempCoupon[,1]^2+tempCoupon[,2]^2)
      
      bootRadius[i,] <- tempRadius
      
      
    }
    
    
    
  } #end of while loop
  
} #end of for loop



#saveRDS(bootRadius, "bootRadiusTheta.rds")
#saveRDS(bootNlsCoef, "bootNlsCoefTheta.rds")

bootNlsCoef <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/uncertaintyData/bootNlsCoefTheta.Rds")


boot95coef <- apply(bootNlsCoef, 2, function(x) quantile(x, prob = seq(0, 1, by = 0.05)))


bootRadius <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/uncertaintyData/bootRadiusTheta.Rds")

nsamples=1000

## create pivot for pore 1
pivMatrix <- matrix(ncol = ncol(bootRadius), nrow = nsamples, NA)

for (i in 1:ncol(bootRadius)){

pivMatrix[,i] <- (bootRadius[,i]-nlsR[i]) / sd(bootRadius[,i])


}

hist(pivMatrix[,1])


pivotCIpore <- matrix(nrow = 2, ncol = ncol(bootRadius), NA)

for (i in 1:ncol(pivMatrix)){

pivotCIpore[,i] <- mean(bootRadius[,i]) - sd(bootRadius[,i])*quantile(pivMatrix[,i],c(.025,.975))

}


boot95radius <- apply(bootRadius, 2, function(x) quantile(x, c(0.025,.975)))


boot95radius[,1:10] #pretty much the same
pivotCIpore[,1:10]  #pretty much the same


## how often do we actually capture the true position?
isIn <- nlsR

trueR <- sqrt(trueCoupon[,1]^2+trueCoupon[,2]^2)

for( i in 1:length(nlsR)){
 isIn[i] <- (trueR[i] >= boot95radius[1,i]  & trueR[i] <= boot95radius[2,i])
}

table(isIn)
1-154/3871


isInPiv <- nlsR
for( i in 1:length(nlsR)){
  isInPiv[i] <- (trueR[i] >= pivotCIpore[2,i]  & trueR[i] <= pivotCIpore[1,i])
}

table(isInPiv)
1-142/3883


##
## plotting stuff
## ----------------------------------------------------

draw <- sample(1:4025,5, replace=F)

plot(nlsTheta[draw], boot95radius[2,][draw], ylim = c(790,805),
     col = "white", 
     main = "15 example pores and their 95% confidence intervals",
     xlab = "theta",
     ylab = "radius")
points(boot95radius[20,draw], col = "white")
#yline(c(mean(boot95radius[2,]), mean(boot95radius[20,])), lwd=4, lty = 2, col="steelblue")

arrows(x0=nlsTheta[draw], y0=boot95radius[1,draw], 
       x1=nlsTheta[draw], y1=boot95radius[2,draw], 
       code=3, angle=90, length=0.1, lty = 2, col = "steelblue", lwd=2)

arrows(x0=nlsTheta[draw], y0=pivotCIpore[1,draw], 
       x1=nlsTheta[draw], y1=pivotCIpore[2,draw], 
       code=3, angle=90, length=0.1, lty = 2, col = "tomato")


plot(nlsTheta[draw], sqrt(nlsCoupon[,1]^2 + nlsCoupon[,2]^2)[draw], pch=16, col="darkorange", cex=1)
points(nlsTheta[draw], sqrt(trueCoupon[,1]^2 + trueCoupon[,2]^2)[draw], pch = 20, col= "black")


legend("topright", c("nls coupon", "true coupon"), pch = c(16,20), col = c("darkorange", "black"))







