##
##
##
##
##
##
## fit the surface using re-parameterized nls (5 parameters)
## --------------------------------------------------------

library(rgl)

source("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/functionsForBootstrapping.R")


# start wih coupon 4 (E22)

surfaceCoordF3 <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/surfaces/nlsSurfaceData/F3/surfaceCoordsF3.rds")

# downsample (every 20th from every 5th layer orig read in)
z <- surfaceCoordF3[,3]

zSmall <- seq(min(z), max(z), by = 20)

smallSurface <- surfaceCoordF3[surfaceCoordF3[,3] %in% zSmall,]

plot3d(smallSurface[,1], smallSurface[,2], smallSurface[,3], type="s", size=0.45)


##
## fit it
## --------------------------

surfaceNlsObjE22 <- nlsAxisFit(smallSurface)

nlsCoeff <- coef(surfaceNlsObjE22)

nlsCoupon <- newCoupon(smallSurface, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])


## --------------------------

##
## theta-based MC for surface
## ---------------------------

## model output from first round of nls
nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
nlsZ <- nlsCoupon[,3]
nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)

## begin the bootstrap for loop ------------

nsamples = 1:1000

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



## create pivot for nls coeff confidence intervals

pivMatrixCoef <- matrix(ncol = length(nlsCoeff), nrow = length(nsamples), NA)

bootCoefCI <- matrix(nrow = 2, ncol = length(nlsCoeff), NA)

for (i in 1:length(nlsCoeff)){
  
  pivMatrixCoef[,i] <- (bootNlsCoef[,i]-nlsCoeff[i]) / sd(bootNlsCoef[,i])
  
  bootCoefCI[,i] <- mean(bootNlsCoef[,i]) - sd(bootNlsCoef[,i])*quantile(pivMatrixCoef[,i],c(.025,.975))
  
}

name <- "F3"

save(bootNlsCoef, bootCoefCI, bootRadius, nlsCoeff, nlsCoupon, file= paste0(name, "surfaceData.rda"))
