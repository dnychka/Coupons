##
##
##
##
##
## Monte Carlo, sampling theta and z from bivariate KDE
## using real coupon data
### -------------------------------------------------------------------

setwd("C:/Users/barna/Documents/Coupons/datasets")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty")
source("getKDEfunction.R")
source("getAccelerationFunction.R")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("nlsFunctions.R") #lolol don't set wd in functions


##
## get the data
## --------------------------------------------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

couponList <- list.files(full.names = TRUE)

niceCoupon <- c(which(couponCov$polarAngle==0), which(couponCov$polarAngle==45))

couponList <- couponList[niceCoupon]


for(n in 4){ # begin the coupon for loop
  
  load(paste0("./nlsCoupon", n, ".rda"))
  
  # nlsCoupon already fit by nls (old nls...iterative median method)
  
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
      # thetas and z values from the bivariate KDE of (z, theta)
      #bootTheta <- runif( length(nlsCoupon[,1]),0,2*pi)
      
      kdeSample <- getKDEzt(nlsCoupon, 0.07, 900, 300, length(nlsZ))
      
      bootTheta <- kdeSample[[1]]
      bootZ <- kdeSample[[2]]
      
      
      # convert new dataset to cylindrical coordinates
      bootX <- nlsR*cos(bootTheta)
      bootY <- nlsR*sin(bootTheta)
      
      bootCoupon <- cbind(bootX, bootY, bootZ)  
      
      
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
        
        bootNlsCoef[i,] <- coef(nlsObj)
        
        tempNlsCoeff <- coef(nlsObj)
        
        tempCoupon <- newCoupon(bootCoupon, tempNlsCoeff["centroidX"], tempNlsCoeff["centroidY"], 
                                tempNlsCoeff["axisVectorX"], tempNlsCoeff["axisVectorY"])
        tempRadius <- sqrt(tempCoupon[,1]^2+tempCoupon[,2]^2)
        
        bootRadius[i,] <- tempRadius
        
        break
        
      }
      
      j <- j+1
      
    } #end of while loop
    
  } #end of bootstrap for loop





  ## ---------------------------------------------------------------------
  ## compute bias-corrected and accelerated CI
  ## ---------------------------------------------------------------------
  setwd("C:/Users/barna/Documents/Coupons/datasets")
  poreData <- readRDS("porosityData.rds")
  
  poreCoordinates <- cropCoupon(n, poreData)
  
  aHat <- getAcceleration(poreCoordinates)
  zHat <- vector()
  for(i in 1:5){
    zHat[i] <- qnorm(sum(ifelse(bootNlsCoef[,i] < nlsCoeff[i], 1, 0))/length(nsamples))
  }

  
  
  paramCI <- matrix(nrow = 2, ncol = length(nlsCoeff), NA)
  
  alpha = 0.025 #significance level
  
  lo <- vector(); high <- vector()
  for(i in 1:length(nlsCoeff)){
    lo <- pnorm(zHat[i] + (zHat[i] + qnorm(alpha)) / (1 - aHat[i] * (zHat[i] + qnorm(alpha)) ))
    high <- pnorm(zHat[i] + (zHat[i] + qnorm(1-alpha)) / (1 - aHat[i] * (zHat[i] + qnorm(1-alpha)) ))
    
    paramCI[,i] <- quantile(bootNlsCoef[,i], c(lo, high))
  }
  
 

} #end of coupon reading-in for loop




par(mfrow = c(2,2))
for(i in 1:4){ hist(bootNlsCoef[,i]); xline(paramCI[,i], col = "violetred1", lwd=2) }








