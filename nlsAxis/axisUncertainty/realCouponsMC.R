##
##
##
##
##
## Monte Carlo, sampling theta and z from bivariate KDE
## using real coupon data
### -------------------------------------------------------------------

setwd("Z:/COUPONS/Coupons-master/datasets")
nPores <- readRDS("0and45degreePoreNumbers.rds")
position <- readRDS("buildPlatePosition.rds")

setwd("Z:/COUPONS/Coupons-master/nlsAxis/axisUncertainty")
source("getKDEfunction.R")


setwd("Z:/COUPONS/Coupons-master/nlsAxis")
source("nlsFunctions.R") #lolol don't set wd in functions


##
## get the data
## --------------------------------------------------------------------

setwd("Z:/COUPONS/Coupons-master/nlsAxis/porosity/nlsPorosityData/cropped/named")

couponList <- list.files(full.names = TRUE)

names(couponList) <- position[order(position)]

couponList <- couponList[which(names(couponList) %in% names(nPores))]


poreInd <- 1
for(n in couponList){ # begin the coupon for loop
  
  setwd("Z:/COUPONS/Coupons-master/nlsAxis/porosity/nlsPorosityData/cropped/named")
  
  load(n)
  
  ## model output from first round of nls
  nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
  nlsZ <- nlsCoupon[,3]
  nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)
  
  
  ## begin the bootstrap for loop ------------
  
  nsamples = 1000
  
  simNlsCoef <- matrix(ncol = 5, nrow = nsamples, NA)
  simRadius <- matrix(ncol = length(nlsCoupon[,1]), nrow = nsamples, NA)
  
  
  for(i in 1:nsamples){
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
        
        simNlsCoef[i,] <- coef(nlsObj)
        
        tempNlsCoeff <- coef(nlsObj)
        
        tempCoupon <- newCoupon(bootCoupon, tempNlsCoeff["centroidX"], tempNlsCoeff["centroidY"], 
           tempNlsCoeff["axisVectorX"], tempNlsCoeff["axisVectorY"])
        tempRadius <- sqrt(                     tempCoupon[,1]^2+tempCoupon[,2]^2)
        
        simRadius[i,] <- tempRadius
        
        break
        
      }
      
      j <- j+1
      
    } #end of while loop
    
  } #end of bootstrap for loop



  
 
  setwd("Z:/COUPONS/newporeuncertaintydata")
  save(simNlsCoef, simRadius, file = paste0("coupon", names(nPores[poreInd]), "poreCI.rda"))

  poreInd <- poreInd + 1

} #end of coupon reading-in for loop












