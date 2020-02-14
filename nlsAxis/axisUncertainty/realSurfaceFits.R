##
##
##
##
##
## Monte Carlo, sampling theta and z from bivariate KDE
## using surface coupon data
### -------------------------------------------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty")
source("getKDEfunction.R")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("nlsFunctions.R") #lolol don't set wd in functions


##
## get the data
## --------------------------------------------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/surfaces/nlsSurfaceData")

couponList <- list.files(full.names = TRUE)

for(n in couponList){ # begin the coupon for loop
  
  setwd("C:/Users/barna/Documents/Coupons/nlsAxis/surfaces/nlsSurfaceData")
  
  load(n)
  
  print(n)
  
  ## model output from first round of nls
  nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
  nlsZ <- nlsCoupon[,3]
  nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)
  
  
  ## begin the bootstrap for loop ------------
  
  nsamples = 1000
  
  nSurfacePxls = 2000
  
  simNlsCoef <- matrix(ncol = 5, nrow = nsamples, NA)
  simRadius <- matrix(ncol = nSurfacePxls, nrow = nsamples, NA)
  
  
  for(i in 1:nsamples){
    print(i)
    
    j = 1
    while(j < 25){ #try nls up to 25 times
      
      bootCoupon <- oldCoupon[sample(oldCoupon, 2000, replace = F),]
      
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
        tempRadius <- sqrt(tempCoupon[,1]^2+tempCoupon[,2]^2)
        
        simRadius[i,] <- tempRadius
        
        break
        
      }
      
      j <- j+1
      
    } #end of while loop
    
  } #end of bootstrap for loop
  
  cNum <- gsub(".rda", "", gsub("./nlsSurfaceCoupon", "", n))
  
  
  
  setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData/surfaces")
  save(simNlsCoef, simRadius, nlsCoupon, file = paste0("coupon", cNum, "surfaceCI.rda"))
  
} #end of coupon reading-in for loop












