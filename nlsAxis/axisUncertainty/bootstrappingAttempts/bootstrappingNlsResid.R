



getDistance <- function(poreCoordinates, centroidX, centroidY, axisVectorX, axisVectorY, r){
  
  centroidGuess <- c(centroidX,centroidY,0)
  
  axisGuess <- c(axisVectorX,axisVectorY,1)
  
  newCoupon <- alignCoupon(poreCoordinates, centroidGuess, axisGuess)
  
  newRadius <- sqrt(newCoupon[,1]^2+newCoupon[,2]^2)
  
  dist <- (newRadius-rep(r,length(newRadius)))^2
  
  return(dist)
}


getInitialParameters <- function(poreCoordinates){
  
  ## get centers of mass for upper and lower half of the coupon to 
  ## compute axisVector, which is the direction vector for 
  ## the initial axis estimate
  deciles <- quantile(poreCoordinates[,3], prob = seq(0, 1, length = 11), type = 5)
  
  upperHalf <- poreCoordinates[,3] >= deciles[6]
  
  lowerHalf <- poreCoordinates[,3] <= deciles[6]
  
  axisVector <- colMeans(poreCoordinates[upperHalf,]) - colMeans(poreCoordinates[lowerHalf,])
  
  axisVector <- axisVector / sqrt(sum(axisVector^2)) # make it unit length
  
  
  ## project the centroid into the x-y plane along axisVector.
  ## allows us the parameterize the centroid using two parameters
  ## (x,y,0) rather than three (x,y,z)
  centroid <- colMeans(poreCoordinates)
  
  xyCentroid <- c( (axisVector[1]*-centroid[3])/axisVector[3] + centroid[1],
                   (axisVector[2]*-centroid[3])/axisVector[3] + centroid[2],
                   0)
  
  r = 1000 #hard coded for testing
  
  startingValues <- c(as.vector(axisVector), as.vector(xyCentroid), r)
  
  return(startingValues)
  
}

alignCoupon <- function(poreCoordinates, xyCentroid, axisVector){
  
  V1 <- axisVector*sign(axisVector[3])
  
  U1<- V1/ sqrt(sum( V1^2))
  
  I3<- diag( 1,3)
  
  V<- qr.qy(qr( cbind( U1, c(-1,0,0), c( 0,-1,0) )),I3) 
  
  V<- V[,c(2,3,1)]
  
  newXYZ <- t(V) %*% (t(poreCoordinates)-xyCentroid)
  
  newXYZ <- t(newXYZ)
  
  return(newXYZ)
  
}





nlsAxisFit <- function(poreCoordinates){  

  startValues <- getInitialParameters(poreCoordinates)
  
  N <- length(poreCoordinates[,1])

  
  distTarget <- rep(0,N)
    
  centerAxis <- nls(distTarget~getDistance(poreCoordinates, 
                                                 centroidX, centroidY, 
                                                 axisVectorX, axisVectorY, r),
                      start = list(centroidX = startValues[4], 
                                   centroidY = startValues[5],
                                   axisVectorX = startValues[1], 
                                   axisVectorY = startValues[2],
                                   r = startValues[7]),
                      control =  nls.control(minFactor = 1/10000000000))
    
    
  nlsCoeff <- coef(centerAxis)
  nlsCoeff
  
  
  return(centerAxis)
    
}



setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("newCoupon.R")

###
##
### -------------------------------------------------------------------

myDir <- "C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty"

# load the sample simulated coupon
trueCoupon <- readRDS(file.path(myDir, "/simulatedCouponInternalJitter.rds")) #has NOT been nls-ed

# knock it off its axis
aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}
theta = -8*pi/180

Ry <- aboutY(theta)

poreCoordinates <- trueCoupon %*% Ry 

rm(myDir, aboutY, theta, Ry)

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


## get resiudals from model and change them to cartesian coords
nlsResid <- residuals(nlsObj) * -1

nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
nlsZ <- nlsCoupon[,3]

#residualToDist <- sqrt(nlsResid) + nlsCoeff[5]

residualToDist <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)


## begin the bootstrap for loop ------------

nsamples = 1:10

bootNlsCoef <- matrix(ncol = 5, nrow = length(nsamples), NA)
bootRadius <- matrix(ncol = length(nlsCoupon[,1]), nrow = length(nsamples), NA)


for(i in nsamples){
  print(i)
  
  j = 1
  while(j < 10){ #try nls up to 10 times
    
    # resample from nlsRes using non-parametric sample to get
    # new error values for every pore location
    bootR <- sample(residualToDist,size=length(residualToDist),replace=TRUE)
    
    
    
    # fabricate new data by adding the sampled errors to the 
    # model output
    bootX <- bootR*cos(nlsTheta)
    bootY <- bootR*sin(nlsTheta)

    bootCoupon <- cbind(bootX, bootY, nlsZ)    
  
    
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



# saveRDS(bootRadius, "bootRadius.rds")
# saveRDS(bootNlsCoef, "bootNlsCoef.rds")

bootNlsCoef <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/bootNlsCoef.Rds")


boot95coef <- apply(bootNlsCoef, 2, function(x) quantile(x, prob = seq(0, 1, by = 0.05)))



 bootRadius <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/bootRadius.Rds")




boot95radius <- apply(bootRadius, 2, function(x) quantile(x, prob = seq(0, 1, by = 0.05)))
boot95radius[c(2,20),1:10]

draw <- sample(1:4015,15, replace=F)

plot(nlsTheta[draw], boot95radius[2,][draw], 
     col = "white", ylim = c(740,890), 
     main = "15 example pores and their 95% confidence intervals",
     xlab = "theta",
     ylab = "radius")
points(boot95radius[20,draw], col = "white")
#yline(c(mean(boot95radius[2,]), mean(boot95radius[20,])), lwd=4, lty = 2, col="steelblue")

arrows(x0=nlsTheta[draw], y0=boot95radius[2,draw], 
       x1=nlsTheta[draw], y1=boot95radius[20,draw], 
       code=3, angle=90, length=0.1, lty = 2, col = "steelblue")


points(nlsTheta[draw], sqrt(nlsCoupon[,1]^2 + nlsCoupon[,2]^2)[draw], pch=16, col="darkorange", cex=1.6)
points(nlsTheta[draw], sqrt(trueCoupon[,1]^2 + trueCoupon[,2]^2)[draw], pch = 20, col= "black")


legend("topright", c("nls coupon", "true coupon"), pch = c(16,20), col = c("darkorange", "black"))



