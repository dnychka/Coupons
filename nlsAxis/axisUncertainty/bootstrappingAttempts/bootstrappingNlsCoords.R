##
##
##
##
##
## bootstrap off the xyz positions
## instead of the residuals
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



## add a pore ID column to nlsCoupon
nlsCouponWithID <- cbind(nlsCoupon, seq_along(nlsCoupon[,1]))



##
## data based bootstrapping for nlsObj
## -----------------------------------------



## begin the bootstrap for loop ------------

nsamples = 1:1000

bootNlsCoefXYZ <- matrix(ncol = 5, nrow = length(nsamples), NA)
bootRadiusXYZ <- matrix(ncol = length(nlsCoupon[,1]), nrow = length(nsamples), NA)

bootRadiusByPoreID <- list()
for(i in length(poreCoordinates[,1])){
  bootRadiusByPoreID[[i]] <- vector()
}
names(bootRadiusByPoreID) <- nlsCouponWithID[,4]


for(i in nsamples){
  print(i)
  
  j = 1
  while(j < 10){ #try nls up to 10 times
    
    
    ## the bootstrap step
    bootCoupon <- nlsCouponWithID[sample(nrow(nlsCouponWithID),
                                         size=length(nlsCouponWithID[,1]),replace=TRUE),]   
    
    
    ## save ID and remove it from coordinate data
    bootID <- bootCoupon[,4]
    bootCoupon <- bootCoupon[,-4]
    
    
    # estimate the initial parameter values for each new fabricated data set
    N <- length(bootCoupon[,1])
  
    startValues <- getInitialParameters(bootCoupon)
    
    
    
    ## run nls
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
      
      bootNlsCoefXYZ[i,] <- coef(nlsObj)
      
      tempNlsCoeff <- coef(nlsObj)
      
      tempCoupon <- newCoupon(bootCoupon, tempNlsCoeff["centroidX"], tempNlsCoeff["centroidY"], 
                              tempNlsCoeff["axisVectorX"], tempNlsCoeff["axisVectorY"])
      tempRadius <- sqrt(tempCoupon[,1]^2+tempCoupon[,2]^2)
      
      bootRadiusXYZ[i,] <- tempRadius
      
      
      if(j == 1){
      
          k=1
          
          for (m in bootID){
            
            bootRadiusByPoreID[[m]] <- c(bootRadiusByPoreID[[m]],bootRadiusXYZ[i,k])
            
            k=k+1
            
          }
      }
      
      
      j <- j + 1
    }
    
    
    
  } #end of while loop
  
} #end of for loop


#saveRDS(bootRadiusByPoreID, "bootRadiusByPoreID.rds")

bootRadiusByPoreID <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/uncertainty/bootRadiusByPoreID.rds")


poreCIs <- matrix(nrow = 21, ncol = length(nlsCoupon[,1]), NA)

for(i in 1:length(nlsCoupon[,1])){
  poreCIs[,i] <- quantile(bootRadiusByPoreID[[i]], prob = seq(0, 1, by = 0.05))
}

max(poreCIs[2,])

##
## plot it
## ----------------------------------


nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
nlsZ <- nlsCoupon[,3]
nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)

trueR <- sqrt(trueCoupon[,1]^2+trueCoupon[,2]^2)

draw <- sample(1:4025,1500, replace=F)

plot(nlsTheta[draw], nlsR[draw], 
     col = "white", ylim = c(810,850), 
     main = expression(paste("95% confidence intervals, bootstrapping off of (r,", theta, ", z)")),
     xlab = "theta",
     ylab = "radius")

arrows(x0=nlsTheta[draw], y0=poreCIs[2,draw], 
       x1=nlsTheta[draw], y1=poreCIs[20,draw], 
       code=3, angle=90, length=0.1, lty = 2, col = "steelblue", lwd=1.5)

points(nlsTheta[draw], nlsR[draw], pch=16, col="darkorange", cex=1.6)
points(nlsTheta[draw], trueR[draw], pch = 20, col= "black")


legend("topright", c("nls coupon", "true coupon"), pch = c(16,20), col = c("darkorange", "black"))


## ? how often do we caputre the true value?

CIdiff <- abs(poreCIs[2,]-poreCIs[20,])

