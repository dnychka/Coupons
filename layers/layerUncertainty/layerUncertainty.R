##
##
##
##
##
##
## layer testing with MC sampling
## -------------------------------------------------------------


library(astsa)
library(fields)
library(FSA)

# load required function files
setwd("C:/Users/barna/Documents/Coupons/layers/layerFunctions")
load("couponCov.rda")
source("calculateTestStatistic.R")
source("rotations.R")


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty")
source("functionsForBootstrapping.R")



## ------------------------------
## 0 degree coupons, real data
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

inputFiles <- list.files(full.names = TRUE)

zeros <- which(couponCov$polarAngle == 0)

N <- length(inputFiles[zeros]) #number of coupons in the sample size

nsamples = 100

TSZreal <- matrix(nrow = nsamples, ncol = N, NA) #test statistic (TS), zero degree, not layered

l = 1 #dummy index to increment coupon number

for(coupon in inputFiles[zeros]){
  
  load(coupon)
  
  k = 1 #dummy index to increment TS
  
  ## model output from first round of nls
  nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
  nlsZ <- nlsCoupon[,3]
  nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)
  
  
  ## begin the bootstrap for loop ------------
  
  
  bootNlsCoef <- matrix(ncol = 5, nrow = length(nsamples), NA)
  bootRadius <- matrix(ncol = length(nlsCoupon[,1]), nrow = length(nsamples), NA)
  
  
  for(i in 1:nsamples){
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
      M <- length(bootCoupon[,1])
      
      startValues <- getInitialParameters(bootCoupon)
      
      
      distTarget <- rep(0,M)
      
      nlsObj <- try(nls(distTarget~getDistance(bootCoupon, 
                                               centroidX, centroidY, 
                                               axisVectorX, axisVectorY, r),
                        start = list(centroidX = startValues[4], 
                                     centroidY = startValues[5],
                                     axisVectorX = startValues[1], 
                                     axisVectorY = startValues[2],
                                     r = startValues[7]),
                        control =  nls.control(minFactor = 1/10000000000)))
      
      
      
      
      if(class(nlsObj) != "try-error"){ # since nls needs babysitting
        
        
        if(j == 1){
          
          tempNlsCoeff <- coef(nlsObj)
          
          tempCoupon <- newCoupon(bootCoupon, tempNlsCoeff["centroidX"], tempNlsCoeff["centroidY"], 
                                  tempNlsCoeff["axisVectorX"], tempNlsCoeff["axisVectorY"])
          
          TSZreal[k,l] <- getTestStatistic(tempCoupon)
          k = k + 1
        }
        
        j <- j + 1
      }
      
      
      
    } #end of while loop
    
  } #end of bootstrap for loop
  
  # calculate the test statistic
  l = l + 1
  
}



## ------------------------------
## 0 degree coupons, layered
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/noLayers")

inputFiles <- list.files(full.names = TRUE)

N <- length(inputFiles) #number of coupons in the sample size

nsamples = 10

TSZlayers <- matrix(nrow = nsamples, ncol = 200, NA) #test statistic (TS), zero degree, not layered

l = 1 #dummy index to increment coupon number

for(coupon in inputFiles[1:200]){
  
  load(coupon)
  
  print(paste('*******', l, '*******'))
  
  k = 1 #dummy index to increment TS
  
  ## model output from first round of nls
  nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
  nlsZ <- nlsCoupon[,3]
  nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)
  
  
  ## begin the bootstrap for loop ------------
  
  
  bootNlsCoef <- matrix(ncol = 5, nrow = length(nsamples), NA)
  bootRadius <- matrix(ncol = length(nlsCoupon[,1]), nrow = length(nsamples), NA)
  
  
  for(i in 1:nsamples){
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
      M <- length(bootCoupon[,1])
      
      startValues <- getInitialParameters(bootCoupon)
      
      
      distTarget <- rep(0,M)
      
      nlsObj <- try(nls(distTarget~getDistance(bootCoupon, 
                                               centroidX, centroidY, 
                                               axisVectorX, axisVectorY, r),
                        start = list(centroidX = startValues[4], 
                                     centroidY = startValues[5],
                                     axisVectorX = startValues[1], 
                                     axisVectorY = startValues[2],
                                     r = startValues[7]),
                        control =  nls.control(minFactor = 1/10000000000)))
      
      
      
      
      if(class(nlsObj) != "try-error"){ # since nls needs babysitting
        
        
        if(j == 1){
          
          tempNlsCoeff <- coef(nlsObj)
          
          tempCoupon <- newCoupon(bootCoupon, tempNlsCoeff["centroidX"], tempNlsCoeff["centroidY"], 
                                  tempNlsCoeff["axisVectorX"], tempNlsCoeff["axisVectorY"])
          
          TSZlayers[k,l] <- getTestStatistic(tempCoupon)
          k = k + 1
        }
        
        j <- j + 1
      }
      
      
      
    } #end of while loop
    
  } #end of bootstrap for loop
  
  # calculate the test statistic
  l = l + 1
  
}


saveRDS(TSZlayers, "TSZnolayersUncertainty2.rds")


hist(TSZlayers)





TSZrealPtEst <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZreal.rds")
#TSZnoLayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZnolayers.rds")

TSZlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerUncertainty/TSZlayersUncertainty.rds")
TSZnolayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerUncertainty/TSZnolayersUncertainty.rds")
TSZnolayersTwo <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerUncertainty/TSZnolayersUncertainty2.rds")

# points(TSZrealPtEst, col = "violetred1", pch=20)
# 
# boxplot(as.vector(TSZreal))
# boxplot(TSZrealPtEst, add=TRUE, col = alpha("grey20", 0.3))
# 
# rhist<- hist(as.vector(TSZnoLayers), freq=F, xlim = c(0,0.8))
# lhist <- hist(TSZlayers, freq=F, add=T)
# xline(TSZrealPtEst, col = "mediumturquoise", lty=1, lwd=2)
# 
# plot(density(as.vector(TSZreal)), xlim = c(0,0.8))
# lines(density(TSZlayers))

length(as.vector(TSZreal))


TSZnoLayerDat <- cbind(c(as.vector(TSZnolayers), as.vector(TSZnolayersTwo)), rep(2, 4000))

TSZlayerDat <- cbind(as.vector(TSZlayers), rep(1,2950))

dat <- rbind(TSZlayerDat, TSZnoLayerDat)

dat <- data.frame(dat)
names(dat) <- c("signal", "population")

dat$type <- as.factor(dat$type)

library(ggplot2)



p <- ggplot(dat, aes(signal, fill=  population, col = population)) +
  ggtitle(label = "Distribution of Signal Strength Test Statistic, 0 degreee Coupons") +
  labs(x = "signal strength", y = "density") +
  geom_density(alpha=0.2, kernel="gaussian", adjust=1.5) +
  theme_bw()+ theme(legend.justification = c(1,1), 
                    legend.position = c(.95, .95), 
                    legend.background = element_blank(),
                    legend.box.background = element_rect(colour = "black"),
                    plot.title = element_text(size = 15),
                    axis.text=element_text(size=12),
                    axis.title.x = element_text(size=12),
                    axis.title.y = element_text(size=12),
                    legend.text=element_text(size=12),
                    legend.title = element_text(size=12)) +
  scale_color_manual(labels = c("layered", "not layered"), 
                     values = c("darkorange", 'cornflowerblue'), 
                     aesthetics = c("colour", "fill"),
                     name = "simulated populations") +
  geom_vline(data = data.frame(TSZrealPtEst),  aes(xintercept = TSZrealPtEst),
             linetype="dotted", color = "maroon3", alpha = 0.7, lwd=1) 
p


hist(TSZlayers, freq=F)
hist(TSZnoLayers, freq=F, add=F, breaks=30)

dat2 <- data.frame(TSZnoLayerDat)

p2 <- ggplot(dat2, aes(dat2$TSZnoLayers)) + stat_density()
p2
