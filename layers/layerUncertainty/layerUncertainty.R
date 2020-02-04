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
source("getKDEfunction.R")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("nlsFunctions.R")



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
  
  nsamples = 100
  
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
        
        tempNlsCoeff <- coef(nlsObj)
        
        tempCoupon <- newCoupon(bootCoupon, tempNlsCoeff["centroidX"], tempNlsCoeff["centroidY"], 
                                tempNlsCoeff["axisVectorX"], tempNlsCoeff["axisVectorY"])
        
        TSZreal[k,l] <- getTestStatistic(tempCoupon)
        k = k + 1
        
        break
        
      }
      
      j <- j+1
      
    } #end of while loop
    
  } #end of bootstrap for loop
  
  
  l = l + 1
  
}

saveRDS(TSZreal, "TSZrealUncertainty.rds")

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



## -----------------------------
## GRAPHICS
## ------------------------------


TSZrealPtEst <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZreal.rds")

TSZreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZrealUncertainty.rds")

TSZlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerUncertainty/TSZlayersUncertainty.rds")
TSZnolayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerUncertainty/TSZnolayersUncertainty.rds")
TSZnolayersTwo <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerUncertainty/TSZnolayersUncertainty2.rds")

TSZnoLayers <- cbind(TSZnolayers, TSZnolayersTwo)


dat <- rbind( cbind(as.vector(TSZlayers), rep(1,length(as.vector(TSZlayers)))),
              cbind(as.vector(TSZreal), rep(2,length(as.vector(TSZreal)))),
              cbind(as.vector(TSZnoLayers), rep(3, length(TSZnoLayers))) )


dat <- data.frame(dat)
names(dat) <- c("signal", "population")

dat$population <- as.factor(dat$population)

library(ggplot2)
library(extrafont)
font_import()
loadfonts(device = "win")

p <- ggplot(dat, aes(signal, fill=  population, col = population)) +
  ggtitle(label = "Distribution of Signal Strength Test Statistic, 45 degreee Coupons") +
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
                    legend.title = element_text(size=12),
                    text=element_text(size=16,  family="Palatino Linotype")) +
  scale_color_manual(labels = c("layered", "measured", "not layered"), 
                     values = c("darkorange", "violetred1", 'cornflowerblue'), 
                     aesthetics = c("colour", "fill"),
                     name = "coupon populations")  +
  annotate("point", x = TSZrealPtEst, y = rep(-.17, 16), pch = "|", col = "maroon3", cex = 2.2)
p

ggsave("font_ggplot.pdf", plot=p,  width=4, height=4)
## ------------------------------
## 45 degree coupons, real data
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

inputFiles <- list.files(full.names = TRUE)

ff <- which(couponCov$polarAngle == 45)

N <- length(inputFiles[ff]) #number of coupons in the sample size

nsamples = 100

angleSeq <- seq(0, 2*pi, length.out = 100) #increment by which we change rotation

tempTS <- rep(NA, 100) #store TS for each rotation

TSFreal <- matrix(nrow = nsamples, ncol = N, NA) #test statistic (TS), zero degree, not layered

l = 1 #dummy index to increment coupon number

for(coupon in inputFiles[ff]){
  
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
      
      kdeSample <- getKDEzt(nlsCoupon, 0.07, 900, 300, length(nlsZ))
      
      bootTheta <- kdeSample[[1]]
      bootZ <- kdeSample[[2]]
      
      
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
          
          ## where we put the angle rotations
        
          m = 1
          for(w in angleSeq){
            
            rotatedCoupon <- rotateFortyFive(w, tempCoupon) #center and rotate the coupon
            
            tempTS[m] <- getTestStatistic(rotatedCoupon) #calc TS for rotated coupon
            
            m=m+1
          }
          
          TSFreal[k,l] <- max(tempTS)
          
      
          k = k + 1
        }
        
        j <- j + 1
      }
      
      
      
    } #end of while loop
    
  } #end of bootstrap for loop
  
  # calculate the test statistic
  l = l + 1
  
}

saveRDS(TSFreal, "TSFreal.rds")

## 


## ------------------------------
## 45 degree coupons, LAYERED
## ------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/spacing40to60")

inputFiles <- list.files(full.names = TRUE)

N <- length(inputFiles) #number of coupons in the sample size

nsamples = 10

angleSeq <- seq(0, 2*pi, length.out = 100) #increment by which we change rotation

tempTS <- rep(NA, 100) #store TS for each rotation

TSFlayers <- matrix(nrow = nsamples, ncol = N, NA) #test statistic (TS), zero degree, not layered


i = 1 #dummy index to increment TS

for(coupon in inputFiles){
  
  load(coupon)
  
  print(i)
  
  m = 1
  for(w in angleSeq){
    
    rotatedCoupon <- rotateFortyFive(w, nlsCoupon) #center and rotate the coupon
    
    tempTS[m] <- getTestStatistic(rotatedCoupon) #calc TS for rotated coupon
    
    m=m+1
  }
  
  TSFlayers[i] <- max(tempTS)
  i = i + 1
}

saveRDS(TSFlayers, "TSFlayers.rds")

## 



## -----------------------------
## GRAPHICS 45 degree
## -----------------------------

TSFrealPtEst <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFreal.rds")
TSFreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerUncertainty/TSFreal.rds")

TSFlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFlayers.rds")
TSFlayers <- TSFlayers[ , colSums(is.na(TSFlayers)) == 0]

TSFnoLayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFnoLayers.rds")



dat <- rbind( cbind(as.vector(TSFlayers), rep(1,length(as.vector(TSFlayers)))),
              cbind(as.vector(TSFreal), rep(2,length(as.vector(TSFreal)))),
              cbind(as.vector(TSFnoLayers), rep(3, length(TSFnoLayers))) )
            

dat <- data.frame(dat)
names(dat) <- c("signal", "population")

dat$population <- as.factor(dat$population)

library(ggplot2)
library(extrafont)
font_import()
loadfonts(device = "win")

p <- ggplot(dat, aes(signal, fill=  population, col = population)) +
  ggtitle(label = "Distribution of Signal Strength Test Statistic, 45 degreee Coupons") +
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
                    legend.title = element_text(size=12),
                    text=element_text(size=16,  family="Palatino Linotype")) +
  scale_color_manual(labels = c("layered", "measured", "not layered"), 
                     values = c("darkorange", "violetred1", 'cornflowerblue'), 
                     aesthetics = c("colour", "fill"),
                     name = "coupon populations")  +
  annotate("point", x = TSFrealPtEst, y = rep(-.17, 24), pch = "|", col = "maroon3", cex = 2.2)
p

ggsave("font_ggplot.pdf", plot=p,  width=4, height=4)

