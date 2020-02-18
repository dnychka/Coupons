## 07/08/2019
##
##
##
##
##
## read in selected tif files and convert the surface
## data to coords
##
## -------------------------------------------------------------------------


library(rgl)
library(raster)
library(rgdal)
library(gtools)


## -----------------------
## which coupons do we want
## -----------------------
p <- readRDS("C:/Users/barna/Documents/Coupons/datasets/preparationData.rds")
location <- rep(NA, 58)
for(i in 1:58){
  location[i] <- paste(unlist(p[[i]]$scalars[3:4]), collapse = "")
}

location

location[order(location)]

load("C:/Users/barna/Documents/Coupons/datasets/couponCov.rda")
location[which(couponCov$polarAngle==0)]
location[which(couponCov$polarAngle==45)]
location[which(couponCov$polarAngle==90)][order(location[which(couponCov$polarAngle==90)])]

location[17]

## -----------------------------------
## functions
## -----------------------------------
## extract only pixels with ones values
getSurface <- function(ind) {
  allPts <- rasterToPoints(rasterStack[[ind]], function(x){x==255})[,-3]
  n <- nrow(allPts)
  surfacePts <- allPts[sample(n, 0.05*n, replace=FALSE),]
  return(surfacePts)}


## ----------------------------------
## Construct and save surface data
## ----------------------------------

setwd("D:/server-files/SURFACES")

folderList <- list.files()

for(folder in folderList){
  
  setwd(paste0("D:/server-files/SURFACES/",folder))
  
  print(folder)
  
  startTime <- Sys.time()
  
  rasterList <- mixedsort(list.files(full.names = TRUE))
  
  numTif <- length(rasterList)
  
  rasterList <- rasterList[c(-(1:100),-(numTif:(numTif-100)))] #take off XCT artifacts at top and bottom
  
  choice <- seq(1,length(rasterList),by=10) #subset the slices of the coupon
  
  rasterStack <- stack(rasterList[choice])
  
  
  
  X <- vector()
  Y <- vector()
  Z <- vector()
  
  
  for(i in 1:dim(rasterStack)[3]){
    
    slice <- getSurface(i)
    
    X <- c(X,slice[,1])
    Y <- c(Y,slice[,2])
    Z <- c(Z, rep(choice[i], length(slice[,1])))
  }
  endTime <- Sys.time()
  
  print(endTime-startTime)
  
  resolution <- 4.5 # multiply by XCT resolution (roughly)
  X <- X*resolution; Y <- Y*resolution; Z <- Z*resolution
  surfaceCoupon <- cbind(X,Y,Z)
  
  
  saveRDS(surfaceCoupon, paste0(folder,".rds"))
  
}

