##
##
##
##
##
## read in .npy files with connected components indices
## ----------------------------------------------------------


library(reticulate)
library(fields)
library(sp)
library(car)
library(raster)
library(rgdal)
library(gtools)
library(scatterplot3d)


## ----------------------------------------------------------
## read in volume python files
## ----------------------------------------------------------

setwd("D:/recon/processed")


np <- import("numpy")


# data reading
matConda1 <- np$load("testSmall1.npy")

xyzg <- matrix(NA, nrow = 1, ncol = 4)

numTif <- dim(matConda1)[1]

for(ind in 1:numTif){

  babyConda <- raster(matConda1[ind,,])
  
  
  whitePts <- rasterToPoints(babyConda, function(x){x!=0})
  
  if(dim(whitePts) != 0){
    
    print(ind)
    
    ridX <- (names(which(table(whitePts[,1]) > 100)))
    porePts <- whitePts[-which(whitePts[,1] %in% ridX),]

  
    connectedPix <- extract(babyConda, SpatialPoints(porePts[,1:2]), sp = T)
  
    numPts <- length(connectedPix@data$layer)
  
    xyzg <- rbind(xyzg, cbind(connectedPix@coords, rep(ind, numPts), connectedPix@data$layer))
  }
  
}

rm(matConda1)
gc()

matConda2 <- np$load("testSmall2.npy")


for(ind in 1:dim(matConda2)[1]){
  
  babyConda <- raster(matConda2[ind,,])
  
  
  whitePts <- rasterToPoints(babyConda, function(x){x!=0})
  
  if(dim(whitePts) != 0){
    
    print(ind+numTif)
    
    ridX <- (names(which(table(whitePts[,1]) > 100)))
    porePts <- whitePts[-which(whitePts[,1] %in% ridX),]
    
    
    connectedPix <- extract(babyConda, SpatialPoints(porePts[,1:2]), sp = T)
    
    numPts <- length(connectedPix@data$layer)
    
    xyzg <- rbind(xyzg, cbind(connectedPix@coords, rep(ind + numTif, numPts), connectedPix@data$layer))
  }
  
}

rm(matConda2)
rm(babyConda)
rm(porePts, whitePts, connectedPix)
gc()



## ----------------------------------------------------------
## read in surface python files
## ----------------------------------------------------------


setwd("D:/recon/processed/outlined")


np <- import("numpy")


# data reading
condaMatrix1 <- np$load("reconSurfacePixels1.npy")

xyzg <- matrix(NA, nrow = 1, ncol = 4)

numTif <- dim(condaMatrix1)[1]

for(ind in 1:numTif){
  
  babyConda <- raster(condaMatrix1[ind,,])
  
  
  whitePts <- rasterToPoints(babyConda, function(x){x!=0})
  
  if(dim(whitePts) != 0){
    
    print(ind)
    
    ridX <- (names(which(table(whitePts[,1]) > 100)))
    porePts <- whitePts[-which(whitePts[,1] %in% ridX),]
    
    
    connectedPix <- extract(babyConda, SpatialPoints(porePts[,1:2]), sp = T)
    
    numPts <- length(connectedPix@data$layer)
    
    xyzg <- rbind(xyzg, cbind(connectedPix@coords, rep(ind, numPts), connectedPix@data$layer))
  }
  
}

rm(condaMatrix1)
gc()

condaMatrix2 <- np$load("reconSurfacePixels2.npy")


for(ind in 1:dim(condaMatrix2)[1]){
  
  babyConda <- raster(condaMatrix2[ind,,])
  
  
  whitePts <- rasterToPoints(babyConda, function(x){x!=0})
  
  if(dim(whitePts) != 0){
    
    print(ind+numTif)
    
    ridX <- (names(which(table(whitePts[,1]) > 100)))
    porePts <- whitePts[-which(whitePts[,1] %in% ridX),]
    
    
    connectedPix <- extract(babyConda, SpatialPoints(porePts[,1:2]), sp = T)
    
    numPts <- length(connectedPix@data$layer)
    
    xyzg <- rbind(xyzg, cbind(connectedPix@coords, rep(ind + numTif, numPts), connectedPix@data$layer))
  }
  
}

rm(condaMatrix2)
rm(babyConda)
rm(porePts, whitePts, connectedPix)
gc()




saveRDS(xyzg, file = "reconSurfacePixels.rds")

























ctab = color.scale(xyzg[,4])

scatter3d(xyzg[,1], xyzg[,3], xyzg[,2], surface = F, point.col = ctab)


# shringake pores

fewPores <- which(xyzg[,4] == 62 | xyzg[,4] == 1)

scatter3d(xyzg[fewPores,1], xyzg[fewPores,2], xyzg[fewPores,3], surface = F, point.col="black")



# crop out small small pores

bigPoreNames<-names(which(table(xyzg[,4]) > 100))

bigPoreNames <- bigPoreNames[-which(bigPoreNames=="118" | bigPoreNames =="188")] # for now, crop out the false pore which is the biggest also

bigPores <- which(xyzg[,4] %in% bigPoreNames)

ctab <- color.scale(xyzg[bigPores,4])

scatter3d(xyzg[bigPores,1], xyzg[bigPores,2], xyzg[bigPores,3], surface = F, point.col = ctab)


## looking at the false pores...really need to figure out how to crop
## those verticle white bars BEFORE python processing

biggestPore <- which(xyzg[,4] == 118)

scatter3d(xyzg[biggestPore,1], xyzg[biggestPore,2], xyzg[biggestPore,3], surface = F, point.col = "black")

which.max(table(xyzg[bigPores,4]))

secondBiggestPore <- which(xyzg[,4]==188)

scatter3d(xyzg[secondBiggestPore,1], xyzg[secondBiggestPore,2], xyzg[secondBiggestPore,3], surface = F, point.col = "black")

