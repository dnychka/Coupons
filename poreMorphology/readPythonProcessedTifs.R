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
## read in volume python files (recon)
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





## ----------------------------------------------------------
## read in volume python files (scoutscan)
## ----------------------------------------------------------

setwd("D:/scoutscan/processed")


np <- import("numpy")


# data reading
matConda <- np$load("scoutscanVolumePixels.npy")

xyzg <- matrix(NA, nrow = 1, ncol = 4)

numTif <- dim(matConda)[1]

for(ind in 1:numTif){
  
  babyConda <- raster(matConda[ind,,])
  
  
  whitePts <- rasterToPoints(babyConda, function(x){x!=0})
  
  if(dim(whitePts)[1] != 0){
    
    print(ind)
    
    # ridX <- (names(which(table(whitePts[,1]) > 100)))
    # porePts <- whitePts[-which(whitePts[,1] %in% ridX),]
    # 
    
    connectedPix <- extract(babyConda, SpatialPoints(whitePts[,1:2]), sp = T)
    
    numPts <- length(connectedPix@data$layer)
    
    xyzg <- rbind(xyzg, cbind(connectedPix@coords, rep(ind, numPts), connectedPix@data$layer))
  }
  
}

rm(matConda)
gc()

xyzg <- xyzg[-1,]



## ----------------------------------------------------------
## read in surface python files (scoutscan)
## ----------------------------------------------------------

setwd("D:/scoutscan/processed")


np <- import("numpy")


# data reading
matConda <- np$load("scoutscanSurfacePixels.npy")

xyzg <- matrix(NA, nrow = 1, ncol = 4)

numTif <- dim(matConda)[1]

for(ind in 1:numTif){
  
  babyConda <- raster(matConda[ind,,])
  
  
  whitePts <- rasterToPoints(babyConda, function(x){x!=0})
  
  if(dim(whitePts)[1] != 0){
    
    print(ind)
    
    # ridX <- (names(which(table(whitePts[,1]) > 100)))
    # porePts <- whitePts[-which(whitePts[,1] %in% ridX),]
    # 
    
    connectedPix <- extract(babyConda, SpatialPoints(whitePts[,1:2]), sp = T)
    
    numPts <- length(connectedPix@data$layer)
    
    xyzg <- rbind(xyzg, cbind(connectedPix@coords, rep(ind, numPts), connectedPix@data$layer))
  }
  
}

rm(matConda)
gc()

xyzg <- xyzg[-1,]


ctab <- color.scale(xyzg[,4])

scatter3d(xyzg[,1], xyzg[,2], xyzg[,3], surface = F, point.col = ctab)
