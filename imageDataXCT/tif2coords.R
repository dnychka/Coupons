##
## 5 June 2019
##
##
## find coords of .tif binary image
## ------------------------------------------------------------------------------


library(dplyr)
library(imager)
library(rgl)

setwd("C:/Users/barna/Documents/Coupons/imageDataXCT")


##
## test single image
## ------------------------------------------------------
fish <- load.image("C:/Users/barna/Documents/Coupons/imageDataXCT/fishBW.tif")
plot(fish)

fishDF <- as.data.frame(fish)

unique(fishDF$value) #check it's actually bianry image

outline <- which(fishDF$value==0)

outlineCoords <- fishDF[outline,]

matplot(outlineCoords$x, outlineCoords$y, type = "p", pch = 20, cex = 0.4)

zcoords <- seq(0,1,length.out = 10)

xyzcoords <- cbind(rep(outlineCoords$x, 10), rep(outlineCoords$y,10), rep(zcoords, each = 1750))

plot3d(xyzcoords[,1], xyzcoords[,2], xyzcoords[,3], type = "s", size = 0.45)

## ------------------------------------------------------




##
## read in stacked .tif
## ------------------------------------------------------

fish <- load.image("C:/Users/barna/Documents/Coupons/imageDataXCT/stackedTIFF/fishBW1.tif")
plot(fish)

fishDF <- as.data.frame(fish)

unique(fishDF$value) 

outline <- which(fishDF$value==0)

outlineCoords <- fishDF[outline,]

surfaceCoords <- cbind(outlineCoords$x, outlineCoords$y, rep(1, length(outlineCoords$x)))

for (i in 2:10){
  
  fish <- load.image(paste0("C:/Users/barna/Documents/Coupons/imageDataXCT/stackedTIF/fishBW",i,".tif"))
  
  fishDF <- as.data.frame(fish)
  
  outline <- which(fishDF$value==0)
  
  outlineCoords <- fishDF[outline,]
  
  surfaceCoords <- rbind(surfaceCoords, cbind(outlineCoords$x, outlineCoords$y, rep(i, length(outlineCoords$x))))
}


plot3d(surfaceCoords[,1], surfaceCoords[,2], surfaceCoords[,3], type = "s", size = 0.35)


