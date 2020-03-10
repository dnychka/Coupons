##
##
##
##
##
##
## ---------------------------------------------

library(fields)
library(car)
library(scatterplot3d)
library(dplyr)

setwd("C:/Users/barna/Documents/Coupons/poreMorphology/morphologyData")
vPoresAll <- readRDS("reconVolumePixels.rds")
sPoresAll <- readRDS("reconSurfacePixels.rds")


## -----------------
## subset to only big pores
## -----------------

## do this by volume: anything with less than 100 pixels of volume is chucked out
bigPoreNames<-names(which(table(vPoresAll[,4]) > 100))

# for now, crop out the false pore which is the biggest (vertical line artifact)
bigPoreNames <- bigPoreNames[-which(bigPoreNames=="118" | bigPoreNames =="188")] 
bigPores <- which(vPoresAll[,4] %in% bigPoreNames)


vPores <- vPoresAll[bigPores,]

sPores <- matrix(NA, nrow = 1, ncol = 4)


for(i in bigPoreNames){
  
  print(i)
  
  getOne <- which(vPores[,4] == i)
  
  surfClumps <- which(sPoresAll[,1] %in% vPores[getOne,1] & sPoresAll[,2] %in% vPores[getOne,2] & sPoresAll[,3] %in% vPores[getOne,3])
  
  sPoresAll[surfClumps,4] = rep(as.numeric(i), length(sPoresAll[surfClumps,4]))
  
  sPores <- rbind(sPores, sPoresAll[surfClumps,])
  
    
  rm(getOne)
}



sPores <- readRDS("reconSurfacePixels_onlyBigPores.rds")

vsPores <- readRDS("reconVolumePixels_onlyBigPores.rds")


colnames(sPores) <- c("x","y","z","g")
colnames(vsPores) <- c("x","y","z","g")


x <- data.table(sPores[which(sPores[,4] == 1),])
v <- data.table(vsPores[which(vsPores[,4] == 1),])

sDT <- data.table(sPores)
vsDT <- data.table(vsPores)

vsDT <- data.table(cbind(vsDT, 0), key=c(c("x","y","z","g"),paste0("V", 2)))

setnames(vsDT, c(head(names(vsDT), -1L), "found"))
isIn <- vsDT[sDT, list(found=ifelse(is.na(found), 0, 1))] #gives us those obnoxious floaty surface pixels otuside volumes

sPores <- sPores[which(isIn$found == 1),] # take out obnoxious floaty pixels

## now ALL pixels in sPores can be found in vPores


# annnnd a second work-around to get only the volume pixels (without the surface pixels) from each pore

bigMat <- rbind_list(anti_join(vsDT, sDT))

## now, we should have vPores and sPores: vPores has only internal pixels and sPores has only surface pixels



b <- as.matrix(bigMat)
b <- b[,-5]

vPxls <- b
sPxls <- sPores

save(vPxls, sPxls, file = "reconVolumeAndSurfacePixels.rda")


dat <- rbind(b, sPores)

ctab <- c(color.scale(b[,4]) , rep("grey80", length(sPores[,4])))

## check

scatter3d(dat[,1], dat[,2], dat[,3], surface = F, point.col = ctab)

scatter3d(bigMat$x, bigMat$y, bigMat$z, surface = F, point.col = "black")

