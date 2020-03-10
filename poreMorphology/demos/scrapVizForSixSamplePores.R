##
##
##
##
##
##
##
##
## -------------------------------------------------


library(rgl)
library(car)
library(raster)
library(rasterVis)
library(rgdal)
library(gtools)
library(fields)

library(scatterplot3d)
library(scales)



## -----------------------------------
## functions
## -----------------------------------
## extract only pixels with ones values
getSurface <- function(ind) {
  allPts <- rasterToPoints(rasterStack[[ind]], function(x){x==255})[,-3]
  surfacePts <- allPts

  # get rid of x border
  ridX <- as.numeric(names(which(table(surfacePts[,1]) > 50)))
  
  surfacePts <- surfacePts[-which(surfacePts[,1] %in% ridX),]
  
  
  n <- nrow(surfacePts)
  surfacePts <- surfacePts[sample(n, n, replace=FALSE),]
  

  return(surfacePts)}


## ----------------------------------
## Construct and save surface data
## ----------------------------------

fileNam <- "D:/recon/processed/outlined"

  
  setwd(fileNam)
  
  rasterList <- mixedsort(list.files(full.names = TRUE))
  
  rasterStack <- stack(rasterList)
  
  
  
  

  X <- vector()
  Y <- vector()
  Z <- vector()
  
  startTime <- Sys.time()
  
  for(i in 1:dim(rasterStack)[3]){
    
  #for(i in 500:700){
    print(i)
    
    slice <- getSurface(i)
    
    X <- c(X,slice[,1])
    Y <- c(Y,slice[,2])
    Z <- c(Z, rep(i, length(slice[,1])))
  }
  endTime <- Sys.time()
  
  print(endTime-startTime)
  
  
  
  
  
  
  setwd("C:/Users/barna/Documents/Coupons/poreMorphology")
  
  
load("fulDenseAll.rda")

scatter3d(X,Z,Y, surface = F, point.col = "black")

#scatter3d(X,Z,Y, pch = ".", surface=FALSE, point.col = ifelse(Z == 475 & X < 150 & Y > 500 & Y < 600, "tomato","black"))

  
load("sixSamplePoresColorCoded.rda")

ctab <- color.scale(allDat[,4], col = c('#FD0E35', '#FF00CC','#66FF66', '#5946B2', '#58427C',  '#2243B6'))

scatter3d(allDat[,1], allDat[,2], allDat[,3], surface = F, point.col = ctab) 

c(sphRed, sphBlue, sphviolet)


alreadyHaveThese <- which(X %in% allDat[,1] & Y %in% allDat[,2] & Z %in% allDat[,3])

D <- rep(0, length(X[-alreadyHaveThese]))

allWithColors <- rbind(cbind(X[-alreadyHaveThese],Y[-alreadyHaveThese],Z[-alreadyHaveThese],D), allDat)

ctabAll <- color.scale(allWithColors[,4], col = c('black', '#FD0E35', '#FF00CC','#66FF66', '#5946B2', '#58427C',  '#2243B6'))

scatter3d(allWithColors[,1], allWithColors[,2], allWithColors[,3], surface = F, point.col = ctabAll)









## get volume of our six colored sample pores
## -----------------------------------------

keep <- which(X > 100 & X < 300 & Y > 100 & Y < 300)


Xv <- X[keep]; Yv <- Y[keep]; Zv <- Z[keep]

load("volumeSixSamplePores.rda")


scatter3d(X,Y,Z, surface = F, point.col = "black")


# now group the volume coords based on the already clustered surfaces

redCluster <- which(allDat[,4]==1)

redVol <- which(Xv <= max(allDat[redCluster,1]) & Xv >= min(allDat[redCluster,1]) & Yv <= max(allDat[redCluster,2]) & Yv >= min(allDat[redCluster,2]) & Zv <= max(allDat[redCluster,3]) & Zv >= min(allDat[redCluster,3]))

scatter3d(Xv[redVol], Yv[redVol], Zv[redVol], surface = F, point.col = "black")

testingRed <- rbind(allDat[redCluster,], cbind(Xv[redVol], Yv[redVol], Zv[redVol], rep(2, length(Xv[redVol]))))

scatter3d(testingRed[,1], testingRed[,2], testingRed[,3], surface = F, 
          point.col = ifelse(testingRed[,4] == 1, "cornflowerblue", "violetred1"),
          id = list(method= "mahal"))


DveqRed <- 2 * (3 * length(Xv[redVol]) / (4 * pi) )^(1/3)

sphRed <- pi * DveqRed ^ 2 / length(allDat[redCluster,1])

pi^(1/3)*(6*length(Xv[redVol]))^(2/3)/length(allDat[redCluster,1])




blueCluster <- which(allDat[,4]==9)

blueVol <- which(Xv <= max(allDat[blueCluster,1]) & Xv >= min(allDat[blueCluster,1]) & Yv <= max(allDat[blueCluster,2]) & Yv >= min(allDat[blueCluster,2]) & Zv <= max(allDat[blueCluster,3]) & Zv >= min(allDat[blueCluster,3]))

scatter3d(Xv[blueVol], Yv[blueVol], Zv[blueVol], surface = F, point.col = "black")

testingBlue <- rbind(allDat[blueCluster,], cbind(Xv[blueVol], Yv[blueVol], Zv[blueVol], rep(2, length(Xv[blueVol]))))

scatter3d(testingBlue[,1], testingBlue[,2], testingBlue[,3], surface = F, 
          point.col = ifelse(testingBlue[,4] == 9, "cornflowerblue", "violetred1"),
          id = list(method= "mahal"))


DveqBlue <- 2 * (3 * length(Xv[blueVol]) / (4 * pi) )^(1/3)

sphBlue <- pi * DveqBlue / length(allDat[blueCluster,1])



violetCluster <- which(allDat[,4]==6)

violetVol <- which(Xv <= max(allDat[violetCluster,1]) & Xv >= min(allDat[violetCluster,1]) & Yv <= max(allDat[violetCluster,2]) & Yv >= min(allDat[violetCluster,2]) & Zv <= max(allDat[violetCluster,3]) & Zv >= min(allDat[violetCluster,3]))

scatter3d(Xv[violetVol], Yv[violetVol], Zv[violetVol], surface = F, point.col = "black")

testingViolet <- rbind(allDat[violetCluster,], cbind(Xv[violetVol], Yv[violetVol], Zv[violetVol], rep(2, length(Xv[violetVol]))))

scatter3d(testingViolet[,1], testingViolet[,2], testingViolet[,3], surface = F, 
          point.col = ifelse(testingViolet[,4] == 6, "cornflowerblue", "violetred1"),
          id = list(method= "mahal"))


Dveqviolet <- 2 * (3 * length(Xv[violetVol]) / (4 * pi) )^(1/3)

sphviolet <- pi * Dveqviolet / length(allDat[violetCluster,1])



pi^(1/3)*(6*length(Xv[violetVol]))^(2/3) / length(allDat[violetCluster,1])








  
c(sphRed, sphBlue, sphviolet)




















  coords <- data.frame(X[get],Z[get],Y[get])
  
  bigDist <- dist(coords, 'euclidean')
  
  hCLustAvgBig <- hclust(bigDist, method = 'average')
  
  gsBig <- cutree(hCLustAvgBig, 40)
  
  scatter3d(X[get], Z[get], Y[get], surface = F, point.col = gsBig)
  
  dat <- cbind(X[get], Y[get], Z[get])
  
  clusterSize <- table(gsBig)
  
  smalls <- quantile(clusterSize, .75)

  dropThese <- names(clusterSize)[which(clusterSize < smalls)]

  dat <- dat[-which(gsBig %in% dropThese),]
  
  rm(bigDist, hCLustAvgBig)

  distMat <- dist(dat, 'euclidean')  
  
  hClustAvg <- hclust(distMat, method = 'average')
  
  gs <- cutree(hClustAvg, 7)
  
  ctab <- color.scale(gs)
  
  scatter3d(dat[,1], dat[,2], dat[,3], surface = F, point.col = ctab)
  
  dat <- dat[keepDemo,]
  
  distMat <- dist(dat, 'euclidean')
  
  hClustAvg <- hclust(distMat, method = 'average')
  
  gs <- cutree(hClustAvg, 3)
  
  ctab <- color.scale(gs, col = c('#e41a1c','#377eb8','darkorange'))
  
  scatter3d(dat[,1], dat[,2], dat[,3], surface = F, point.col = ctab)
  
  
  roundBlobs <- dat[which(gs == 2),]

  distMat <- dist(roundBlobs, 'euclidean')  
  
  hClustAvg <- hclust(distMat, method = 'average')

  gsRound <- cutree(hClustAvg, 4) + 5

  scatter3d(roundBlobs[,1], roundBlobs[,2], roundBlobs[,3], surface = F, point.col = gsRound)  
  
  
  allDat <- rbind(cbind(dat[-which(gs == 2),],gs[-which(gs==2)]), cbind(roundBlobs, gsRound))

  
  ctab <- color.scale(allDat[,4], col = c('#FD0E35', '#FF00CC','#66FF66',    '#5946B2', '#58427C',  '#2243B6'))
  
  scatter3d(allDat[,1], allDat[,2], allDat[,3], surface = F, point.col = ctab) 

  
  save(allDat, file = "sixSamplePoresColorCoded.rda")
  
    
  
  
  
