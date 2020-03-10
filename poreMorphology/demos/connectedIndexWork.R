##
##
##
##
##
##
##
## test out the connected pixels method
## ------------------------------------

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
## extract only white pixels
getSurface <- function(ind) {
  allPts <- rasterToPoints(rasterStack[[ind]], function(x){x==255})[,-3]
  surfacePts <- allPts
  
  # get rid of x border
  ridX <- as.numeric(names(which(table(surfacePts[,1]) > 50)))
  surfacePts <- surfacePts[-which(surfacePts[,1] %in% ridX),]
  
return(surfacePts)}


## check connectivity index of white pixels
isConnected <- function(ind){
  
  
  getCoordsGetClusters <- function(ind){
  # get rid of x border
  whitePts <- rasterToPoints(rasterStack[[ind]], function(x){x==255})[,-3]
  ridX <- as.numeric(names(which(table(whitePts[,1]) > 50)))
  porePts <- whitePts[-which(whitePts[,1] %in% ridX),]
  
  porePts <- cbind(porePts, rep(255, length(porePts[,1])))
  
  poreRaster <- rasterFromXYZ(porePts)
  
  # identify connected white pixels
  connectedRaster <- clump(poreRaster, directions = 8)
  
  connectedPix <- extract(connectedRaster, SpatialPoints(porePts[,1:2]), sp = T)
  
  return(connectedPix)}
  
  
  
  
  top <- getCoordsGetClusters(ind)
  bottom <- getCoordsGetClusters(ind-1)
  
  
  ## if a top clump shares any coords with a bottom clump
  for(j in 1:length(unique(bottom@data$clumps))){
    
    isNew = 1
    
    for(k in 1:length(unique(top@data$clumps))){
      
      testSharing <- sum(top@coords[top@data$clumps == k,1] %in% bottom@coords[bottom@data$clumps == j,1] & top@coords[top@data$clumps == k,2] %in% bottom@coords[bottom@data$clumps == j,2])
      
      print(paste0("****",k," with ****", j))
      print(testSharing)
      
      if(testSharing != 0){
        
        xyB <- bottom@coords[bottom@data$clumps == j,]
        
        xyT <- top@coords[top@data$clumps == k,]
        
        xyz <- rbind(cbind(xyB, rep(ind-1, length(xyB))), cbind(xyT, rep(ind,length(xyT))))
        

        xyzc <- cbind(xyz, rep(j, length(xyz))) 
      }
    }
  
  
    
  }
  
  plot(top@coords, col = 'violetred1', pch = 15, cex=.5)
  points(bottom@coords, col = alpha('cornflowerblue', .4), pch = 15, cex=1)
  
  
  
}


plot(surfacePts,pch=15, cex=0.3)

idmaker <- function(x){
  max.val = x*100
  count <- nchar(as.character(max.val))                       # find out how many 'numbers' each ID will have after the letter
  size <- paste("%0",count,"d",sep="")                        # set the variable to be fed into 'sprintf' to ensure we have leading 0's
  lets <- toupper(sample(letters,x, replace=T))  
  nums <- sprintf(size, sample(1:max.val)[1])         # randomising the letters 
  ids <- paste(lets,nums,sep="")  

  return(ids)}


## -----------------------------------
## working
## -----------------------------------

fileNam <- "D:/recon/processed"


setwd(fileNam)

rasterList <- mixedsort(list.files(full.names = TRUE))

rasterStack <- stack(rasterList[-1])





X <- vector()
Y <- vector()
Z <- vector()

startTime <- Sys.time()

# for(i in 1:dim(rasterStack)[3]){

for(i in 500:700){
  print(i)
  
  slice <- getSurface(i)
  
  X <- c(X,slice[,1])
  Y <- c(Y,slice[,2])
  Z <- c(Z, rep(i, length(slice[,1])))
}
endTime <- Sys.time()

print(endTime-startTime)















