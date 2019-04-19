## Dani Barna
## 
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## Return radial pore distance as measured by 
## (1) distance to first principle component and 
## (2) vertical axis through centroid
## can use distances to create histograms
##
## also returns  first eigenvector
##
##------------------------------------------------------------------------------------------------------------------------

radDist <- function(comX, comY, comZ){
  
  # caluclate the centroid for the test coupon
  centroid <- c(mean(scaleX), mean(scaleY), mean(scaleZ))
  
  # function to calculate distance
  eucDist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  # caluculate the distance to veritcal centroid axis
  radialDist <- rep(NA, length(scaleZ))
  for(i in 1:length(scaleZ)){
    radialDist[i] <- eucDist(c(centroid[1],centroid[2]), c(scaleX[i], scaleY[i]))
  }
  
  # PCA on the scaled and centered 
  xyz <- data.frame(matrix(c(scaleX, scaleY, scaleZ), byrow = FALSE, ncol = 3))
  pr <- prcomp(xyz)
  
  # find x and y coordinates for the line that makes up the
  # first principle component
  xpr <- pr$rotation[1,1]*(scaleZ - centroid[3])/pr$rotation[3,1] + centroid[1]
  ypr <- pr$rotation[2,1]*(scaleZ - centroid[3])/pr$rotation[3,1] + centroid[2]
  
  # calc distance to first principle component
  radialDistPCA <- rep(NA, length(scaleZ))
  for(i in 1:length(scaleZ)){
    radialDistPCA[i] <- eucDist(c(xpr[i], ypr[i]), c(scaleX[i], scaleY[i]))
  }
  
  distanceInfo <- data.frame(vertDist = radialDist,
                             pcaDist = radialDistPCA)
  
  firstEig = pr$rotation[,1]
  
  info = list(distanceInfo, firstEig)
  
  return(info)
}