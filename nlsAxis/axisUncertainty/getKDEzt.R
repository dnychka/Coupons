##
##
##
##
##
##
##
##
## function to return z, theta smaples from the
## bivariate kde. NOTE: needs accurate bw and grd
## values that are selected from data viz provided 
## in file bivariateKDEviz.R
## ------------------------------------------------

library(KernSmooth)


getKDEzt <- function(nlsCoupon, bw, grX, grY){
  
  
  ## model output from first round of nls
  nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
  nlsZ <- nlsCoupon[,3]
  
  
  ##
  ## estimate the joint distribution of z and theta
  ## ----------------------------------------------
  
  
  ## working from package kernSmooth...
  
  zN <- scale(nlsZ)
  tN <- scale(nlsTheta)
  
  ztMat <- cbind(tN, zN)
  
  ## triple the data to approximate periodicity
  
  shift <- max(tN) - min(tN)
  
  ztMatBig <- rbind( cbind(tN - shift, zN), cbind(tN, zN), cbind(tN + shift, zN))
  
  
  kernd <- bkde2D(ztMatBig, bandwidth = bw, gridsize = c(grX,grY)) # note hardcoded bw & grdsize
  
  ## chuck off the ends
  
  kdeT <- kernd$x1[which(kernd$x1 <= max(tN) & (kernd$x1 >= min(tN))) ]
  kdeZ <- kernd$x2
  
  kdeFhat <- kernd$fhat[which(kernd$x1 <= max(tN) & kernd$x1 >= min(tN)),]
  
  
  # how to sample from this (strategy: Nychka email post Thanksgiving break)
  
  totalDensity <- sum(kdeFhat)
  
  kerndProb <- (kdeFhat) / totalDensity
  
  
  
  numResamp <- 5000
  
  test <- sample(kdeFhat, numResamp, replace = T, prob = kerndProb)
  
  
  indZ <- vector()
  indT <- vector()
  
  for(i in 1:numResamp){
    indT[i] <- which(kdeFhat == test[i], arr.ind=TRUE)[1]
    indZ[i] <- which(kdeFhat == test[i], arr.ind=TRUE)[2]
  }
  
  sampledT <- kdeT[indT]
  sampledZ <- kdeZ[indZ]
  
  # descale Z and theta
  
  sampledT <- sampledT * attr(tN, 'scaled:scale') + attr(tN, 'scaled:center')
  sampledZ <- sampledZ * attr(zN, 'scaled:scale') + attr(zN, 'scaled:center')
  
  
  return(list(sampledT, sampledZ))
}  
  
  