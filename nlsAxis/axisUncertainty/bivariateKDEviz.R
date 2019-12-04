##
##
##
##
##
##
##
## implements bivariate kernel density estimate
## for sampling from joint distribution of theta, z
## ------------------------------------------------------

#library(ks) #inly needed for hpi bandwidth (not fully understood)

library(scales)

library(KernSmooth)

library(fields)



load("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData/coupon4poreCI.rda")

rm(bootCoefCI, bootNlsCoef, bootRadius, bootRadiusCI)



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

image.plot(kernd$x1, kernd$x2, kernd$fhat,
           main = "kde for z and theta, bw = hpi, grd = 1000",
           xlab = "theta (normalized)", ylab = "z (normalized)")

 xline(c(max(tN), min(tN)), col = "white", lty=2)

  ## chuck off the ends
  
  kdeT <- kernd$x1[which(kernd$x1 <= max(tN) & (kernd$x1 >= min(tN))) ]
  kdeZ <- kernd$x2
  
  kdeFhat <- kernd$fhat[which(kernd$x1 <= max(tN) & kernd$x1 >= min(tN)),]

# check it

image.plot(kdeT, kdeZ, (kdeFhat))

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



par(mfrow=c(1,2))
image.plot(kdeT, kdeZ, kdeFhat,
           main = "",
           xlab = "theta (normalized)", ylab = "z (normalized)")

plot(kdeT[indT], kdeZ[indZ], pch = 20, col = alpha("black", 0.2))



  

## source code for kde smoothing from kernSmooth
# function (x, bandwidth, gridsize = c(51L, 51L), range.x, truncate = TRUE) 
# {
#   if (!missing(bandwidth) && min(bandwidth) <= 0) 
#     stop("'bandwidth' must be strictly positive")
#   n <- nrow(x)
#   M <- gridsize
#   h <- bandwidth
#   tau <- 3.4
#   if (length(h) == 1L) 
#     h <- c(h, h)
#   if (missing(range.x)) {
#     range.x <- list(0, 0)
#     for (id in (1L:2L)) range.x[[id]] <- c(min(x[, id]) - 
#                                              1.5 * h[id], max(x[, id]) + 1.5 * h[id])
#   }
#   a <- c(range.x[[1L]][1L], range.x[[2L]][1L])
#   b <- c(range.x[[1L]][2L], range.x[[2L]][2L])
#   gpoints1 <- seq(a[1L], b[1L], length = M[1L])
#   gpoints2 <- seq(a[2L], b[2L], length = M[2L])
#   gcounts <- linbin2D(x, gpoints1, gpoints2)
#   L <- numeric(2L)
#   kapid <- list(0, 0)
#   for (id in 1L:2L) {
#     L[id] <- min(floor(tau * h[id] * (M[id] - 1)/(b[id] - 
#                                                     a[id])), M[id] - 1L)
#     lvecid <- 0:L[id]
#     facid <- (b[id] - a[id])/(h[id] * (M[id] - 1L))
#     z <- matrix(dnorm(lvecid * facid)/bbh[id])
#     tot <- sum(c(z, rev(z[-1L]))) * facid * h[id]
#     kapid[[id]] <- z/tot
#   }
#   kapp <- kapid[[1L]] %*% (t(kapid[[2L]]))/n
#   if (min(L) == 0) 
#     warning("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'")
#   P <- 2^(ceiling(log(M + L)/log(2)))
#   L1 <- L[1L]
#   L2 <- L[2L]
#   M1 <- M[1L]
#   M2 <- M[2L]
#   P1 <- P[1L]
#   P2 <- P[2L]
#   rp <- matrix(0, P1, P2)
#   rp[1L:(L1 + 1), 1L:(L2 + 1)] <- kapp
#   if (L1) 
#     rp[(P1 - L1 + 1):P1, 1L:(L2 + 1)] <- kapp[(L1 + 1):2, 
#                                               1L:(L2 + 1)]
#   if (L2) 
#     rp[, (P2 - L2 + 1):P2] <- rp[, (L2 + 1):2]
#   sp <- matrix(0, P1, P2)
#   sp[1L:M1, 1L:M2] <- gcounts
#   rp <- fft(rp)
#   sp <- fft(sp)
#   rp <- Re(fft(rp * sp, inverse = TRUE)/(P1 * P2))[1L:M1, 1L:M2]
#   rp <- rp * matrix(as.numeric(rp > 0), nrow(rp), ncol(rp))
#   list(x1 = gpoints1, x2 = gpoints2, fhat = rp)
# }
# 
# 
# 
# 
# 
# 
# 
# 
