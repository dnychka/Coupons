##
##
##
##
##
##
##
## attempt to implement bivariate kernel density estimate
## for sampling from joint distribution of theta, z
## ------------------------------------------------------

library(ks)

library(KernSmooth)

myDir <- "C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData"

setwd(myDir)

# load the sample simulated coupon
trueCoupon <- readRDS(file.path(myDir, "/simulatedCouponInternalJitter.rds")) #has NOT been nls-ed

# knock it off its axis
aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}
theta = -13*pi/180

Ry <- aboutY(theta)

poreCoordinates <- trueCoupon %*% Ry 

rm(myDir, aboutY, theta, Ry)


setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty")
source("functionsForBootstrapping.R")




##
## fit it
## -----------------------------------------


nlsObj <- nlsAxisFit(poreCoordinates)

nlsCoeff <- coef(nlsObj)

nlsCoupon <- newCoupon(poreCoordinates, nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                       nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])



##
## get the parameters from the nls obj
## -----------------------------------------

load("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData/coupon4poreCI.rda")


## model output from first round of nls
nlsTheta <- atan2(nlsCoupon[,2], nlsCoupon[,1])
nlsZ <- nlsCoupon[,3]
nlsR <- sqrt(nlsCoupon[,1]^2+nlsCoupon[,2]^2)



##
## estimate the joint distribution of z and theta
## ----------------------------------------------


## working from package kernSmooth...

zN <- scale(nlsZ)
tN <- scale(nlsTheta)

ztMat <- cbind(tN, zN)

bw = Hpi(ztMat)

kernd <- bkde2D(ztMat, bandwidth=c(bw[1,1], bw[2,2]), gridsize = c(1000,1000))

image.plot(kernd$x1, kernd$x2, kernd$fhat,
           main = "kde for z and theta, bw = hpi, grd = 1000",
           xlab = "theta (normalized)", ylab = "z (normalized)")


#compare to ks

Hpi1 <- Hpi(x = ztMat)
 
fhatPi1 <- kde(x=ztMat, H=Hpi1)

image.plot(fhatPi1$eval.points[[1]], fhatPi1$eval.points[[2]], fhatPi1$estimate) #um what


# how to sample from this (kernSmooth)




















function (x, bandwidth, gridsize = c(51L, 51L), range.x, truncate = TRUE) 
{
  if (!missing(bandwidth) && min(bandwidth) <= 0) 
    stop("'bandwidth' must be strictly positive")
  n <- nrow(x)
  M <- gridsize
  h <- bandwidth
  tau <- 3.4
  if (length(h) == 1L) 
    h <- c(h, h)
  if (missing(range.x)) {
    range.x <- list(0, 0)
    for (id in (1L:2L)) range.x[[id]] <- c(min(x[, id]) - 
                                             1.5 * h[id], max(x[, id]) + 1.5 * h[id])
  }
  a <- c(range.x[[1L]][1L], range.x[[2L]][1L])
  b <- c(range.x[[1L]][2L], range.x[[2L]][2L])
  gpoints1 <- seq(a[1L], b[1L], length = M[1L])
  gpoints2 <- seq(a[2L], b[2L], length = M[2L])
  gcounts <- linbin2D(x, gpoints1, gpoints2)
  L <- numeric(2L)
  kapid <- list(0, 0)
  for (id in 1L:2L) {
    L[id] <- min(floor(tau * h[id] * (M[id] - 1)/(b[id] - 
                                                    a[id])), M[id] - 1L)
    lvecid <- 0:L[id]
    facid <- (b[id] - a[id])/(h[id] * (M[id] - 1L))
    z <- matrix(dnorm(lvecid * facid)/bbh[id])
    tot <- sum(c(z, rev(z[-1L]))) * facid * h[id]
    kapid[[id]] <- z/tot
  }
  kapp <- kapid[[1L]] %*% (t(kapid[[2L]]))/n
  if (min(L) == 0) 
    warning("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'")
  P <- 2^(ceiling(log(M + L)/log(2)))
  L1 <- L[1L]
  L2 <- L[2L]
  M1 <- M[1L]
  M2 <- M[2L]
  P1 <- P[1L]
  P2 <- P[2L]
  rp <- matrix(0, P1, P2)
  rp[1L:(L1 + 1), 1L:(L2 + 1)] <- kapp
  if (L1) 
    rp[(P1 - L1 + 1):P1, 1L:(L2 + 1)] <- kapp[(L1 + 1):2, 
                                              1L:(L2 + 1)]
  if (L2) 
    rp[, (P2 - L2 + 1):P2] <- rp[, (L2 + 1):2]
  sp <- matrix(0, P1, P2)
  sp[1L:M1, 1L:M2] <- gcounts
  rp <- fft(rp)
  sp <- fft(sp)
  rp <- Re(fft(rp * sp, inverse = TRUE)/(P1 * P2))[1L:M1, 1L:M2]
  rp <- rp * matrix(as.numeric(rp > 0), nrow(rp), ncol(rp))
  list(x1 = gpoints1, x2 = gpoints2, fhat = rp)
}

























## example from web

# set.seed(8192)
# samp <- 200
# mus <- rbind(c(-2,2), c(0,0), c(2,-2))
# Sigmas <- rbind(diag(2), matrix(c(0.8, -0.72, -0.72, 0.8), nrow=2), diag(2))
# cwt <- 3/11
# props <- c((1-cwt)/2, cwt, (1-cwt)/2)
# x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)
# 
# 
# ## test on my data
# 
# Znormed <- scale(nlsZ)
# ThetaNormed <- scale(nlsTheta)
# 
# 
# # build cov matrix by hand
# 
# ztMat <- cbind(Znormed, ThetaNormed)
# 
# sig <- cov(ztMat)
# 
# n <- nrow(ztMat) #number of pores
# 
# #create means for each column
# M_mean <- matrix(data=1, nrow=n) %*% cbind(mean(Znormed),mean(ThetaNormed)) 
# 
# #creates a difference matrix
# D <- ztMat - M_mean
# 
# #creates the covariance matrix
# C <- (n-1)^-1 * t(D) %*% D     ## yep same as sig
# 
# 
# ## transform my data
# 
# ztTransform <- solve(sqrt(sig)) %*% ztMat
# 
# 
# # okay so maybe we don't transform maybe we just use unconstrained bandwidth
# # parameterization instead.
# 
# Hpi1 <- Hpi(x = ztMat)
# 
# fhatPi1 <- kde(x=ztMat, H=Hpi1)
# 
# plot(fhatPi1)
# 
# points(ztMat[,1], ztMat[,2])
# 
# # sample from the kde
# 
# zNew <- sample(fhatPi1$x[1,], 1, replace = F) # sample from kde pt
# 
# kernValZ <- rnorm(1, zNew, fhatPi1$H[1,1]) # draw from the kernel associated with the pt
# 
# 
# 
# 
# ztMatUnscaled <- cbind(nlsZ, nlsTheta)
# 
# Hpi2 <- Hpi(x = ztMatUnscaled)
# 
# fHatPi2 <- kde(x=ztMatUnscaled, H = Hpi2)
# 
# plot(fHatPi2)







