

library(fractaldim)
library(RandomFields)
library(fields)
library(raster)
library(fda)

# 2d random fields
n <- 128 
rf2d <- GaussRF(x = c(0,1, 1/n), y = c(0,1, 1/n), 
                model = "stable", grid = TRUE, gridtriple = TRUE, 
                param = c(mean=0, variance=1, nugget=0, scale=1, kappa=1)) 
par(mfrow=c(1,3)) 
fd.estim.isotropic(rf2d, p.index = 1, direction='hv', plot.loglog = TRUE, plot.allpoints = TRUE)
fd.estim.squareincr(rf2d, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE) 
fd.estim.filter1(rf2d, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE)

dev.off()

par(mfrow=c(1,2))

smoothSurf <- "D:/server-files/SURFACES/O4.tif"

smoothest <- raster(smoothSurf)

smoothestMat <- as.matrix(smoothest)[250:750,250:750]

fd.estim.squareincr(smoothestMat, p.index = 1, 
                   plot.loglog = T, plot.allpoints = T)


roughSurf <- "D:/server-files/SURFACES/O22.tif"

roughest <- raster(roughSurf)

roughestMat <- as.matrix(roughest)[250:750,250:750]

fd.estim.squareincr(roughestMat, p.index = 1, 
                   plot.loglog = T, plot.allpoints = T)


