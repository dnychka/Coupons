## Dani Barna
## 
##
## Details: Test coupons, Inconel 718 Plate 02-Build 01
##
##
## Working document for creating histograms from function radDist
## in file "radialDistance.R" 
##
## Also plots test coupons for manual look
##
##
##------------------------------------------------------------------------------------------------------------------------

library(rgl)
library(fields)
library(plyr)
library(dplyr)
library(stats)


setwd("C:/Users/barna/Documents/Inconel718_73")

poreData <- readRDS("porosityData.rds")

prepData <- readRDS("preparationData.rds")


## As of 4/18, we just want to look at test coupons 
## that have the 0-0 orientation. These are indexed at
## 4 12 15 17 22 30 31 37 41 43 44 45 51 54 56 57

## ---------------------------------------------
## Not croppping the coupon yet
## ---------------------------------------------

n = 57

# plot coupon
plot3d(poreData[[n]]$comX, 
       poreData[[n]]$comY, 
       poreData[[n]]$comZ, 
       type = "s", size = 0.35,
       xlab = " ", ylab = " ", zlab = " ")


# scale and center the pore coordinates 
scaleX <- scale(poreData[[n]]$comX, scale = TRUE, center = TRUE)
scaleY <- scale(poreData[[n]]$comY, scale = TRUE, center = TRUE)
scaleZ <- scale(poreData[[n]]$comZ, scale = TRUE, center = TRUE)

centroid <- c(mean(scaleX), mean(scaleY), mean(scaleZ))

# then call radDist
source("radialDistance.R")
histData <- radDist(scaleX, scaleY, scaleZ)

# create the vertical axis histogram
h <- hist(histData[[1]]$vertDist, freq = FALSE,breaks = seq(0,3.5,by = 0.05),
          main = "",
          xlab = "scaled radial distance",
          col = "grey35")
mtext("radial pore distance to vertical centroid line", side=3, adj=0, line=1.6, cex=1.4, font=1)
mtext(paste0("test coupon index ", n, " 0-0"), side=3, adj=0, line=0.3, cex=1, font=1)


# create the principle component axis histogram
h <- hist(histData[[1]]$pcaDist, freq = FALSE,breaks = seq(0,3.5,by = 0.05),
          main = "",
          xlab = "scaled radial distance",
          col = "grey35")
mtext("radial pore distance to first principle component", side=3, adj=0, line=1.6, cex=1.4, font=1)
mtext(paste0("test coupon index ", n, " 0-0"), side=3, adj=0, line=0.3, cex=1, font=1)

# plot both the vertical centroid line and the PC
plot3d(scaleX,
       scaleY,
       scaleZ,
       type = "s", size = 0.35,
       xlab = " ", ylab = " ", zlab = "z")
lines3d(3*c(centroid[1], centroid[1] + histData[[2]][1]), 
        3*c(centroid[2], centroid[2] + histData[[2]][2]), 
        3*c(centroid[3], centroid[3] + histData[[2]][3]), 
        col="tomato", lwd = 3)
spheres3d(centroid[1],centroid[2], centroid[3], radius = 0.05, col = "magenta")
lines3d(-3*c(centroid[1], centroid[1] + histData[[2]][1]), 
        -3*c(centroid[2], centroid[2] + histData[[2]][2]), 
        -3*c(centroid[3], centroid[3] + histData[[2]][3]), 
        col="tomato", lwd = 3)
lines3d(c(centroid[1],centroid[1]),
        c(centroid[2],centroid[2]),
        c(centroid[3]-2.5,centroid[3]+2),
        col = "cornflowerblue", lwd = 3)


