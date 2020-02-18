##
##
##
##
##
##
##
## KDE visualizations
## ------------------------------------------------------

library(scales)

library(KernSmooth)

library(fields)

library(rgl)

library(reshape2)

library(ggplot2)

library(extrafont)
loadfonts(device = "win")


##
## implement KDE sampling
## ----------------------------------------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped")

n = 11 #coupon F3

load(paste0("./nlsCoupon", n, ".rda"))
rm(oldCoupon, nlsCoeff)


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

  bw = 0.04
  grX = 900
  grY = 300

  kernd <- bkde2D(ztMatBig, bandwidth = bw, gridsize = c(grX,grY)) # note hardcoded bw & grdsize
  

          # image.plot(kernd$x1, kernd$x2, kernd$fhat,
          #            main = "kde for z and theta, bw = hpi, grd = 1000",
          #            xlab = "theta (normalized)", ylab = "z (normalized)")
          # 
          #  xline(c(0, shift), col = "white", lty=2)
          #  xline(c(min(tN), max(tN)), col = "white")

  ## chuck off the ends
  
  kdeT <- kernd$x1[which(kernd$x1 <= shift & (kernd$x1 >= 0)) ]
  kdeZ <- kernd$x2
  
  kdeFhat <- kernd$fhat[which(kernd$x1 <= shift & kernd$x1 >= 0),]
  
  
## -------------------------------
## make the ggplot
## -------------------------------
  
longData <- melt(kdeFhat)
  
names(longData) <- c("t", "z", "density")
  
ctab = tim.colors(dim(longData)[1])

tLabels <- c(quantile(kdeT, probs = c(0, 0.5, 1)))
tLabels <- round(tLabels, 1)
tLabels <- as.character(tLabels)

zLabels <- c(quantile(kdeZ, probs = c(0, 0.5, 1))) 
zLabels <- round(zLabels, 1)
zLabels <- as.character(zLabels)
  
ggplot(longData, aes(x = t, y = z)) + 
  geom_raster(aes(fill=density)) + 
  scale_fill_gradientn(colors = ctab)+
  theme_minimal() +
  scale_x_continuous(name = "theta (normalized)",
                   breaks = c(0, 150, 300),
                   labels = tLabels)+
  scale_y_continuous(name = "z (normalized)",
                     breaks = c(0, 150, 300),
                     labels = zLabels) +
  ggtitle("Kernel Density Estimate for Coupon F3")+
  theme(plot.title = element_text(size = 20),
                            axis.text=element_text(size=16),
                            axis.title.x = element_text(size=16),
                            axis.title.y = element_text(size=16),
                            legend.text=element_text(size=16),
                            legend.title = element_text(size=16),
                            text=element_text(size=16,  family="Palatino Linotype"))


  
## ------------------------------





ggplot( NULL ) + 
  geom_raster( data = kdf , aes( t , z , fill = density ) ) +
  scale_fill_gradientn(colors = ctab)+
  ggtitle("Surface Plot for Coupon F3")+
  labs(x = "theta (normalized)", y = "z (normalized)") +
  theme_minimal()   + theme(plot.title = element_text(size = 20),
                            axis.text=element_text(size=16),
                            axis.title.x = element_text(size=16),
                            axis.title.y = element_text(size=16),
                            legend.text=element_text(size=16),
                            legend.title = element_text(size=16),
                            text=element_text(size=16,  family="Palatino Linotype"))




