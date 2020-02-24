##
##
##
##
##
## uncertainty in positioning
##
## construct populations of unit vectors for both pore
## and surface data sets and associated graphics
## --------------------------------------------------------
library(ggplot2)
library(fields)
library(bootstrap)

library(extrafont)
loadfonts(device = "win")

setwd("C:/Users/barna/Documents/Coupons/datasets")
nPores <- readRDS("0and45degreePoreNumbers.rds")

unitVecPore <- matrix(ncol = length(nPores), nrow = 1000, NA)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData/porosity")

i=1
for(n in names(nPores)){
 
  load(paste0("coupon", n, "poreCI.rda"))
  
  load(paste0("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped/nlsCoupon", n, ".rda"))

  
  ## -----------
  
  print(n)
  
  for(j in 1:1000){
   unitVecPore[j,i] <- sqrt( sum( (c(simNlsCoef[j,3:4],1) - c(nlsCoeff[3:4],1))^2 ) )
  }

   i = i + 1
  
}


colnames(unitVecPore) <- names(nPores)


## ----------

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData/surfaces")

unitVecSurface <- matrix(ncol = length(nPores), nrow = 1000, NA)

i=1
for(n in names(nPores)){
  
  load(paste0("coupon", n, "surfaceCI.rda"))
  
  load(paste0("C:/Users/barna/Documents/Coupons/nlsAxis/surfaces/nlsSurfaceData/nlsSurfaceCoupon", n, ".rda"))
  
  
  ## -----------
  
  print(n)
  
  for(j in 1:1000){
    unitVecSurface[j,i] <- sqrt( sum( (c(simNlsCoef[j,3:4],1) - c(nlsCoeff[3:4],1))^2 ) )
  }
  
  i = i + 1
  
}

rm(bigCoupon, nlsCoeff, nlsCoeffBig, nlsCoupon, oldCoupon)

colnames(unitVecSurface) <- names(nPores)



bplot(unitVecSurface)

part1 <- cbind(as.vector(unitVecPore), rep(1,length(as.vector(unitVecPore))))
part2 <- cbind(as.vector(unitVecSurface), rep(2,length(as.vector(unitVecSurface))))

part2 <- part2[-which(part2[,1] > 0.3),]
      
dat <- data.frame(rbind(part1, part2))

dat$X2 <- as.factor(dat$X2)
      


ggplot(data = dat, aes(X1, fill = X2)) +
  ggtitle(label = "Distance From Original Center Axis") +
  labs(x = "distance (2-norm)", y = "density") +
  geom_density(alpha=0.5, kernel="gaussian", adjust=1) +
  theme_bw()+ theme(legend.justification = c(1,1), 
                    legend.position = c(.95, .95), 
                    legend.background = element_blank(),
                    legend.box.background = element_rect(colour = "black"),
                    #legend.key = element_rect(size = 5),
                    #legend.key.size = unit(1.5, 'lines'),
                    plot.title = element_text(size = 20),
                    axis.text=element_text(size=16),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    text=element_text(size=16,  family="Palatino Linotype"),
                    legend.title = element_blank())+
  scale_color_manual(labels = c("porosity", "surface"), 
                     values = c("royalblue", 'red'), 
                     aesthetics = c("colour", "fill"))

+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
           linetype=guide_legend(keywidth = 3, keyheight = 1),
           colour=guide_legend(keywidth = 3, keyheight = 1))




packCI <- bcanon(as.vector(unitVecPore), 1000, median)
median(as.vector(unitVecPore))
surfCI <- bcanon(as.vector(unitVecSurface), 1000, median)
median(as.vector(unitVecSurface))
  
# ## ----------
# ## which coupon is the most uncertain?
# ## ----------
# 
# which.min(colMeans(unitVec))
# 
# bigcolmean <- colMeans(unitVec)[rev(order(colMeans(unitVec)))][1:3]
# 
# 
# ## get the number of pores each coupon has
# numPre <- vector()
# j=1
# for(i in colnames(unitVec)){
#        load(paste0("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped/nlsCoupon", i, ".rda"))
# 
#        numPre[j] <- length(nlsCoupon[,1])
#        j=j+1
# }
# 
# names(numPre) <- colnames(unitVec)
# 
# lowNumpre <- numPre[order(numPre)][c(1, 15, 3)]
# 
# par(family = "Palatino Linotype", cex.lab = 1.2)
# plot(colMeans(unitVec), (numPre), 
#      ylab = "number of pores",
#      xlab = "distance (2-norm)",
#      main = "",
#      pch = 16, frame = T)
# mtext("Number of Pores per Coupon Against Average Distance", 
#       side=3, adj=0, line=0.6, cex=1.6, font=2, family = "Palatino Linotype")
# text(bigcolmean, lowNumpre, 
#      labels = c("coupon 56", "coupon 17", "coupon 31"), pos = c(2,3,3), cex = 0.8)
# 
# i=17
# load(paste0("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped/nlsCoupon", i, ".rda"))
# 
# plot3d(nlsCoupon[,1], nlsCoupon[,2], nlsCoupon[,3], size = 0.45, type = "s")
# 
# boxplot(unitVec[,which(colnames(unitVec) == 17)])
# 
# ## ---------------
# ## GRAPHICS - most uncertain coupon
# ## ---------------
# 
# par(family = "Palatino Linotype", cex.lab = 1.2)
# b<-boxplot(colMeans(unitVec), horizontal = T, boxwex = 0.5,
#         frame = F, xlab = "distance (2-norm)", pch = 20)
# mtext("Distance from Original Axis, Coupon Averages", 
#       side=3, adj=0, line=0.6, cex=1.6, font=2, family = "Palatino Linotype")
# 
# text(b$out, b$group, labels = c("coupon 17", "coupon 31", "coupon 56"), pos = c(1,3,3), cex = 0.8)
# 
# 
# 
# ## -------------
# ## get confidence interval for distance
# ## -------------
# 
# longUnitVec <- as.vector(unitVec)
# 
# packCI <- bcanon(longUnitVec, 1000, median)
# 
# numDelta <- length(longUnitVec)
# loT <- mean(longUnitVec) - qnorm(.975)*sd(longUnitVec)/sqrt(numDelta) 
# upT <- mean(longUnitVec) + qnorm(.975)*sd(longUnitVec)/sqrt(numDelta) 
# 
# 
# med <- mean(longUnitVec)
# 
# dat <- data.frame(cbind(as.vector(unitVec), 
#                         rep(1, numDelta), 
#                         rep(medCI, ( numDelta / 2) ), 
#                         rep(med, numDelta) ) )
#                   
#                   
# names(dat) <- c("deviance", "population", "ci", "med")
# 
# dat$population <- as.factor(dat$population)
# 
# library(extrafont)
# font_import()
# loadfonts(device = "win")
# 
# 
# p <- ggplot(data = dat, aes(deviance)) +
#   ggtitle(label = "Distance From Original Center Axis") +
#   labs(x = "distance (2-norm)", y = "density") +
#   geom_vline(aes(xintercept=med,  color = "med"),
#              lwd=1.3, linetype = "longdash") +
#   geom_density(fill = "cornflowerblue", color = "cornflowerblue", 
#                alpha=0.2, kernel="gaussian", adjust=1) +
#   scale_color_manual(labels = "mean", values = c("tomato")) +
#   #scale_linetype_manual(values=c(3)) +
#   theme_bw()+ theme(legend.justification = c(1,1), 
#                     legend.position = c(.95, .95), 
#                     legend.background = element_blank(),
#                     legend.box.background = element_rect(colour = "black"),
#                     legend.title=element_blank(),
#                     plot.title = element_text(size = 20),
#                     axis.text=element_text(size=12),
#                     axis.title.x = element_text(size=16),
#                     axis.title.y = element_text(size=16),
#                     legend.text=element_text(size=16),
#                     text=element_text(size=20,  family="Palatino Linotype")) +
#   guides(fill = guide_legend(keywidth = 1, keyheight = 1),
#            linetype=guide_legend(keywidth = 3, keyheight = 1),
#            colour=guide_legend(keywidth = 3, keyheight = 1))
# 
# p 
# 
# 
# 
