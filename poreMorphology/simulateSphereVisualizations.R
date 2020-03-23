##
##
##
##
##
##
## create visualizations that show sphericity accuracy of calculated
## perfect spheres -- voxel approximation
## -----------------------------------------------------------------

library(ggplot2)
library(car)
library(scatterplot3d)
library(fields)
library(fda)
library(scales)
library(dplyr)


setwd("C:/Users/barna/Documents/Coupons/poreMorphology/morphologyData")

source('naturalCubicSplineR.R')
source('naturalSplineBasis.R')

load("perfectSpheresData.rda")

perfectSpheres <- readRDS("perfectSphereDataRadius3-52.rds")

setwd("C:/Users/barna/Documents/Coupons/poreMorphology/morphologyData/spheres")

attach(perfectSpheres)

ggplot(perfectSpheres[which(perfectSpheres$sphVolume<10000),], aes(x=sphVolume, y=sph)) +
  geom_point(col = "black", size = 2) +
  geom_hline(aes(yintercept= 1, linetype = "ideal sphericity"), color = "cornflowerblue", size = 1.5, show.legend = T) +
  geom_vline(aes(xintercept = 1000))+
  ylim(0,1.2) + ggtitle("sphericity of perfect spheres using voxel approximation")+
  labs(x = "volume", y = "sphericity") +
  theme(plot.title = element_text(size = 20),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title = element_text(size=16),
        legend.position = "bottom")+
  scale_linetype_manual(name = "", values = c(1,1))





## -----------------------------
## put cubic smoothing spline on

perfectSpheres <- perfectSpheres[order(as.numeric(perfectSpheres$poreID)),]

perfectSpheres$equivRadius <- seq(3,52,by=0.25)[-197]

attach(perfectSpheres)

sGrid<-seq( 0,max(equivRadius),length.out=500) # grid for prediction locations
objGCV<-Tps( equivRadius, sph) # omitting smoothing parameter or df means
# Tps(.) finds values by GCV
fHatGCV<-predict( objGCV, sGrid)

plot(scoutscanSphericity$equivRad, scoutscanSphericity$sph, pch = 20, col = alpha("violetred1", .5), ylim = c(0,1))
points(reconSphericity$equivRad, reconSphericity$sph, pch = 20, col = alpha("darkorange", .5))
points(equivRadius, sph,pch=20)
lines( sGrid, fHatGCV, col=alpha("mediumturquoise",1), lwd=3)

# first get the covariance for the predicted values.
# This involves creating a basis matrix for cubic smoothing splines
# to calculate the cubicSpline hat matrix
Phi <-naturalSplineBasis(sGrid,sGrid)
R <-naturalCubicSplineR(sGrid)
lambda <- objGCV$lambda
H <- Phi%*%solve(t(Phi)%*%Phi+lambda*R)%*%t(Phi)
V <- objGCV$shat.GCV^2 * H %*% t(H) # covariance matrix for predicted values
# then standard error for fHat
fHatSE <- sqrt(diag(V))
# create upper and lower confidence limits
zB<-qnorm( .025/length( sGrid), lower.tail = FALSE)
upper<- fHatGCV+zB*fHatSE
lower<- fHatGCV-zB*fHatSE
# polygon band with semitranparent color.
xg<-c( sGrid,rev( sGrid))
yg<-c( lower,rev( upper))

polygon( xg, yg, border="grey",col=alpha("cyan",.2) )

yline(1, lty = 2, col = "darkorange", lwd=2)

xline(6, lty=3)


datapoly <- data.frame(x = xg, y = yg)

dataSpline <- data.frame(yspline = fHatGCV, xspline = sGrid)


ggplot(perfectSpheres, aes(x=equivRadius, y=sph)) +
  geom_segment(aes(x=0,xend=max(equivRadius),y=1,yend=1, linetype = "ideal sphericity", 
               color = "ideal sphericity"), size = 1.5, show.legend = T)+
  geom_polygon(data = datapoly, aes(x = x, y = y), fill = "cyan", alpha = .2)+
  geom_line(data = dataSpline, aes(x = xspline, y = yspline, linetype = "voxel approximated sphericity", 
            color = "voxel approximated sphericity"), lwd = 2) +
  geom_point(col = "grey50", size = 2) +
  ggtitle("sphericity of perfect spheres using voxel approximation")+
  labs(x = "sphere diameter (voxels)", y = "sphericity") + 
  theme(plot.title = element_text(size = 20),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title = element_text(size=16),
        legend.position = "bottom")+
  scale_linetype_manual(name = "", values = c(1,1)) + 
  scale_color_manual(name = "", values = c("cornflowerblue", "mediumturquoise"))
  












## -----------------------------
## which pores > 1000 vol have the largest sph?

vPxls <- vPxls %>% mutate_at("z",rescale)

sphericity$vol <- as.vector(vInd)

namesLessThan1000 <- as.numeric(names(vInd[which(vInd < 1000 & vInd > 10)]))

littlePoresSphericity <- sphericity[which(sphericity$poreID %in% namesLessThan1000),]

themLittle <- which(perfectSpheres$sphVolume < 1500)

ps <- data.frame(xspline(perfectSpheres[themLittle,c("sphVolume", "sph")], shape=1, lwd=2, draw=F))

psF <- approxfun(x=ps$x, y=ps$y)

above <- which(littlePoresSphericity$sph > psF(littlePoresSphericity$vol))

ggplot(littlePoresSphericity, aes(vol, sph)) +
  geom_point(col = "grey70") +
  geom_point(data = littlePoresSphericity[above,], col = "darkorange")+
  geom_path(data=ps, aes(x, y), col="black", size = 1.4)+
  ylim(0,1) +
  labs(x = "volume", y = "sphericity") +
  ggtitle("sphericity for scoutscan pores") +
  theme(plot.title = element_text(size = 20),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title = element_text(size=16))




namesAbove <- littlePoresSphericity$poreID[above]

shapesAbove <- vPxls[which(vPxls$poreID %in% namesAbove),]

with(shapesAbove, scatter3d(x,y,z, surface = F, point.col = "black",
                            axis.scales = F, axis.ticks = F))

save(shapesAbove, file = "shapesAboveShiny.rda")


library(rsconnect)
rsconnect::deployApp('C:/Users/barna/Documents/Coupons/poreMorphology')


## -----------------------------------
## plot a perfect sphere alongside some of these for fun
## -----------------------------------


# sphere/pore volume of about 500

t <- vInd[order(vInd)]

ssCompare <- as.numeric(names(t[which(t > 400 & t < 600)]))

getThese <- which(vPxls$poreID %in% ssCompare)


vPxls <- vPxls %>% mutate_at("z",rescale)

smallSetVPxls <- vPxls[getThese,]


smallSetVPxlsPlusSphere %>% group_by(poreID) %>% summarize(hoeMany = length(unique(z)),
                                                           lowZ = range(z)[2],
                                                           upZ = range(z)[1],
                                                           lowZup = range(z)[2]-range(z)[1],
                                                           zDenisty = hoeMany / lowZup)




myrescale <- function(x){x / 100}

scaledSphere <- vol500 %>% mutate_at(c("x","y","z"), myrescale)

mydensity <- function(x){rescale(x, to = c(0.628, 0.608))}

scaledSphere <- scaledSphere %>% mutate_at(c("x","y","z"), mydensity)

with(scaledSphere, scatter3d(x, y, z, surface = F, point.col = "black"))

smallSetVPxlsPlusSphere <- bind_rows(smallSetVPxls, scaledSphere)

with(smallSetVPxlsPlusSphere, scatter3d(x, y, z, surface = F, 
                                        point.col = color.scale(sph, zlim = c(0,1)),
                                        axis.col = "white"))




## --------------------------------------------
## get a picture of the simulated spheres
## --------------------------------------------

rad <- c("3", "5","10","13","25")
spacing <- seq(0,150, length.out = 5)

spacing <- c(0, 7, 20)

perfMat <- matrix(NA, nrow = 1, ncol = 4)

for(i in 1:3){
  
  tempVol <- readRDS(paste0("r",rad[i],"pixels.rds"))
  tempVol <- tempVol + spacing[i]
  
  tempVol <- cbind(tempVol, perfectSpheres$sph[which(perfectSpheres$poreID == as.numeric(rad[i]))])
  
  perfMat <- rbind(perfMat, tempVol)
  
}

perfMat <- perfMat[-1,]

perfMat <- as_tibble(perfMat)

with(perfMat, scatter3d(X, Y, Z, surface = F, point.col = "black"))


with(perfMat, scatterplot3d(X, Y, Z, box = F,
                            pch = 15, highlight.3d = T))



testMe <- as_tibble(readRDS("r10pixels.rds"))

with(testMe, 
     scatterplot3d(X, Y, Z, angle = 10, box = F, 
                   pch = 15, highlight.3d = T,
                   main = "",
                   cex.main = 2,
                   xlab = "", ylab = "", zlab = "",
                   tick.marks = F))
text(x = 0, y = 60, "Y-axis")


irregularPore <- sPxls[which(sPxls$poreID==1),]

subIrregPore <- sample_frac(irregularPore, 3/5)

with(subIrregPore, 
     scatterplot3d(x, y, z,
                   angle = 10, box = F, 
                   pch = 15, highlight.3d = T,
                   main = "",
                   cex.main = 2,
                   xlab = "", ylab = "", zlab = "",
                   tick.marks = F))


sphericity[which(sphericity$sa > 1000 & sphericity$sa < 1500),]
     