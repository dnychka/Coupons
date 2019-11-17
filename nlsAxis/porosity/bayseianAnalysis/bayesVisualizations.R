##
##
##
##
##
##
##
## Bayseian uncertainty FIGURES
## -------------------------------------------------------

library(fields)
library(MASS)
library(FSA)
library(rgl)
library(scales)


setwd("C:/Users/barna/Documents/Coupons/datasets")
poreData <- readRDS("porosityData.rds")
load("couponCov.rda")

setwd("C:/Users/barna/Documents/Coupons/nlsAxis")
source("initialParameters.R")
source("rotateCoupon.R")
source("getRadius.R")
source("newCoupon.R")
source("cropCoupon.R")
source("nlsAxisFit.R")


##--------------------------------------------------------------------
## RESIDUAL PLOTS, interior pores and no interior pores
## (sample coupon, one at a time)
##--------------------------------------------------------------------

m=58

## generate test coupon to show residual deviance
##------------------------------------------------
poreCoordinates <- cropCoupon(m, poreData)
firstFit <- nlsAxisFit(poreCoordinates)

rseRealTest <- summary(firstFit[[1]])$sigma
firstFitRes <- residuals(firstFit[[1]])


param <- coef(firstFit[[1]])
surfaceCoupon <- newCoupon(poreCoordinates, param[1], param[2], param[3], param[4])

#interiorLow <- median(firstFitRes)-sd(firstFitRes)
interiorHigh <- median(firstFitRes)+sd(firstFitRes)
discard <- c(which(firstFitRes > interiorHigh))

## polar for surface plot
##-----------------------------------------------

X <-surfaceCoupon[,1] - coef(firstFit[[1]])["axisVectorX"]
Y <-surfaceCoupon[,2] - coef(firstFit[[1]])["axisVectorY"]

r <- sqrt(X^2 + Y^2)
theta <- atan2(Y, X)
result <- data_frame(r = r, theta = theta, x = X, y = Y)
result$theta <- result$theta + (result$y < 0) * 2 * pi

greys <- grey.colors(firstFitRes)

ctab <- color.scale(firstFitRes[discard], col=two.colors(start="#ffffb2", end="#bd0026", 
                                                         middle = "#fd8d3c"))


windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)

plot(result$theta, surfaceCoupon[,3], 
     col = alpha(greys,0.8),  pch = 16,
     main = "",
     xlab = "angle (radians)", ylab = "z axis",
     cex=2, axes = F)
points(result$theta[discard], surfaceCoupon[,3][discard],
       col = ctab, pch = 16, cex=2)


axis(1, at=seq(0,6, by=1), labels=seq(0,6, by=1))
axis(2, at=seq(-3300,-900, by=400), labels=seq(0,2400, by=400))

mtext("Surface Plot with Interior Pores Colored by Residual Value", 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")

cVal <- firstFitRes[discard]
vertical.image.legend(zlim = range(cVal), col= ctab[order(ctab)])

##--------------------------------------------------------------------
## attempts at boxplot, sample coupon residuals
##--------------------------------------------------------------------

surfB <- boxplot(firstFitRes, plot=FALSE)

bxplotCtab <- color.scale(surfB$out, col=two.colors(start="#ffffb2", end="#bd0026", 
                                                                       middle = "#fd8d3c"))

boxplot(firstFitRes, outcol = "black", axes = F, boxwex = 0.4)
points(rep(1,length(surfB$out)),surfB$out, bg = bxplotCtab, col = "grey40", pch = 21)

hist(firstFitRes)

##--------------------------------------------------------------------
## dot plot, sample coupon residuals
##--------------------------------------------------------------------


vec <- seq(1,length(firstFitRes))

# greys <- grey.colors(firstFitRes)
# 
# plot(vec,firstFitRes, col = greys, pch = 16)
# points(vec[discard],firstFitRes[discard], col = ctab, pch = 16)

interior <- median(firstFit[[2]])-sd(firstFit[[2]])
discard <- which(firstFit[[2]] < interior)
keep <- which(firstFit[[2]] >= interior)

outerPores <- poreCoordinates[-discard,]
nlsObj <- nlsAxisFit(outerPores)

rseOuterTest <- summary(nlsObj[[1]])$sigma
nlsObjRes <- residuals(nlsObj[[1]])

windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)

plot(firstFitRes, col = "white", pch = 8, frame = F,
     xlab = "pore ID", ylab = "residual value",
     axes = F)
points(keep[discard], firstFitRes[discard], pch = 8, col = "grey60")
points(keep, nlsObjRes, col = "grey30", pch = 16)
mtext("Residual Values for Two Models fit to Coupon 58", 
      side=3, adj=0, line=1.1, cex=1.4, font=2, family = "A")
axis(1, at=seq(0,1300, by=100), labels=seq(0,1300, by=100))
axis(2, at=seq(-200,810, by=100), labels=seq(-200,810, by=100))

legend(x=-75, y = 900, c("model fit with interior pores", "model fit without interior pores"), 
       xpd = T, horiz = TRUE, 
       adj=0, cex=1, bty = "n", 
       pch = c(8,16), col = c("grey60", "grey30"))

dev.off()


##--------------------------------------------------------------------
## box plot, layer uncertainty results 0 degree coupons
##--------------------------------------------------------------------

TSuncertainty <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/bayseianAnalysis/bayesData/porosityUncertainty.rds")
TSZlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZlayers.rds")
TSZreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZreal.rds")
TSZuncertain <- TSuncertainty[1:16,]


uncertainZeros <- cbind(as.vector(TSZuncertain),TSZlayers)
colnames(uncertainZeros) <- c("experimental", "layered")

windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)

boxplot(uncertainZeros, 
        col = alpha("grey60", 0.5), border = alpha("grey60", 0.5),
        ylab = "test statistic", xlab = "coupon populations",
        boxwex=0.7, pch=20,
        ylim = c(0,1))

boxplot(TSZreal, boxwex = 0.5, pch=20, add=TRUE, xlab = "")
boxplot(2, TSZlayers, add=TRUE, border = alpha("black", 1), xlab = "")
mtext(expression(paste("Uncertainty in Layer Detection, 0", degree, " coupons")), 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")


#--------------------------------------------------------------------
## box plot, layer uncertainty results 45 degree coupons
##--------------------------------------------------------------------

TSuncertainty <- readRDS("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/bayseianAnalysis/bayesData/porosityUncertainty.rds")
TSFlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFlayers.rds")
TSFreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFreal.rds")
TSFuncertain <- TSuncertainty[17:40,]


uncertainF <- cbind(as.vector(uncertaintyInSignals),TSFlayers)
colnames(uncertainF) <- c("experimental", "layered")

windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)

boxplot(uncertainF, 
        col = alpha("grey60", 0.5), border = alpha("grey60", 0.5),
        ylab = "test statistic", xlab = "coupon populations",
        boxwex=0.7, pch=20,
        ylim = c(0,1))

boxplot(TSFreal, boxwex = 0.5, pch=20, add=TRUE, xlab = "")
boxplot(2, TSFlayers, add=TRUE, border = alpha("black", 1), xlab = "")
mtext(expression(paste("Uncertainty in Layer Detection, 45", degree, " coupons")), 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")

