##
##
##
##
##
## create visualizations for layers chapter in thesis
## ------------------------------------------------------

library(scales)
library(FSA)
library(astsa)
library(fields)

## ---------------------------------------------------------------------------------
## INITIAL RESULTS
## ---------------------------------------------------------------------------------


## -------------------------------
## TS boxplot for 0 degree coupons
## --------------------------------
TSZlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZlayers.rds")
TSZnolayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZnolayers.rds")
TSZreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZreal.rds")

TSZreal <- c(TSZreal, rep(NA, 100-16))

zerosAll <- cbind(TSZnolayers, TSZreal, TSZlayers)
colnames(zerosAll) <- c("no layers", "experimental", "layered")

windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)
boxplot(zerosAll, 
        ylab = "test statistic", xlab = "coupon populations", 
        cex.lab = 1.2, family = "A",
        pch = 20, boxlwd = 2, boxwex = 0.5,
        ylim = c(0,1))
mtext("Distribution of Test Statistic, 0 degree Coupons", 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")


## --------------------------------
## TS boxplot for 45 degree coupons
## --------------------------------
TSFlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFlayers.rds")
TSFnolayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFnolayers.rds")
TSFreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/45degreeData/TSFreal.rds")

TSFreal <- c(TSFreal, rep(NA, 100-24))

ffAll <- cbind(TSFnolayers, TSFreal, TSFlayers)
colnames(ffAll) <- c("no layers", "experimental", "layered")

windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)
boxplot(ffAll, 
        ylab = "test statistic", xlab = "coupon populations", 
        cex.lab = 1.2, family = "A",
        pch = 20, boxlwd = 2, boxwex = 0.5,
        ylim = c(0,1))
mtext("Distribution of Test Statistic, 45 degree Coupons", 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")




## ---------------------------------------------------------------------------------
## NOISE (SENSITIVITY) ANALYSIS
## ---------------------------------------------------------------------------------

## -------------------------------
## TS noise for 0 degree coupons
## --------------------------------

TSZlayers <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZlayers.rds")
TSZnoise <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/noisySpacing40to60/TSZnoise.rds")
TSZreal <- readRDS("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/TSZreal.rds")


allNoise <- cbind(TSZlayers, TSZnoise, c(TSZreal, rep(NA, 100-16)))

boxplot(allNoise)


## ---------------------------------------------------------------------------------
## CREATE EXAMPLE PERIODOGRAMS
## ---------------------------------------------------------------------------------
setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/spacing40to60")

inputFiles <- list.files(full.names = TRUE)

which(inputFiles=="./54.141spacingSyn0.rda")

for(coupon in inputFiles[99]){
  
  load(coupon)
  
  h <- hist(~nlsCoupon[,3], w=5, plot = FALSE) 
  
  counts <- h$counts
  coords <- h$mids
  
  
  # detrend the series using a spline 
  xGrid <- seq(min(coords), max(coords), length.out = length(coords))
  highpass <- splint(coords, counts, xGrid, lambda = 50) #lambda set to be very smooth
  
  detrendedData <- counts - highpass
  
  
  # calculate the dft and compute raw periodogram using mvspec(.)
  windowsFonts(A = windowsFont("Times New Roman"))
  par(family = "A", cex.lab = 1.2)
  P <- mvspec(detrendedData, log='no', main = "", col = "black", sub="", ylim = c(0,51.5))
  mtext("Periodogram for Layered Coupon", 
        side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")
  # size: 679x522
  
}



setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/noLayers")

inputFiles <- list.files(full.names = TRUE)

for(coupon in inputFiles[8]){
  
  load(coupon)
  
  h <- hist(~nlsCoupon[,3], w=5, plot = FALSE) 
  
  counts <- h$counts
  coords <- h$mids
  
  
  # detrend the series using a spline 
  xGrid <- seq(min(coords), max(coords), length.out = length(coords))
  highpass <- splint(coords, counts, xGrid, lambda = 50) #lambda set to be very smooth
  
  detrendedData <- counts - highpass
  
  
  # calculate the dft and compute raw periodogram using mvspec(.)
  windowsFonts(A = windowsFont("Times New Roman"))
  par(family = "A", cex.lab = 1.2)
  P <- mvspec(detrendedData, log='no', main = "", col = "black", sub="", ylim = c(0,51.5))
  mtext("Periodogram for Non-Layered Coupon", 
        side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")
  
}

## ---------------------------------------------------------------------------------
## ILLUSTRATE THE PEAK DETECTION PROCESS
## ---------------------------------------------------------------------------------

setwd("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/spacing40to60")

inputFiles <- list.files(full.names = TRUE)

for(coupon in inputFiles){
  
  load(coupon)
  
  load("C:/Users/barna/Documents/Coupons/layers/layerData/0degreeData/spacing40to60/56.162spacingSyn0.rda")

#create series data (count values for hist bins)
h <- hist(~nlsCoupon[,3], w=5, plot = FALSE) 

counts <- h$counts
coords <- h$mids


# detrend the series using a spline 
xGrid <- seq(min(coords), max(coords), length.out = length(coords))
highpass <- splint(coords, counts, xGrid, lambda = 50) #lambda set to be very smooth

detrendedData <- counts - highpass


# calculate the dft and compute smoothed periodogram using mvspec(.)
P <- mvspec(detrendedData, kernel('daniell',2), log='no', main = "", col = "grey40")
mtext("Smoothed Periodogram with Peak Detection", 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")

#P <- mvspec(detrendedData)

# calculate lower bounds of 95% confidence intervals for each
# periodogram ordinate
df = P$df
L = qchisq(.975, df) 
U = qchisq(.025, df)
n = length(P$freq)

lowerBound <- rep(NA, n)
upperBound <- rep(NA, n)
for(i in 1:n){
  lowerBound[i] <- df*P$spec[i]/L
  upperBound[i] <- df*P$spec[i]/U
}

# points(P$freq, upperBound, pch = 20, col = "tomato")
# lines(P$freq, lowerBound, col = "cornflowerblue", lwd = 2)

points(P$freq, lowerBound, pch = 20, col = "cornflowerblue")

# establish baseline for periodogram using linear model
baseline <- lm(P$spec ~ P$freq)

lines(P$freq, baseline$fitted.values, col = "tomato", lwd = 2)

N = length(P$freq)

# get standard errors at each fitted value
pred <- predict(baseline, se.fit = TRUE)

SE <- pred$se.fit
zB <- qnorm(.025/N, lower.tail = FALSE)

xg <- c(P$freq,P$freq[N:1]) 
yg <- c(baseline$fitted.values - zB*SE,(baseline$fitted.values + zB*SE)[N:1])


polygon( xg, yg, border="grey",col=alpha("magenta",.1))
polygon( xg, yg, border="grey",col=alpha("tomato",.3))


legend("topright", legend = "baseline",
       col = c("tomato"), lty = c(1),lwd = 3,inset = 0.05)

# find peaks in periodogram
peaks <- which(lowerBound > baseline$fitted.values)


# calculate test statistic
print(sum(P$spec[peaks])/sum(P$spec))

}
