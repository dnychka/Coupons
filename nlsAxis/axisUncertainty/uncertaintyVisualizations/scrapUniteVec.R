##
##
##
##
##
##
##
## exploratory analysis for the coupon pore uncertainty
## --------------------------------------------------------
library(ggplot2)

setwd("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData")

fileList <- list.files("C:/Users/barna/Documents/Coupons/nlsAxis/axisUncertainty/uncertaintyData")


unitVec <- matrix(ncol = length(fileList), nrow = 1000, NA)
cNum <- vector()

i=1
for(couponCI in fileList){
  
  load(couponCI)
  
  cNum[i] <- gsub("poreCI.rda", "", gsub("coupon", "", couponCI)) # a little workaround to get coupon number
  
  print(cNum[i])
  
  load(paste0("C:/Users/barna/Documents/Coupons/nlsAxis/porosity/nlsPorosityData/cropped/nlsCoupon", cNum[i], ".rda"))

  
  ## -----------
  
  
  for(j in 1:1000){
    unitVec[j,i] <- sqrt( sum( (c(simNlsCoef[j,3:4],1) - c(nlsCoeff[3:4],1))^2 ) )
  }
  
  i = i + 1
  
}

colnames(unitVec) <- cNum

## -------------
## get confidence interval for distance
## -------------

longUnitVec <- as.vector(unitVec)


sampMean <- function(x, i) {
  mean(x[i])
}

bootUnitVec <- boot(longUnitVec, sampMean, R = 1000)


zHat <- qnorm(sum(ifelse(bootUnitVec$t > bootUnitVec$t0, 1, 0))/bootUnitVec$R)

jackknifeValues <- vector()
for(i in 1:bootUnitVec$R){
  jackknifeValues[i] <- mean(bootUnitVec$t[-i])
}

aHat <- sum( (jDot - jackknifeValues)^3 ) / (6 * ( sum( (jDot - jackknifeValues)^2 ) )^(3/2))

alpha = 0.025 #significance level

lo <- pnorm(zHat + (zHat + qnorm(alpha)) / (1 - aHat * (zHat + qnorm(alpha)) ))
high <- pnorm(zHat + (zHat + qnorm(1-alpha)) / (1 - aHat * (zHat + qnorm(1-alpha)) ))


hist(longUnitVec)
probs <- c(0.025, 0.975 )
quantiles <- quantile(as.vector(unitVec), prob=probs)
xline(quantiles, col = "cornflowerblue", lwd=2)
xline(packCI$confpoints[c(1,8)], col = "violetred1", lwd = 2, lty=2)


packCI <- bcanon(longUnitVec, 1000, median)

numDelta <- length(longUnitVec)
loT <- mean(longUnitVec) - qnorm(.975)*sd(longUnitVec)/sqrt(numDelta) 
upT <- mean(longUnitVec) + qnorm(.975)*sd(longUnitVec)/sqrt(numDelta) 

c(loT, upT)
medCI <- c(packCI$confpoints[1,2], packCI$confpoints[8,2])

mean(longUnitVec)

MeanCI(longUnitVec, level = 0.975, method = "classic")


hist(longUnitVec, xlim = c(0.0240, 0.025))
xline(mean(longUnitVec))
xline(c(loT, upT), col="violetred1", lty=2)
xline(c(packCI$confpoints[1,2], packCI$confpoints[8,2]), col = "darkorange")

par(mfrow=c(3,3))
apply(unitVec, 2, function(x) hist(x))

dev.off()

hist(as.vector(unitVec))

setwd("C:/Users/barna/Documents/Coupons/datasets")
load("couponCov.rda")

zeros <- c(which(couponCov$polarAngle==0))

numCnum <- as.numeric(cNum) 

zInd <- which(numCnum %in% zeros == T)

par(mfrow=c(1,2))
hist(as.vector(unitVec[,zInd]))

hist(as.vector(unitVec[,-zInd]))

N <- length(as.vector(unitVec))

probs <- c(0.025, 0.975 )
quantiles <- quantile(as.vector(unitVec), prob=probs)

med <- mean(longUnitVec)

dat <- data.frame(cbind(as.vector(unitVec), rep(1, numDelta), rep(medCI, ( numDelta / 2) ), rep(med, numDelta) ) )
                  
                  
names(dat) <- c("deviance", "population", "ci", "med")

dat$population <- as.factor(dat$population)

library(extrafont)
font_import()
loadfonts(device = "win")


p <- ggplot(data = dat, aes(deviance)) +
  ggtitle(label = "Distance From Original Center Axis") +
  labs(x = "distance (2-norm)", y = "density") +
  geom_vline(aes(xintercept=med,  color = "med"),
             lwd=1.3, linetype = "longdash") +
  geom_density(fill = "cornflowerblue", color = "cornflowerblue", 
               alpha=0.2, kernel="gaussian", adjust=1) +
  scale_color_manual(labels = "mean", values = c("tomato")) +
  #scale_linetype_manual(values=c(3)) +
  theme_bw()+ theme(legend.justification = c(1,1), 
                    legend.position = c(.95, .95), 
                    legend.background = element_blank(),
                    legend.box.background = element_rect(colour = "black"),
                    legend.title=element_blank(),
                    plot.title = element_text(size = 20),
                    axis.text=element_text(size=12),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    legend.text=element_text(size=16),
                    text=element_text(size=20,  family="Palatino Linotype")) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
           linetype=guide_legend(keywidth = 3, keyheight = 1),
           colour=guide_legend(keywidth = 3, keyheight = 1))

p 



