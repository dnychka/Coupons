##
##
##
##
##
## create comparative visualizations for 45 degree coupons
## ------------------------------------------------------

## maybe junk
layeredCoupons <- apply(layeredCoupons,2,max)
nullCoupons <- c(apply(nullSignals[-1,], 2, max), rep(NA,
                                                      (length(layeredCoupons)-length(apply(nullSignals[-1,],2,max)))))
realCoupons <- c(apply(realSignals,2,max), 
                 rep(NA, (length(nullCoupons)-length(apply(realSignals,2,max)))))

dat <- cbind(nullCoupons, realCoupons)

dev.off()
boxplot(dat, main = "relative peak strength, 45 degree coupons",
        pch = 16, boxwex = 0.7)


dat <- cbind(dat, layeredCoupons)

boxplot(dat, main = "relative peak strength, 45 degree coupons",
        pch = 16, boxwex = 0.7)

t.test(nullCoupons, realCoupons)
t.test(layeredCoupons, realCoupons)



## what do we have saved as of 9/26/2019

setwd("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData")
nullSignals <- readRDS("45sNullSignals.rds")
layeredCoupons <- readRDS("45sLayeredSignals.rds")
realCoupons <- readRDS("45sRealSignals.rds")

nullSignals <- c(nullSignals, rep(NA, 45))
realCoupons <- c(realCoupons, rep(NA, 145-24))


ffAll <- cbind(nullSignals, realCoupons, layeredCoupons)
colnames(ffAll) <- c("no layers", "experimental", "layered")

windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)
boxplot(ffAll, 
        ylab = "test statistic", xlab = "coupon populations", 
        cex.lab = 1.2, family = "A",
        pch = 20, boxlwd = 2, boxwex = 0.5,
        ylim = c(0.1,1))
mtext("Distribution of Test Statistic, 45 degree Coupons", 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")

