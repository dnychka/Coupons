##
##
##
##
##
## create comparative visualizations for 0 degree coupons
## ------------------------------------------------------




# what do we have saved as of 9/26/2019

zerosLayered <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosLayeredSignals.rds")
zerosNotLayered <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosNullSignals.rds")
zerosReal <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosRealSignals.rds")


len <- length(zerosLayered)

#pad with zeros
zerosNotLayered <- c(zerosNotLayered, rep(NA, 50))
zerosReal <- c(zerosReal, rep(NA,(150-16)))

zerosAll <- cbind(zerosNotLayered, zerosReal, zerosLayered)
colnames(zerosAll) <- c("no layers", "experimental", "layered")

windowsFonts(A = windowsFont("Times New Roman"))
par(family = "A", cex.lab = 1.2)
boxplot(zerosAll, 
        ylab = "test statistic", xlab = "coupon populations", 
        cex.lab = 1.2, family = "A",
        pch = 20, boxlwd = 2, boxwex = 0.5,
        ylim = c(0.1,1))
mtext("Distribution of Test Statistic, 0 degree Coupons", 
      side=3, adj=0, line=1.6, cex=1.4, font=2, family = "A")


