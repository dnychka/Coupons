##
##
##
## comparison graphics
##
## -------------------------------------------------------------

## -------------------------------------------------------------
## test signal strength for different layerings with set
## window of fundFreq < 0.15
## -------------------------------------------------------------

zeroLayeredSignal <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosLayeredSignals.rds")
zeroVeryLayeredSignal <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosFineLayeredSignals.rds")
zeroNullSignal <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosNullSignals.rds")

zeroVeryLayeredSignal <- c(zeroVeryLayeredSignal, rep(NA, (length(zeroLayeredSignal)-length(zeroVeryLayeredSignal))))
zeroNullSignal <- c(zeroNullSignal, rep(NA, (length(zeroLayeredSignal)-length(zeroNullSignal))))


datZ <- cbind(zeroLayeredSignal, zeroVeryLayeredSignal, zeroNullSignal)
boxplot(datZ, main = "signal strength, different layerings, zeros", pch = 20)

## --------------------------

ffLayeredSignal <- readRDS("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/45sLayeredSignals.rds")
ffVeryLayeredSignal <- readRDS("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/45sFineLayeredSignals.rds")
ffVeryLayeredSignal <- apply(ffVeryLayeredSignal, 2, max)
ffNullSignal <- readRDS("C:/Users/barna/Documents/Coupons/layers/45degreeSlice/45degreeData/45sNullSignals.rds")
ffNullSignal <- apply(ffNullSignal[-1,], 2, max)

ffVeryLayeredSignal <- c(ffVeryLayeredSignal, rep(NA, (length(ffLayeredSignal)-length(ffVeryLayeredSignal))))
ffNullSignal <- c(ffNullSignal, rep(NA, (length(ffLayeredSignal)-length(ffNullSignal))))


datFF <- cbind(ffLayeredSignal, ffVeryLayeredSignal, ffNullSignal)
boxplot(datFF, main = "signal strength, different layerings, forty-fives", pch = 20)




## -------------------------------------------------------------
## test signal strength for differing amounts of noise
## in the 0 degree coupons
## -------------------------------------------------------------
noNoise <-  readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosLayeredSignals.rds")
low <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/noisySpacing40to60/noisy01spacing40to60signals.rds")
medLow <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/noisySpacing40to60/noisy02spacing40to60signals.rds")
med <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/noisySpacing40to60/noisy03spacing40to60signals.rds")
medHigh <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/noisySpacing40to60/noisy04spacing40to60signals.rds")
high <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/noisySpacing40to60/noisy05spacing40to60signals.rds")

realC <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosRealSignals.rds")
realC <- c(realC, rep(NA, (length(noNoise)-length(realC))))

nulls <- readRDS("C:/Users/barna/Documents/Coupons/layers/0degreeSlice/0degreeData/zerosNullSignals.rds")


datNoise <- data.frame(noNoise, low, medLow, med, medHigh, high, nulls, realC)

boxplot(datNoise, pch = 20, boxwex = 0.5,
        names = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "nulls", "reals"),
        ylab = "signal strength", xlab = "noise as percentage of layer spacing",
        border = ifelse(names(datNoise)=="realC", "firebrick1", 
                        ifelse(names(datNoise)=="nulls", "cornflowerblue", "black")),
        main = "effect of noise on signal strength, zero degree coupons")

t.test(noNoise,realC)
t.test(low, realC)
t.test(medLow, realC)
t.test(med,realC)
t.test(med, nulls)
t.test(medHigh,realC) #t-test at alpha=0.05 fails here
t.test(medHigh, nulls)
t.test(high,realC)
t.test(high,nulls)
t.test(nulls,realC)

xline(4.5, col = "grey", lty=3)
