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



