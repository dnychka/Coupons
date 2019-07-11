##
## 
## function to identify spikes in periodogram 
##
## takes a weighted moving average to measure deviation from 
## a baseline signal strength. Deviation is set as some threshold
## of standard deviations above the moving mean
##
## y - periodogram signal
## bandwidth - width of moving average
## threshold - number of standard deviations
## w - weight which to give positive signals (peaks)
##
## ---------------------------------------------------------------------------

findFreq <- function(y,bandwidth,threshold,w) {

n = length(y)
signal <- rep(0, n)
baseline <- y[1:bandwidth]

meanBaseline <- rep(NA,n)
sdBaseline <- rep(NA,n)

meanBaseline[bandwidth] <- mean(baseline)
sdBaseline[bandwidth] <- sd(baseline)

for( i in (bandwidth+1):n ){
  if( (y[i]-meanBaseline[i-1])  >  threshold*sdBaseline[i-1]){
          signal[i] = 1
          baseline[i] = w*y[i] + (1-w)*baseline[i-1]
  } else {
    signal[i] = 0
    baseline[i] = y[i]
  }
  
  meanBaseline[i] = mean(baseline[(i-bandwidth):i])
  sdBaseline[i] = sd(baseline[(i-bandwidth):i])
}

return(list("signals"=signal,"baselineAvg"=meanBaseline,"baselineStd"=sdBaseline))

}




