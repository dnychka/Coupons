##
##
##
##
##
##
##
## compute the acceleration of each nls parameter
## by the jackknife method. 
## ---------------------------------------------------------

getAcceleration <- function(poreCoordinates){

  numPores <- length(poreCoordinates[,1])

  jackknifeValues <- matrix(nrow = numPores, ncol = 5, NA)

  for(i in 1:numPores){
  
  
    # remove one pore at a time
    jackCoords <- poreCoordinates[-i,]
    
    # estimate the parameter values for each new jackknife data set
    N <- length(jackCoords[,1])
    
    startValues <- getInitialParameters(jackCoords)
    
    
    distTarget <- rep(0,N)
    
    nlsObj <- try(nls(distTarget~getDistance(jackCoords, 
                                             centroidX, centroidY, 
                                             axisVectorX, axisVectorY, r),
                      start = list(centroidX = startValues[4], 
                                   centroidY = startValues[5],
                                   axisVectorX = startValues[1], 
                                   axisVectorY = startValues[2],
                                   r = startValues[7]),
                      control =  nls.control(minFactor = 1/10000000000)))
  
    
    jackknifeValues[i,] <- coef(nlsObj)
  
  
  }


  # calculate acceleration
  jDot <- vector()
  jDot <- apply(jackknifeValues, 2, function(x) sum(x/numPores))
  
  a <- vector()
  for(i in 1:5){
    a[i] <- sum( (jDot[i] - jackknifeValues[,i])^3 ) / (6 * ( sum( (jDot[i] - jackknifeValues[,i])^2 ) )^(3/2))
  }
  
  return(a)

}




