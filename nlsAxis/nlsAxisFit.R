##
## function to fit the model by nonlinear least squares
##
## relies on functions 'getInitialParameters' and 'radiusAligned'
##
## poreCoordinates - a tall matrix containing x, y, and z coordinates of pores
##                   (should already be cropped)   
##
## returns nls object
##
## ---------------------------------------------------------------------------


nlsAxisFit <- function(poreCoordinates){  
  
  ## initial radius guess
  r = 1000 # based on ideal coupon radius of 1000 micor-meters
  ## tolerance for radius diff
  tol = 0.001
  diffMedians = 100
  
  
  startValues <- getInitialParameters(poreCoordinates)
  
  N <- length(poreCoordinates[,1])
  
  
  while (diffMedians > tol){
    radiusTarget <- rep(r, N) 
    
    centerAxis <- nls(radiusTarget~radiusAligned(poreCoordinates, 
                                                 centroidX, centroidY, 
                                                 axisVectorX, axisVectorY),
                      start = list(centroidX = startValues[2,1], 
                                   centroidY = startValues[2,2],
                                   axisVectorX = startValues[1,1], 
                                   axisVectorY = startValues[1,2]))
    
    
    nlsCoeff <- coef(centerAxis)
    
    ## run 'getRadius' with the optimal parameter values found by nls
    radiusFinal <- radiusAligned(poreCoordinates, 
                                 nlsCoeff["centroidX"], nlsCoeff["centroidY"], 
                                 nlsCoeff["axisVectorX"], nlsCoeff["axisVectorY"])
    
    
    # update radius and compare to last iteration
    diffMedians <- abs(r - median(radiusFinal))
    r <- median(radiusFinal)

  } # end of while loop
  
  
  return(centerAxis)
  
  
}





