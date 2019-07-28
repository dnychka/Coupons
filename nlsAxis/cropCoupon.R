##
## function to crop the coupon to give to function nlsAxisFit
##
## n - index of coupon of interest
## poreData - a tall matrix containing x, y, and z coordinates of pores
##            (taken straight from unwrapped json files)  
##
## returns poreCoordinates, a tall matrix of x, y, and z coords of pores
## from the cropped coupon
##
## ---------------------------------------------------------------------------

cropCoupon <- function(n, poreData){
  ordered <- order(poreData[[n]]$comZ)
  
  comX <- poreData[[n]]$comX[ordered]
  comY <- poreData[[n]]$comY[ordered]
  comZ <- poreData[[n]]$comZ[ordered]
  
  oldCoupon <- cbind( comX,
                      comY,
                      comZ)
  
  cropSections <- quantile(comZ, prob = seq(0, 1, length = 11), type = 5)
  
  ## subset the coupon to avoid the weld/support material remnants
  good <- (comZ >= cropSections[3] & comZ <= cropSections[9])
  
  poreCoordinates <- cbind( comX[good],
                            comY[good],
                            comZ[good])
  
  return(poreCoordinates)
}