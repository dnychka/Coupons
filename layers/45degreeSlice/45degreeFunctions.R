## 6/25/2019
##
##
##
## functionalize 45degreeSlice
##
## this function rotates the coupon to find proper rotational
## orientation. 
## ------------------------------------------------------


## --------------rotation matrices------------------------
aboutY <- function(phi) {matrix( c(cos(phi), 0, sin(phi), 0, 1, 0, -sin(phi), 0, cos(phi)), 3, 3)}

aboutZ <- function(theta) {matrix( c(cos(theta), sin(theta), 0, -sin(theta), cos(theta), 0, 0, 0, 1), 3, 3 )}
## ------------------------------------------------------


## -------------rotation function------------------------
FFrotate <- function(angleSeq, centerCoupon){

histSave <- list()

i = 1

for(theta in angleSeq){
  
  Rz <- aboutZ(theta)
  
  RzCoupon <- centerCoupon %*% Rz
  
  
  phi <- 45*pi/180 #45 degree tilt
  
  Ry <- aboutY(phi)
  
  RzyCoupon <- RzCoupon %*% Ry
 
  h<-hist(~RzyCoupon[,3], w=5, plot = FALSE)
  
  histSave[[i]] <- rbind(h$counts, h$mids)
  
  i = i + 1
  
}

return(histSave)

}
## ------------------------------------------------------


## -----------get relative peak strengths----------------
makePeriodogram <- function(angleSeq, histSave){
  
  findTheta <- rep(NA,length(angleSeq))
  frDt <- list()
  Dtlen <- rep(NA, length(angleSeq))
  
  relSignal <- rep(NA, length(angleSeq))
  
  threshold <- 10
  influence <- 0
  
  for(m in 1:length(angleSeq)){
    
    xGrid <- seq(min(histSave[[m]][2,]), max(histSave[[m]][2,]), length.out = length(histSave[[m]][1,]))
    
    lowpass <- splint(histSave[[m]][2,], histSave[[m]][1,], xGrid, lambda = 50)
    
    frDt[[m]] <- fft(histSave[[m]][1,]-lowpass)
    
    Dtlen[m] <- length(histSave[[1]][1,])
    
    Fr <- (1:Dtlen[m]/Dtlen[m])[1:(Dtlen[m]/2)]
    P <- (Mod(2*frDt[[m]]/Dtlen[m])^2)[1:(Dtlen[m]/2)]
    
    lag <- length(Fr)/10
    
    peak <- findFreq(P,lag,threshold,influence)
    
    ind <- which(peak$signals==1)
    
    findHarmonic <- vector()
    
    if(length(ind) != 0){# check to make sure there's actually peaks
        
        Ppeaks <- P[ind]
        Frpeaks <- Fr[ind]
        
        origFreq <- Frpeaks[which.max(Ppeaks)]
        
        for(j in c(1,1/2,1/3,1/4,1/5)){
          fundFreq <- origFreq*j
          if(fundFreq <= .15){break}
        }
        
        ##find the harmonics
        for(j in 1:5){
          for(i in 1:length(Fr)){
            
            if(isTRUE(all.equal(fundFreq*j, Fr[i], tol = 0.05))){
              if(Fr[i] %in% findBand){
                NA #don't want to repeat values
              } else {findHarmonic <- c(findHarmonic, P[i])
              findBand <- c(findBand, Fr[i])}
            } else {NA}
            
          }
        }
 
    } else {
      fundFreq <- Fr[which.max(P[which(Fr<=0.15)])]
      ##find the harmonics
      for(j in 1:5){
        for(i in 1:length(Fr)){
          
          if(isTRUE(all.equal(fundFreq*j, Fr[i], tol = 0.05))){
            if(Fr[i] %in% findBand){
              NA #don't want to repeat values
            } else {findHarmonic <- c(findHarmonic, P[i])
            findBand <- c(findBand, Fr[i])}
          } else {NA}
          
        }
      }
    }
    
    relSignal[m] <- sum(findHarmonic)/sum(P)
    
  }
  
  return(relSignal)
  
}
## ------------------------------------------------------




