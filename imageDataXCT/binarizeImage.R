##
## 5 June 2019
##
##
## use imageMagick to get a binary .tif file
## ------------------------------------------------------------------------------


library(magick)

setwd("C:/Users/barna/Documents/Coupons/imageDataXCT")

## 
## get a binary .tif file for testing. This step unecessary for actual coupon data
## -------------------------------------------------------------------------------
fishMagick <- image_read("C:/Users/barna/Documents/Coupons/imageDataXCT/fish.tif")

plot(fishMagick)

fishGrey <- image_convert(fishMagick, colorspace = "Gray")

plot(fishGrey)

fishBW <- fishGrey %>%
  image_threshold(type = "white", threshold = "50%") %>%
  image_threshold(type = "black", threshold = "50%")

plot(fishBW)

image_write(fishBW, 'fishBW.tif', format = 'tif')
## -------------------------------------------------------------------------------









