##
##
##
##
##
##
##
##
##
## ---------------------------------------------------

library(raster)
library(fields)


testImg <- raster("D:/NADCA/processed/D0276.tif")

plot(testImg)

ti2 <- crop(testImg, extent(testImg, rowFromY(testImg, 4000), rowFromY(testImg, 3700), 200, 900))

plot(ti2)

contour(ti2, add=T)

xline(c(200, 900))
yline(c(3000,5000))

rowFromY(testImg, 5000)
