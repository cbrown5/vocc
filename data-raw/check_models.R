# Run VoCC models for a test layer
# CJ Brown 18 Apr 2017

rm(list = ls())
library(RColorBrewer)

library(raster)
library(reshape)
library(CircStats)

devtools::load_all("~/Code/geoengineering/vocc")
source('~/Code/geoengineering/vocc/data-raw/Velocity_functions.R')
# ---------------
# Create raster brick
# ---------------
set.seed(42)
nr <- 5
ncell <- nr ^ 2

nyrs <- 20
rx <- NULL
for (i in 1:nyrs){
	x <- raster(matrix(rep(20+(1:nr), each = nr), nrow = nr) + rnorm(ncell, sd = 0.2))
	x[10:12] <- NA
	rx <- c(rx, list(x+i))
}

sst <- stack(rx)

plot(sst)

# ---------------
# Run scripts
# ---------------
allyears <- rep(1, nlayers(sst))
mnsst <- stackApply(sst, indices = allyears, fun = mean) # Calculate mean SST over all years

slopedat <- calcslope(sst)
spatx <- spatialgrad(mnsst)
velodf <- calcvelocity(spatx, slopedat)

slopedat2 <- calcslope2(sst)
spatx2 <- spatialgrad2(mnsst)
velodf2 <- calcvelocity2(spatx2, slopedat2)

# ---------------
# Plots to check results
# ---------------

rvocc <- raster(sst)
rvocc[velodf$icell] <- velodf$velocity

plot(rvocc, col = rev(brewer.pal(9, 'Reds')))

plot(velodf2$velocity, velodf$velocity)
abline(0,1)
