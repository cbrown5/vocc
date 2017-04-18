# Compare G3 and 4.5 models
# CJ Brown 18 Apr 2017

rm(list = ls())
library(RColorBrewer)

devtools::load_all("~/Code/geoengineering/vocc")

setwd("~/Databases/geoengineering/hadley_model")
# ---------------
# Load data
# ---------------
sst_rcp <- stack("Hadley_rcp45.grd")
sst_g3 <- stack("Hadley_g3.grd")

yrs <- paste0("X",as.character(2017:(2017 +49)))
sst_rcp2 <- subset(sst_rcp, which(names(sst_rcp) %in% yrs))
sst_g32 <- subset(sst_g3, which(names(sst_g3) %in% yrs))

# ---------------
# Run scripts
# ---------------

vocc_raster <- function(r, maxval = 500){
	allyears <- rep(1, nlayers(r))
	mnr <- stackApply(r, indices = allyears, fun = mean) # Calculate mean SST over all years

	slopedat <- calcslope(r)
	spatx <- spatialgrad(mnr)
	velodf <- calcvelocity(spatx, slopedat)
	rvocc <- raster(r)
	rvocc[velodf$icell] <- velodf$velocity
	rvocc[rvocc[] > maxval] <- maxval
	rvocc[rvocc[] < -maxval] <- -maxval
	return(rvocc)
}


vocc_rcp <- vocc_raster(sst_rcp2, maxval = 500)
vocc_g3 <- vocc_raster(sst_g32, maxval = 500)

# ---------------
# Plots
# ---------------
landcol <- "grey"

dev.new(width = 10, height = 4)
par(mfrow = c(1,2))
plot(vocc_rcp, col = rev(brewer.pal(9, 'RdBu')), colNA = landcol)
plot(vocc_g3, col = rev(brewer.pal(9, 'RdBu')), colNA = landcol)
