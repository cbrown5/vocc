# Compare new VoCC calculations to old ones
# CJ Brown 18 Apr 2017
# TODO:
# There are some differences in cul-de-sacs, due to gradients being slightly different.
# My new code tends to calculate much greater values of vocc in these places.
#Overall this won't affect maps much, but would be good to know why difference occurs.

rm(list = ls())

library(raster)
library(reshape)
library(CircStats)

library(RColorBrewer)

library(dplyr)
library(tidyr)

rbase <- raster(nrow = 180, ncol = 90)
rbase[] <- NA

devtools::load_all("~/Code/geoengineering/vocc")

source('~/Code/geoengineering/vocc/data-raw/Velocity_functions.R')

# ---------------
# Run scripts
# ---------------
# Functions appended '2' are the old ones  we are checking against
data(sst)

system.time(slopedat <- calcslope(sst))

allyears <- rep(1, nlayers(sst))
mnsst <- stackApply(sst, indices = allyears, fun = mean) # Calculate mean SST over all years

system.time(spatx <- spatialgrad(mnsst))
system.time(spatx2 <- spatialgrad2(mnsst))

system.time(velodf <- calcvelocity(spatx, slopedat))
system.time(velodf2 <- calcvelocity2(spatx2, slopedat))

# ---------------
# Plots to check results
# ---------------

r <- raster(sst)

#Old code - already checked trend, pretty confident I have it right.
spvel <- velodf2
coordinates(spvel) <- ~ x + y
icell <- cellFromXY(r, spvel)
rgrad2 <- r
rgrad2[icell] <- spvel$NSgrad
rvocc2 <- r
rvocc2[icell] <- spvel$velocity

#new code
rgrad <- r
rgrad[spatx$icell] <- spatx$NS
rvocc <- r
rvocc[velodf$icell] <- velodf$velocity

#Check
rdiffg <- rgrad -rgrad2
plot(rdiffg, col = rev(brewer.pal(9, "Reds")))
sum(rdiffg[] != 0, na.rm = T)

rdiffv <- rvocc - rvocc2
plot(rdiffv, col = rev(brewer.pal(9, "Reds")))
zoom(rdiffv, col = rev(brewer.pal(9, "Reds")))
sum(rdiffv[] != 0, na.rm = T)

#Check if data is lost
sum(!is.na(rvocc[]))
sum(!is.na(rvocc2[]))
# ---------------
# Attempt to speed up Dave's code for slopes, but twice as slow!
# ---------------



getslope <- function(x, y){
	desmat <- as.matrix(cbind(rep(1, length(x)), x))
	mod <- .lm.fit(desmat, y)
	mod$coef[2]
}

mnthmean <- sst
calcslope2 <- function(mnthmean, times = NULL){
	if (is.null(times)) times <- data.frame(times = 1:nlayers(mnthmean), colnames = names(mnthmean))
	dat <- data.frame(mnthmean[]) %>%
		mutate(cellnum = 1:ncell(mnthmean)) %>%
		gather(key = colnames, value = temp, -cellnum) %>%
		left_join(times) %>%
		mutate(int = 1) %>%
		filter(!is.na(temp)) %>%
		group_by(cellnum) %>%
		summarize(slope = getslope(times, temp))

	return(dat)
}

system.time(slopedat <- calcslope(sst))
system.time(slopedat2 <- calcslope2(sst))


slopedat
