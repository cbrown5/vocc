## ----fig.width = 10, fig.height = 7, warning=FALSE, message = FALSE------

library(vocc)
data(sst)
plot(sst, 1:4)

## ------------------------------------------------------------------------
slopedat <- calcslope(sst)

## ------------------------------------------------------------------------
allyears <- rep(1, nlayers(sst))
mnsst <- stackApply(sst, indices = allyears, fun = mean)
spatx <- spatialgrad(mnsst)

## ------------------------------------------------------------------------
velodf <- calcvelocity(spatx, slopedat)

## ----fig.width = 10, fig.height = 3--------------------------------------
rtrend <- rgrad <- rvocc <- raster(sst)
rgrad[spatx$icell] <- spatx$NS
rtrend[slopedat$icell] <- slopedat$slope
rvocc[velodf$icell] <- velodf$velocity

par(mfrow = c(1,3))
plot(rtrend)
plot(rgrad)
plot(rvocc)


## ------------------------------------------------------------------------
citation("vocc")

