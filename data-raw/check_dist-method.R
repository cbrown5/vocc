# Run VoCC distance method
# CJ Brown 2018-01-12
# Based on Hamann et al. 2015 GCB

rm(list = ls())
library(RColorBrewer)
library(raster)
library(sf)
library(dplyr)

devtools::load_all("~/Code/geoengineering/vocc")
# devtools::build_vignettes("~/Code/geoengineering/vocc")
data(sst)
proj4string(sst) <- sp::CRS("+init=epsg:4326")

#
# Distance method functions
#

.getmindist <- function(inum, dat, distfun = st_distance, ...){
	dtemppre <- dat %>% dplyr::filter(prefact == inum)
	dtemppost <- dat %>% dplyr::filter(postfact == inum)

	if (nrow(dtemppost) == 0){
		dtemppre$dist <- Inf
	} else {
		dmat <- distfun(dtemppre, dtemppost, ...)
		dtemppre$dist <- apply(dmat, 1, min)
	}
	dtemppre
}

distvocc <- function(pre, post, tdiff, tol, denom = 1, distfun = st_distance, ...){

	names(pre) <- "prevar"
	names(post) <- "postvar"

	datpre <- raster::rasterToPoints(pre, spatial = T) %>%
		sf::st_as_sf()
	datpost <- raster::rasterToPoints(post, spatial = T) %>%
	sf::st_as_sf()

	dat <- sf::st_join(datpre, datpost)
	ymin <- min(c(datpre$prevar, datpre$postvar))
	ymax <- max(c(datpre$prevar, datpre$postvar))
	breaks <- seq(ymin-tol, ymax+tol, by = tol)

	dat$prefact <- cut(dat$pre, breaks = breaks, labels = FALSE)
	dat$postfact <- cut(dat$post, breaks = breaks, labels = FALSE)

	unq <- unique(dat$prefact)
	xout <- purrr::map(unq, ~.getmindist(.x, dat, distfun = st_distance, ...))

	dvocc <- do.call("rbind", xout)
	dvocc$vocc <- (dvocc$dist/denom)/tdiff
	spvocc <- as(dvocc, "Spatial")
	r <- raster::rasterize(spvocc, pre, field = "vocc")
	return(r)
}

#
# Implement
#

data(sst)
tol <- 0.5

pre <- raster(sst,1)
post <- raster(sst, 50)
tdiff <- 50
units <- 1000 #convert to km

rvocc <- distvocc(pre, post, tdiff, tol, denom = units)
plot(rvocc)
