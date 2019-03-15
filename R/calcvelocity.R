#' Calculate velocity from the spatial gradient and temporal trends
#'
#' @param grad \code{data.frame} returned from \code{spatialgrad}
#' @param slope \code{data.frame} returned from \code{calcslope}
#' @param y_dist \code{numeric} giving units (degrees latitude) for the coordinate
#'
#' @return \code{data.frame} giving velocity.
#'
#' @author Christopher J. Brown
#' @examples
#' data(sst)
#' slopedat <- calcslope(sst)
#' allyears <- rep(1, raster::nlayers(sst))
#' mnsst <- raster::stackApply(sst, indices = allyears, fun = mean)
#' spatx <- spatialgrad(mnsst)
#' velodf <- calcvelocity(spatx, slopedat)
#' r <- raster::raster(sst)
#' r[velodf$icell] <- velodf$velocity
#' raster::plot(r)
#' @rdname calcvelociy
#' @export

calcvelocity <- function(grad, slope, y_dist = 111.325) {
  slope$w <- y_dist * cos(CircStats::rad(slope$y))
  grd <- data.frame(NSold = grad$NS, WEold = grad$WE)
  grd$NS <- ifelse(is.na(grd$NSold) == TRUE, 0, grd$NSold)
  grd$WE <- ifelse(is.na(grd$WEold) == TRUE, 0, grd$WEold)
  grd$NAsort <- ifelse(abs(grd$NS) + abs(grd$WE) == 0, NA, 1)
  grd$Grad <- grd$NAsort * sqrt((grd$WE^2) + (grd$NS^2))
  velocity <- data.frame(
    x = slope$x, y = slope$y, temporal_trend = slope$slope,
    spatial_gradient = grd$Grad, NSgrad = grad$NS, WEgrad = grad$WE,
    angle = grad$angle, w = slope$w, icell = grad$icell
  )
  velocity$velocity <- velocity$temporal_trend / velocity$spatial_gradient
  velocity
}
