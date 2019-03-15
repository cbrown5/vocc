#' Calculate trend in temperature for pixels in a raster brick
#'
#' @param rx \code{brick} of rasters giving temperature at years.
#' @param divisor \code{numeric} giving divisor for rate.
#' @param na.rm \code{logical} should NAs be removed?
#'
#' @return A \code{data.frame} giving the slope at each pixel and number of
#' points in each pixel.
#' @importFrom raster xyFromCell getValues ncell
#'
#' @details Calculates the trend in temperature.
#' \code{divisor} is for getting correct units, e.g. if data is annual then
#' \code{divisor = 10} gives degrees C per decade.
#' This method is pretty fast. It
#' takes about 2.5 seconds on a 3Ghz laptop to proces 64800 cells over
#' 80 years of annual data.
#'
#' @author David Schoeman, Christopher Brown
#' @examples
#' data(sst)
#' dat <- calcslope(sst)
#' head(dat$slope)
#' @rdname calcslope
#' @export

calcslope <- function(rx, divisor = 10, na.rm = TRUE) {
  icell <- seq(1, ncell(rx))
  lonlat <- xyFromCell(rx, icell)
  y <- t(getValues(rx))
  x <- row(y)
  x <- x / divisor
  x1 <- y
  x1[!is.na(x1)] <- 1
  N <- apply(x1, 2, sum, na.rm = na.rm)
  x <- x * x1
  rm(x1)
  xy <- x * y
  sxy <- apply(xy, 2, sum, na.rm = na.rm)
  rm(xy)
  x2 <- x * x
  sx2 <- apply(x2, 2, sum, na.rm = na.rm)
  rm(x2)
  sx <- apply(x, 2, sum, na.rm = na.rm)
  sy <- apply(y, 2, sum, na.rm = na.rm)
  slope <- (sxy - (sx * sy / N)) / (sx2 - ((sx^2) / N))
  data.frame(slope = slope, N = N, lonlat, icell)
}
