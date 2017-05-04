#' Calculate trend in temperature for pixels in a raster brick
#'
#' @usage calcslope(x, divisor = 10, na.rm = TRUE)
#'
#' @param rx \code{brick} of rasters giving temperature at years.
#' @param divisor \code{numeric} giving divisor for rate.
#' @param na.rm \code{logical} should NAs be removed?
#' @param output either \code{"data.frame"} or \code{"raster"}, see details.
#'
#' @return A \code{data.frame} giving the slope at each pixel and number of
#' points in each pixel.
#'
#' @details Calculates the trend in temperature.
#' \code{divisor} is for getting correct units, e.g. if data is annual then
#' \code{divisor = 10} gives degrees C per decade.
#'This method is pretty fast. It
#' takes about 2.5 seconds on a 3Ghz laptop to proces 64800 cells over
#' 80 years of annual data.
#' Use \code{output = "data.frame"} to generate a \code{data.frame} that can
#' be used in \code{\link{calcvelocity}} or \code{output = "raster"} to
#' generate a \code{raster} of the trend field.
#'
#' @author David Schoeman, Christopher Brown
#' @examples
#' data(sst)
#' dat <- calcslope(sst)
#' head(dat$slope)
#' @rdname calcslope
#' @export

calcslope<-function(rx, divisor = 10, na.rm = TRUE, output = "data.frame"){
    if (! (output %in% c("data.frame", "raster"))) stop("Invalid term used in
        the output argument. ")
    icell <- 1:ncell(rx)
	lonlat <- xyFromCell(rx, icell)
    y <- t(getValues(rx))
    x <- row(y)
    #x <- row(x)
    x <- x/divisor
    x1 <- y
    x1[is.na(x1) == F] <- 1
    N <- apply(x1, 2, sum, na.rm = na.rm)
    x <- x * x1
    rm(x1)
    xy <- x*y
    sxy <- apply(xy, 2, sum, na.rm = na.rm)
    rm(xy)
    x2 <- x*x
    sx2 <- apply(x2, 2, sum, na.rm = na.rm)
    rm(x2)
    sx <- apply(x, 2, sum, na.rm = na.rm)
    sy <- apply(y, 2, sum, na.rm = na.rm)
    slope <- (sxy-(sx*sy/N))/(sx2-((sx^2)/N))
    if (output == "data.frame"){
        return(data.frame(slope = slope, N = N, lonlat, icell))}
        else if(output == "raster"){
            rtrend <- raster(rx)
        	rtrend[icell] <- slope
            return(rtrend)
        }
	}
