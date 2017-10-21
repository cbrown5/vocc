#' Calculate spatial gradient across a raster layer
#'
#' @usage spatialgrad(rx, y_dist = c(111.325, 111.325), y_diff = 1)
#'
#' @param rx \code{raster} giving raster to calculate spatial gradient on
#' @param y_dist \code{numeric} giving distance for WE and NS gradient for
#' 1 unit of the y-coordinate
#' @param y_diff \code{numeric} giving units (degrees latitude) for the pixel
#' dimensions
#'
#' @return \code{data.frame} giving \code{icell} raster cell numbers,
#'   \code{WE} west-east gradient, \code{NS} north-south gradient,
#'   \code{angle} the angle of the trajectory.
#'
#' @details Assumes unprojected coordinates (ie lon-lat).
#' If you are using projected coordinates, with distance in metres,
#' you should give \code{y_dist = res(rx)} and \code{y_diff = NA}.
#' If \code{y_diff = NA}, no correction will be made for changing
#' size of 1 degree latitude with latitude.
#'
#' @author Christopher J. Brown
#' @examples
#' data(sst)
#' allyears <- rep(1, nlayers(sst))
#' mnsst <- stackApply(sst, indices = allyears, fun = mean)
#' spatx <- spatialgrad(mnsst)
#' head(spatx$NS)
#' @rdname spatialgrad
#' @export

spatialgrad <- function(rx, y_dist = c(111.325, 111.325), y_diff = 1){
    nlats <- nrow(rx)
    nlons <- ncol(rx)
    y <- data.frame(adjacent(rx, 1:ncell(rx), 8))
    y <- y[order(y$from, y$to),]
    y <- na.omit(y)
    y$sst <- getValues(rx)[y$to]
    y$sy <- rowFromCell(rx, y$from)-rowFromCell(rx, y$to)
    y$sx <- colFromCell(rx, y$to)-colFromCell(rx, y$from)
    y$sx[y$sx > 1] <- -1
    y$sx[y$sx < -1] <- 1
    y$code <- paste(y$sx, y$sy)

    y$code1 <- dplyr::recode(y$code,
    `1 0` = "sstE",
    `-1 0` = "sstW",
    `-1 1` = "sstNW",
    `-1 -1` = "sstSW",
    `1 1` = "sstNE",
    `1 -1` = "sstSE",
    `0 1` = "sstN",
    `0 -1` = "sstS")

    y3b <- y %>% dplyr::select(from, code1, sst) %>% tidyr::spread(code1, sst)
    y3b$sstFocal <- getValues(rx)[y3b$from]
    y3b$LAT <- yFromCell(rx, y3b$from)

    if(!is.na(y_diff)){
        y3b <- dplyr::mutate(y3b,
            latpos = cos(CircStats::rad(LAT + y_diff)),
            latneg = cos(CircStats::rad(LAT - y_diff)),
            latfocal = cos(CircStats::rad(LAT)))
        } else {
            y3b <- dplyr::mutate(y3b,
                latpos = 1,
                latneg = 1,
                latfocal = 1)
        }

    y3c <- dplyr::mutate(y3b,
        gradWE1 = (sstN-sstNW)/
            (latpos *  y_dist[1]),
        gradWE2 = (sstFocal - sstW)/(latfocal * y_dist[1]),
        gradWE3 = (sstS-sstSW)/(latneg * y_dist[1]),
        gradWE4 = (sstNE-sstN)/(latpos * y_dist[1]),
        gradWE5 = (sstE-sstFocal)/(latfocal * y_dist[1]),
        gradWE6 = (sstSE-sstS)/(latneg*y_dist[1]),
        gradNS1 = (sstNW-sstW)/y_dist[2],
        gradNS2 = (sstN-sstFocal)/y_dist[2],
        gradNS3 = (sstNE-sstE)/y_dist[2],
        gradNS4 = (sstW-sstSW)/y_dist[2],
        gradNS5 = (sstFocal-sstS)/y_dist[2],
        gradNS6 = (sstE-sstSE)/y_dist[2]
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        WEgrad = .mnwm(gradWE1, gradWE2, gradWE3, gradWE4, gradWE5, gradWE6),
        NSgrad = .mnwm(gradNS1, gradNS2, gradNS3, gradNS4, gradNS5, gradNS6),
        angle = .ang(WEgrad, NSgrad)) %>%
    dplyr::select(icell = from, WE = WEgrad, NS = NSgrad, angle = angle)

    return(y3c)
    }
