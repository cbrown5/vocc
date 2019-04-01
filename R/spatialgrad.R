#' Calculate spatial gradient across a raster layer
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
#' @importFrom raster rowFromCell colFromCell getValues yFromCell
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
#' allyears <- rep(1, raster::nlayers(sst))
#' mnsst <- raster::stackApply(sst, indices = allyears, fun = mean)
#' spatx <- spatialgrad(mnsst)
#' head(spatx$NS)
#' @rdname spatialgrad
#' @export

spatialgrad <- function(rx, y_dist = c(111.325, 111.325), y_diff = 1) {
  if (length(y_dist) != 2L) stop("`y_dist` must be of length 2.", call. = FALSE)

  nlats <- nrow(rx)
  nlons <- ncol(rx)
  y <- data.frame(raster::adjacent(rx, seq(1, ncell(rx)), 8))
  y <- y[order(y$from, y$to), ]
  y <- stats::na.omit(y)
  y$sst <- getValues(rx)[y$to]
  y$sy <- rowFromCell(rx, y$from) - rowFromCell(rx, y$to)
  y$sx <- colFromCell(rx, y$to) - colFromCell(rx, y$from)
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
    `0 -1` = "sstS"
  )

  y3b <- dplyr::select(y, .data$from, .data$code1, .data$sst)
  y3b <- tidyr::spread(y3b, .data$code1, .data$sst)
  y3b$sstFocal <- getValues(rx)[y3b$from]
  y3b$LAT <- yFromCell(rx, y3b$from)

  if (!is.na(y_diff)) {
    y3b <- dplyr::mutate(y3b,
      latpos = cos(CircStats::rad(.data$LAT + y_diff)),
      latneg = cos(CircStats::rad(.data$LAT - y_diff)),
      latfocal = cos(CircStats::rad(.data$LAT))
    )
  } else {
    y3b <- dplyr::mutate(y3b,
      latpos = 1,
      latneg = 1,
      latfocal = 1
    )
  }

  y3c <- dplyr::mutate(y3b,
    gradWE1 = (.data$sstN - .data$sstNW) /
      (.data$latpos * y_dist[1]),
    gradWE2 = (.data$sstFocal - .data$sstW) / (.data$latfocal * y_dist[1]),
    gradWE3 = (.data$sstS - .data$sstSW) / (.data$latneg * y_dist[1]),
    gradWE4 = (.data$sstNE - .data$sstN) / (.data$latpos * y_dist[1]),
    gradWE5 = (.data$sstE - .data$sstFocal) / (.data$latfocal * y_dist[1]),
    gradWE6 = (.data$sstSE - .data$sstS) / (.data$latneg * y_dist[1]),
    gradNS1 = (.data$sstNW - .data$sstW) / y_dist[2],
    gradNS2 = (.data$sstN - .data$sstFocal) / y_dist[2],
    gradNS3 = (.data$sstNE - .data$sstE) / y_dist[2],
    gradNS4 = (.data$sstW - .data$sstSW) / y_dist[2],
    gradNS5 = (.data$sstFocal - .data$sstS) / y_dist[2],
    gradNS6 = (.data$sstE - .data$sstSE) / y_dist[2]
  )
  y3c <- dplyr::rowwise(y3c)
  y3c <- dplyr::mutate(y3c,
    WEgrad = .mnwm(
      .data$gradWE1, .data$gradWE2, .data$gradWE3,
      .data$gradWE4, .data$gradWE5, .data$gradWE6
    ),
    NSgrad = .mnwm(
      .data$gradNS1, .data$gradNS2, .data$gradNS3,
      .data$gradNS4, .data$gradNS5, .data$gradNS6
    ),
    angle = .ang(.data$WEgrad, .data$NSgrad)
  )
  y3c <- dplyr::select(y3c,
    icell = .data$from, WE = .data$WEgrad,
    NS = .data$NSgrad, angle = .data$angle
  )

  y3c
}
