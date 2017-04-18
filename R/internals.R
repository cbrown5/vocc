#' @rdname .mnwm
.mnwm <- function(d1, d2, d3, d4, d5, d6){
    X <- sum(c(d1, d2*2, d3, d4, d5*2, d6), na.rm = T)
    w <- sum(c(1,2,1,1,2,1) * is.finite(c(d1, d2, d3, d4, d5, d6)))
   return(X/w)
}

#' @rdname .ang
.ang <- function(dx, dy){
	ifelse(dy < 0, 180 + CircStats::deg(atan(dx/dy)),
        ifelse(dx < 0, 360 + CircStats::deg(atan(dx /dy )), CircStats::deg(atan(dx/dy))))
	}
