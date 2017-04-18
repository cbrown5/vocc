#functions for calculating annual velocities
#CJ Brown and Dave. Schoeman
#Nov 2011

########################################
# Calculate spatial gradient data frame #
########################################
spatialgrad2 <- function(thisraster){
nlats<-nrow(thisraster)
nlons<-ncol(thisraster)
# Read the data into a raster
x <- thisraster # This is the SST field you want a spatial gradient for
rm(thisraster)
# Manipulate the data so that you get a column of focal and columns of adjacent cells
y <- data.frame(adjacent(x, 1:ncell(x), 8)) # Collect adjacent cells for ALL cells
#y$sstFocal <- getValues(x)[y$from] # Get sst from focal cell
y <- y[order(y$from, y$to),] # Order them in the order of the raster
y <- na.omit(y) # Work with the smallest logical dataset (for speed only)
y$sst <- getValues(x)[y$to] # Insert SST for each target cell
y$sy <- rowFromCell(x, y$from)-rowFromCell(x, y$to) # Set up a column to identify rows in the raster (N = 1, mid = 0, S = -1)
y$sx <- colFromCell(x, y$to)-colFromCell(x, y$from) # Set up a column to identify columns in the raster (E = 1, mid = 0, W = -1)
y$sx[y$sx > 1] <- -1 # Sort out the W-E wrap at the dateline, part I
y$sx[y$sx < -1] <- 1 # Sort out the W-E wrap at the dateline, part II
y$code <- paste(y$sx, y$sy) # Make a unique code for each of the eight neighbouring cells
y$code1 <- ifelse(y$code == "1 0", "sstE", ifelse(y$code == "-1 0", "sstW", ifelse(y$code == "-1 1", "sstNW", ifelse(y$code == "-1 -1", "sstSW", ifelse(y$code == "1 1", "sstNE", ifelse(y$code == "1 -1", "sstSE", ifelse(y$code == "0 1", "sstN", "sstS"))))))) # Code cells with positions
y <- with(y, data.frame(from, code1, sst)) # Make a data frame that cast will work with
y <- cast(y, from ~ code1) # Cast the data so that each focal cell has it's neighbours
y$sstFocal <- getValues(x)[y$from] # Put sstFocal back in
y$LAT <- yFromCell(x, y$from) # The latitude of the focal cell

# Calculate individual spatial temperature gradients []Mike's exact method from email of 12 July 2011]: grads (º/km)
y$gradWE1 <- with(y,(sstN-sstNW)/(cos(rad(LAT+1))*111.325)) # Positive values indicate an increase in sst from W to E (i.e., in line with the Cartesian x axis)
y$gradWE2 <- with(y,(sstFocal-sstW)/(cos(rad(LAT))*111.325))
y$gradWE3 <- with(y,(sstS-sstSW)/(cos(rad(LAT-1))*111.325))
y$gradWE4 <- with(y,(sstNE-sstN)/(cos(rad(LAT+1))*111.325))
y$gradWE5 <- with(y,(sstE-sstFocal)/(cos(rad(LAT))*111.325))
y$gradWE6 <- with(y,(sstSE-sstS)/(cos(rad(LAT-1))*111.325))
y$gradNS1 <- with(y,(sstNW-sstW)/111.325) # Positive values indicate an increase in sst from S to N (i.e., in line with the Cartesian y axis)
y$gradNS2 <- with(y,(sstN-sstFocal)/111.325)
y$gradNS3 <- with(y,(sstNE-sstE)/111.325)
y$gradNS4 <- with(y,(sstW-sstSW)/111.325)
y$gradNS5 <- with(y,(sstFocal-sstS)/111.325)
y$gradNS6 <- with(y,(sstE-sstSE)/111.325)

### *** ### CB added 16 Feb 2016
# Set gradients of zero to a small number (zeros cause NAs in VoCC calculation)
 #igrad <-12:23 #columns for gradients
 #izero <- y[,igrad]==0
 #izero[is.na(izero)] <- FALSE
 #tinynum <- min(abs(y[,igrad][!izero]), na.rm=T)
 #y[,igrad][izero] <- tinynum
### *** ###

# Write a function for to calculate weighted mean using sum products (the weighted.mean call is inaccurate for some reason)
	mnwm <- function (d) # SSTs are d; weights are wgt
	{
	    X <- array()
	    w <- array()
	    for (i in 1:nrow(d)) {
	        X[i] <- sum(c(d[i,1], d[i,2]*2, d[i,3], d[i,4], d[i,5]*2, d[i,6]), na.rm = T)
	        w[i] <- sum(c(d[i,1]/d[i,1], d[i,2]*2/d[i,2], d[i,3]/d[i,3],  d[i,4]/d[i,4], d[i,5]*2/d[i,5], d[i,6]/d[i,6]), na.rm = T) # Should be 1s & 2s, with NAs where sst is NA
	    }
	   return(X/w)
	}

# Calulate NS and WE gradients. NOTE: for angles to work (at least using simple positive and negative values on Cartesian axes), S-N & W-E gradients need to be positive)
y$WEgrad <- mnwm(cbind(y$gradWE1, y$gradWE2, y$gradWE3, y$gradWE4, y$gradWE5, y$gradWE6)) # Temperature gradient, increasing W to E (ºC/km)
y$NSgrad <- mnwm(cbind(y$gradNS1, y$gradNS2, y$gradNS3, y$gradNS4, y$gradNS5, y$gradNS6)) # Temperature gradient, increasing S to N

# Calculate angles of gradients (º) - adjusted for quadrant: NE = angle; SE = angle +180º; SW = angle + 180º; NW = angle + 360º (0º is North)
ang <- function(dx, dy){ # A function that calculates angle for each row (dx = WEgrad, y = NSgrad)
	angle <- array()
	for(i in 1:length(dx)){
		angle[i] <- ifelse(dy[i] < 0, 180+deg(atan(dx[i]/dy[i])), ifelse(dx[i] < 0, 360+deg(atan(dx[i]/dy[i])), deg(atan(dx[i]/dy[i]))))
		}
		return(angle)
	}
y$angle <- ang(y$WEgrad, y$NSgrad) # Angle (º)

# Merge the reduced file back into the main file to effectively undo the na.omit
from <- c(1:(nlats*nlons)) # Make ordered from cells
y <- merge(data.frame(from), y, by.x = "from", by.y = "from", all.x = T, incomparables = NA) # Undo the na.omit
return(y)
# Output the data
data.frame(WE = y$WEgrad, NS = y$NSgrad, angle = y$angle)

} #end of function
########################################
# CALCULATE LINEAR TREND OVER TIME
########################################
calcslope2 <- function(mnthmean){
	#takes a brick of temperatures at time (months/years) as input
	#returns a data frame with slopes and number of data points for each grid cell.
	lonlat<-xyFromCell(mnthmean,1:ncell(mnthmean))

y <- t(getValues(mnthmean)) # Write a transposed dataframe of the data; columns are cells, rows are months; this makes the regression calculations quick and easy (don't need to transpose, but it costs little time)
rm(mnthmean)

x <- y # Make a matrix of xs to match your ys
x <- row(x) # Turn it into a matrix of time references (where time reference runs from t = 1 to t = nrow(x))
x <- x/10 # Turn into fractions of a decade so that slopes are in the correct units (º/decade)
x1 <- y # Make a masking matrix to insert NAs
x1[is.na(x1) == F] <- 1 # Make NAs = NA and all else 1
N <- apply(x1, 2, sum, na.rm = T) # Cound number of valid data points in each cell (= n)
x <- x*x1 # Multiply the time matrix by the 1/NA mask = time OR NA (NA where temperature datum is missing)
# y <- x # Use this line to replace y with x in a test of the code - should provide a slope of 1
rm(x1)

# Do the stuff needed for regression analysis
xy <- x*y # Matrix of xy
sxy <- apply(xy, 2, sum, na.rm = T) # Sum(XY)
rm(xy)
x2 <- x*x # Matrix of x squared
sx2 <- apply(x2, 2, sum, na.rm = T) # Sum(x2)
rm(x2)
sx <- apply(x, 2, sum, na.rm = T) # Sum(x)
sy <- apply(y, 2, sum, na.rm = T) # Sum(y)
rm(y)
rm(x) # CAN keep N and use it to dump out regressions where there are too few data points, but more important for annual than monthly means?

slope <- (sxy-(sx*sy/N))/(sx2-((sx^2)/N)) # Regression slopes (ºC/decade)
# NOTE: This routine provides slopes of 1 for all non-NA cells if you replace y with x, so the method and implementation are correct
data.frame(slope,N,lonlat)
	}

########################################
# CALCULATE VELOCITY OF CLIMATE CHANGE FUNCTION #
#TAKES OUTPUT FROM SPATIAL GRADIENT FUNCTION (GRAD)
#AND TIME DIFFERENCE FUNCTION (SLOPE)
########################################
calcvelocity2 <- function(grad,slope){

#vmap <- raster(nrows=nlats,ncols=nlons) # Make a receiving raster
#slope <- data.frame(xyFromCell(vmap, 1:ncell(vmap)), slope) # Get lat & lon for each cell
slope$w <- 111.325*cos(rad(slope$y)) # Add weights by lat (same as using "area" from raster, but area requires unloading and reloading fBasics and Mass because of a clashing function)

# Next, calculate the vector sum of spatial gradients, then use this in calculating velocity
grd <- data.frame(NSold = grad$NS, WEold = grad$WE) # Make a data frame for the calculation
grd$NS <- ifelse(is.na(grd$NSold) == TRUE, 0, grd$NSold) # Make NAs zeroes for NS, so that they don't make the vector sum fall over resulting in an NA
grd$WE <- ifelse(is.na(grd$WEold) == TRUE, 0, grd$WEold) # Make NAs zeroes for WE, so that they don't make the vector sum fall over resulting in an NA
grd$NAsort <- ifelse(abs(grd$NS)+abs(grd$WE) == 0, NA, 1) # Make sure that the calculation returns NA if BOTH NS AND WE were NA
grd$Grad <- grd$NAsort*sqrt((grd$WE^2)+(grd$NS^2)) # Calculate the vector sum of gradients, ignoring the sign of the component gradients (º/km)

# Calculate velocity and output to file
velocity <- data.frame(x = slope$x, y = slope$y, temporal_trend = slope$slope, spatial_gradient = grd$Grad, NSgrad = grad$NS, WEgrad = grad$WE, angle = grad$angle, w = slope$w) # Make a separate data frame with the useful measures
velocity$velocity <- with(velocity, temporal_trend/spatial_gradient) # Calculate velocity of climate change (km per time unit)

velocity
}
