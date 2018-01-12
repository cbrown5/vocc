# Create a data-set for package
# Hadley data that is already converted to brick. Cannot distribute this with final release.
#Including now for testing
# CJ Brown 18 Apr 2017

# Compare new VoCC calculations to old ones
# CJ Brown 18 Apr 2017
# TODO:
# There are some differences in cul-de-sacs, due to gradients being slightly different.
# My new code tends to calculate much greater values of vocc in these places.
#Overall this won't affect maps much, but would be good to know why difference occurs.

library(raster)
devtools::load_all("~/Code/geoengineering/vocc")

sst <- stack("~/Code/geoengineering/vocc/data-raw/hadley/Hadley_rcp45.gri")
sst <- crop(sst, extent(120, 180, -50, -8))

devtools::use_data(sst, pkg = "~/Code/geoengineering/vocc", overwrite = TRUE)
