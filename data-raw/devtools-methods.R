# A few useful devtools things
# CJ Brown 2018-01-12

devtools::load_all("~/Code/geoengineering/vocc")

#Build vignettes and put them in the right folder
 devtools::build_vignettes("~/Code/geoengineering/vocc")

#do this first, so when we build the vignette it can find teh functions!
devtools::install("~/Code/geoengineering/vocc", build_vignette = FALSE)
devtools::build("~/Code/geoengineering/vocc")
