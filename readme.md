# Functions for Calculating the Velocity of Climate Change

D Schoeman, CJ Brown  18 Apr 2017

Beta release

Functions based on [Burrows MT, Schoeman DS, Buckley LB, Moore P, Poloczanska ES, Brander KM, Brown C, Bruno JF, Duarte CM, Halpern BS, Holding J. The pace of shifting climate in marine and terrestrial ecosystems. Science. 2011 Nov 4;334(6056):652-5.](http://science.sciencemag.org/content/334/6056/652).

The package is designed to work with raster layers of temperature from the R package `raster`.

## CAUTION

Don't take it for granted that these functions will perform correctly! This is a beta release. We have attempted to optimise code from the originals used in our publication, and we are currently performing testing to ensure numerical accuracy before we release this to CRAN later in 2017.

As of 18 April 2017 there were some discrepencies between results from this code and the code from our paper.

Please contact [Chris Brown](christo.j.brown@gmail.com) if you want to use this code, or to provide feed-back.

## Installation

To install this beta version package, open R and type:

    install.packages("devtools")

Then, you can install vocc:

    devtools::install_github("cbrown5/vocc")

You migth like to start with:

    library(vocc)
    vignette("vocc")

And it should load. Let me know if you have troubles
