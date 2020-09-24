setwd("/Users/gnj/Dropbox/R-Projects/map.displ.r/test_data")
setwd("C:/Users/jg/Dropbox/R-Projects/map.displ.r/test_data")

## Load packages and data ###############################
library(raster)

#source raster
dem2012 <- raster("rock_glacier_dem_2012_50cm.tif")
#target raster
dem2017 <- raster("rock_glacier_dem_2017_50cm.tif")

# bunwarpj raw transformation file
f_tx <- "RAW_rock_glacier_hillshade_2012_50cm_direct_transf.txt"


## Determine 3D ground deformations #####################
d_tx <- dem.displacement.mapping(tx_file = f_tx, r_source = dem2012,
                                 r_target = dem2017)

## Scale deformations ###################################
d_est <- scale.displacements(d_tx, scale_factor = 0.13, return_sp = TRUE)
crs(d_est) <- crs(dem2012)

# Convert to grid
library(gstat)
spgdf <- as(dem2012, 'SpatialGridDataFrame')
grd_idw <- gstat::idw(z ~ 1, d_est, newdata=spgdf, idp=2.0, nmax = 20, nmin=12)
r_idw <- raster(grd_idw)


