## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE,
comment = "#>")
knitr::opts_knit$set(root.dir = "/Users/gnj/Dropbox/R-Projects/map.displ.r/inst/extdata/")

## ------------------------------------------------------------------------
library(raster)

dem2012 <- raster("rock_glacier_dem_2012_50cm.tif")
dem2017 <- raster("rock_glacier_dem_2017_50cm.tif")

## ------------------------------------------------------------------------
library(map.displ.r)
f_tx <- "RAW_rock_glacier_hillshade_2012_50cm_direct_transf.txt"
d_tx <- dem.displacement.mapping(tx_file = f_tx, r_source = dem2012, r_target = dem2017)
head(d_tx)

## ----  fig.height=6, fig.width = 6---------------------------------------

# Assign displacement values to an emptry raster to export
r_na <- setValues(dem2012, NA)
disp_3d <- setValues(r_na, d_tx$xyz_disp)
plot(disp_3d)

## ----  fig.height=6, fig.width = 6---------------------------------------
library(ggplot2)
library(metR)
library(ggnewscale)

# Load hillshade for our map of displacements
hs2017 <- raster("rock_glacier_hillshade_2017_50cm.tif")
hs2017_df <- as.data.frame(hs2017, xy=TRUE)
names(hs2017_df) <- c("x", "y", "hs")

# Estimate mean annual surface velocity (m/yr) over the five year period
d_tx$mean_disp <- d_tx$xyz_disp/5

# Make displacement map
map <- ggplot(d_tx, aes(x_source,y_source)) +
  geom_raster(data=hs2017_df, aes(x=x, y=y, fill = hs), show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") + 
  # Allow for multiple scale fills using ggnewscale package
  new_scale("fill") +
  
  geom_raster(data=d_tx, alpha = 0.4, aes(x=x_source, y=y_source, fill = mean_disp)) +
  scale_fill_viridis_c(name = "Mean annual\nsurface velocity\n(m/yr)", direction = 1) +
  
  # Create arrows using metR package
  geom_arrow(data=d_tx, aes(dx = x_disp, dy = y_disp), skip = 30, show.legend = FALSE) +
  
  xlab("Easting (m)") +
  ylab("Northing (m)") +
  coord_fixed() +
  theme(text = element_text(size = 9), axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90))

map

