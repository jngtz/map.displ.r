# PROCESSING BUNWARPJ RAW TRANSFORMATION FOR 3D DISPLACEMENT MAPPING
# WITH DIGITAL ELEVATION MODELS

# These functions are made for mapping 3D displacements of earth-surface
# processes (e.g. creep) from 2D image registration based on elastic
# deformations repsented by B-splines applied using digital elevation models.
# The resulting raw transformation file from the BUnwarpJ plugin in imagej
# is converted into the local coordinate system.

# Purpose: The free form deformation using b-splines can be performed using
#          derivatives of the digital elevation model (DEM) to improve
#          the registration results (e.g. a model of hillshade, slope,
#          surface roughness etc...). To run the bunwarpj plugin in imagej
#          the input source and target images (of equal row & col) must
#          be converted to a non-spaial image format such as a png. To
#          apply the resulting displacement field (i.e. raw trasnformation)
#          of the image registration to our geogrpahic data, the displacement
#          field must be converted to a coordinate reference system (CRS).
#          The main 'bunwarpjDisplacementField' function takes care of this.
#
#          Note the raw transformation file from the bunwarpj function
#          describes the registration by giving the registration coordinates
#          for each pixel of the source image by referencing the new row
#          and column position. This script is designed for use in a local
#          CRS.

# https://imagej.net/BUnwarpJ


## Euclidean distance #########################################

#' Calculate the Euclidean distance between two points.
#'
#' @param x1 x coordinate of point 1.
#' @param x2 x coordinate of point 2.
#' @param y1 y coordinate of point 1.
#' @param y2 y coordinate of point 2.
#' @return The Euclidean distance \code{x} and \code{y}.
#' @examples
#' euclidean.distance(2,3,4,5)

euclidean.distance <- function(x1, x2, y1, y2){
  sqrt( (x1 - x2)^2 + (y1 - y2)^2 )
}

## imageJ Coordinate transform. image to CRS ###########################

#' Define CRS for points from ImageJ/Fiji
#'
#' This function can transform the points saved from ImageJ/Fiji
#'     in image coordinates (col, row) with origin (1,1)
#'     to the  coordinate reference system (CRS) cooresponding
#'     to the image used for marking.
#'
#' @param file_name File name of text file with point coordinates that was
#'     saved from ImageJ
#' @param r A Raster* object

landmarks.to.crs <- function(file_name, r){

  # This function can convert the points saved from imagej (landmarks)
  # in image coordinates (col, row) with origin (1,1)
  # to the  coordinate reference system (CRS) cooresponding
  # to the image that are being used.

  # file_name: file name of text file with point coordinates
  #            that was saved from imageJ
  # r: Raster* object

  # define origin in local coord from top left corner
  x_min <- raster::xmin(r) # cell center at this point (not corner)
  y_max <- raster::ymax(r) # cell center (not corner)
  x_scale <- raster::res(r)[1]
  y_scale <- raster::res(r)[2]

  # Read image coordinates (col, rows) from output of SIFT in imagej
  img_crd_pnts <- read.delim(file_name, header=FALSE)
  names( img_crd_pnts) <- c("col", "row")

  # Shift  the origin of the image from (1,1) to (0,0)
  img_crd_pnts <-  img_crd_pnts-1

  img_crd_pnts$x_scale <-  img_crd_pnts$col * x_scale
  img_crd_pnts$y_scale <-  img_crd_pnts$row * y_scale

  # Calculate position in local CRS
  #   !note: x_scale/2 is set to make the CRS origin the cell corner
  #          instead of center

  img_crd_pnts$x_crs <- ( x_min - x_scale/2 + img_crd_pnts$x_scale)
  img_crd_pnts$y_crs <- ( y_max + y_scale/2 - img_crd_pnts$y_scale)

  return(data.frame(x=img_crd_pnts$x_crs, y=img_crd_pnts$y_crs))
}



## Matrix handling ########################################

#' Matrix to raster
#'
#' Converts a matrix to a raster object using the CRS and extent of a reference raster
#'
#' @param m Matrix.
#' @param r_ref A raster object from {raster} used as a reference
#'     for CRS and extent
#' @return a raster object projected in reference CRS and extent



matrix.to.raster <- function(m, r_ref){
  #require(raster)
  # Converts a matrix to a raster object using crs and
  # extent of a reference raster

  #   !matrix nrow and ncol must match ref raster

  # m: matrix
  # r_ref: a raster object from {raster} used as
  #        referece for CRS and and extent
  r_m <- raster::raster(m, xmn=r_ref@extent@xmin,
                xmx=r_ref@extent@xmax,
                ymn=r_ref@extent@ymin,
                ymx=r_ref@extent@ymax,
                crs=crs(r_ref))
  return(r_m)
}


## Angle handling #########################################

#' Calculate the direction component of the vector field.
#'
#' Calculate the direction component of the vector field
#'     with a geographic reference system - 0 degrees as due
#'     north and 90 degrees as due east (i.e. aspect)
#'
#' @param x x displacement magnitude
#' @param y y displacement magnitude
#' @return The aspect of displacement vectors \code{x} and \code{y}.
#' @examples
#' aspect.calc(3,4)

aspect.calc <- function(x,y){
  # Calculate the direction component of the vector field
  # with a geographic reference system - 0 degrees as due
  # north and 90 degrees as due east (i.e. aspect)

  if(x==0 & y==0){
    asp <- 0
  } else if(x==0){
    asp <- 0
  } else if(x>0 & y>=0){
    asp <- abs(atan(x/y))*180/pi
  } else if(x>0 & y<=0) {
    asp <- 180 - abs(atan(x/y))*180/pi
  } else if(x<0 & y>=0) {
    asp <- 360 - abs(atan(x/y))*180/pi
  } else if (x<0 & y<=0){
    asp <- 180 + abs(atan(x/y))*180/pi
  }
  return(asp)
}

#' Get reverse orientation
#'
#' Calculates direction in the reverse orientation
#'
#' @param angle an angle value in degrees
#' @return An angle in degrees of the opposite direction.
#' @examples
#' reverse.orientation(90)

reverse.orientation <- function(angle){
  # Calculate direction in a reverse orientation

  if(angle+180 >= 360){
    opp_asp <- angle - 180
  } else{
    opp_asp <- angle + 180
  }
  return(opp_asp)
}


## Read BUnwarpJ Raw transformation file ##################

#' 2D image displacement mapping from BUnwarpJ
#'
#' Perform 2D displacement mapping from raster images (or DEMs) using raw BUnwarpJ transformation file
#'
#' @param tx_file The raw transformation file (.txt) from BUnwarpJ.
#' @param r_source The DEM raster filename used as source for image registration.
#' @param r_target The DEM raster filename used as target for image registration.
#' @return Returns a data frame containing the following columns:
#' \item{x_source & y_source}{x & y coordinates from a CRS for the registration}
#' \item{x_target & y_target}{x & y coordinates from a CRS representing the x &
#' y transformation for every location in the source image}
#' \item{x_disp}{x displacement of registered image (from source image)}
#' \item{y_disp}{y displacement of registered image (from source image)}
#' \item{xy_disp}{magnitude of the displacement in 2D (xy)}
#' \item{aspect}{geographic orientation of vector direction with 0 degrees
#' being due north and 90 degrees being due east}


read.bunwarpj.displacements <- function(tx_file, r_source, r_target){
  ## Determine no. of columns and rows from file header
  dimensions <- as.character(read.table(tx_file, header=FALSE, nrows=2)[,1])
  dimensions <- as.numeric(gsub("^.*\\=", "",dimensions))
  n_cols <- dimensions[1]
  n_rows <- dimensions[2]

  ## Read the X and Y displacements from the raw transformation file
  x_raw_disp <- read.table(tx_file, header=FALSE, skip=4, nrows=n_rows)
  y_raw_disp <- read.table(tx_file, header=FALSE, skip=6+n_rows, nrows=n_rows)

  # Convert to matrix
  m_x_disp <- data.matrix(x_raw_disp, rownames.force=FALSE)
  m_y_disp <- unname(data.matrix(y_raw_disp, rownames.force=FALSE))
  remove(x_raw_disp, y_raw_disp, dimensions)

  ## Convert cell to cell description of displacement (from transformation
  #  file) to relative displacements using row and column numbers
  for (i in 1:ncol(m_x_disp)){
    m_x_disp[,i] <- m_x_disp[, i] - (i)
  }

  for (i in 1:nrow(m_y_disp)){
    m_y_disp[i,] <- (m_y_disp[i, ] - (i)) * -1 # Made negative to reflect UTM scale
  }

  ## Load source and target rasters
  #require(raster)
  resolution <- raster::res(r_source)[1]

  ## Scale matrix displacement to spatial resolution or raster cells
  m_x_disp <- m_x_disp * resolution
  m_y_disp <- m_y_disp * resolution

  ## Convert matrix to raster with corresponding extent and CRS
  r_x_disp <- matrix.to.raster(m_x_disp, r_source)
  r_y_disp <- matrix.to.raster(m_y_disp, r_source)

  ## Convert raster to SpatialPointsDataFrame
  # Using x displacements
  spdf <- raster::rasterToPoints(r_x_disp, spatial=TRUE)

  # Add y displacements
  spdf$y_disp <- values(r_y_disp)
  names(spdf@data) <- c('x_disp', 'y_disp')
  d <- as.data.frame(spdf)
  names(d) <- c('x_disp', 'y_disp', 'x_source', 'y_source')

  # Compute target position from x and y displacement values
  d$x_target <- d$x_disp + d$x_source
  d$y_target <- d$y_disp + d$y_source

  # Compute 2D (xy) displacement magnitude
  d$xy_disp <- sqrt( (d$x_disp)^2 + (d$y_disp)^2 )

  # Compute 2D (xy) dispalcement vector angles relative to North
  d$aspect <- mapply(aspect.calc, x=d$x_disp, y=d$y_disp)

  return(d)
}

## Import and process BUnwarpJ raw transformation file ####

#' 3D DEM displacement mapping from BUnwarpJ
#'
#' Perform 3D displacement mapping from DEMs using raw BUnwarpJ transformation file
#'
#' @param tx_file The raw transformation file (.txt) from BUnwarpJ.
#' @param r_source The DEM raster filename used as source for image registration.
#' @param r_target The DEM raster filename used as target for image registration.
#' @param is_inverse (\code{logical}) \code{TRUE} if image registraiton performed
#'      in the inverse direction, otherwise \code{FALSE}.
#' @return Returns a data frame containing the following columns:
#' \item{x_source & y_source}{x & y coordinates from a CRS for the registration}
#' \item{z_source}{the elevation value (intensity) corresponding the source
#' locations}
#' \item{x_target & y_target}{x & y coordinates from a CRS representing the x &
#' y transformation for every location in the source image}
#' \item{z_target}{the elevation value (intensity) corresponding the target
#' locations (from the target image /elevation model)}
#' \item{x_disp}{x displacement of registered image (from source image)}
#' \item{y_disp}{y displacement of registered image (from source image)}
#' \item{z_disp}{z displacement of registered image (from source image)
#' (z_target - z_source)}
#' \item{xy_disp}{magnitude of the displacement in 2D (xy)}
#' \item{xyz_disp}{magnitude of the displacement in 3D (xyz)}
#' \item{aspect}{geographic orientation of vector direction with 0 degrees
#' being due north and 90 degrees being due east}
#' \item{slope}{angle in degrees of displacement in z direction}

dem.displacement.mapping <- function(tx_file, r_source, r_target, is_inverse = FALSE){

  # tx_file: the raw transformation file (.txt) from BUnwarpJ
  # r_source: the DEM raster filename used as source for image registration
  # target_file: the DEM raster filename used as target for image registration

  # Returns a data frame

  ## Determine no. of columns and rows from file header
  dimensions <- as.character(read.table(tx_file, header=FALSE, nrows=2)[,1])
  dimensions <- as.numeric(gsub("^.*\\=", "",dimensions))
  n_cols <- dimensions[1]
  n_rows <- dimensions[2]

  ## Read the X and Y displacements from the raw transformation file
  x_raw_disp <- read.table(tx_file, header=FALSE, skip=4, nrows=n_rows)
  y_raw_disp <- read.table(tx_file, header=FALSE, skip=6+n_rows, nrows=n_rows)

  # Convert to matrix
  m_x_disp <- data.matrix(x_raw_disp, rownames.force=FALSE)

  m_y_disp <- unname(data.matrix(y_raw_disp, rownames.force=FALSE))

  remove(x_raw_disp, y_raw_disp, dimensions)

  ## Convert cell to cell description of displacement (from transformation
  #    file) to relative displacements using row and column numbers
  for (i in 1:ncol(m_x_disp)){
    m_x_disp[,i] <- m_x_disp[, i] - (i)
  }

  for (i in 1:nrow(m_y_disp)){
    m_y_disp[i,] <- (m_y_disp[i, ] - (i)) * -1 # Made negative to reflect UTM scale
  }

  ## Load source and target rasters
  #require(raster)
  #r_source <- raster(source_file)
  #r_target <- raster(target_file)

  # Get raster resolution
  resolution <- raster::res(r_source)[1]

  ## Scale matrix displacement to spatial resolution or raster cells
  m_x_disp <- m_x_disp * resolution
  m_y_disp <- m_y_disp * resolution

  ## Convert matrix to raster with corresponding extent and CRS
  r_x_disp <- matrix.to.raster(m_x_disp, r_source)
  r_y_disp <- matrix.to.raster(m_y_disp, r_source)

  # Reverse direction of displacements
  if(is_inverse == FALSE){
    r_x_disp <- r_x_disp * -1
    r_y_disp <- r_y_disp * -1
  }

  ## Convert raster to SpatialPointsDataFrame
  # Using x displacements
  spdf <- raster::rasterToPoints(r_x_disp, spatial=TRUE)

  # Add y displacements
  spdf$y_disp <- values(r_y_disp)
  names(spdf@data) <- c('x_disp', 'y_disp')

  ## Get elevation values from source DEM
  spdf$z_source <- raster::extract(r_source, spdf)

  ## Convert SpatialPointsDataFrame to data.frame
  d <- as.data.frame(spdf)
  names(d) <- c('x_disp', 'y_disp', 'z_source', 'x_source', 'y_source')

  # Compute target position from x and y displacement values
  d$x_target <- d$x_disp + d$x_source
  d$y_target <- d$y_disp + d$y_source

  # Compute 2D (xy) displacement magnitude
  d$xy_disp <- sqrt( (d$x_disp)^2 + (d$y_disp)^2 )

  # Compute 2D (xy) dispalcement vector angles relative to North
  d$aspect <- mapply(aspect.calc, x=d$x_disp, y=d$y_disp)

  # Get elevation values from target DEM at target positions
  # i.e. elevations at the corresponding points from the displaced
  # displaced source positions

  xy <- data.frame(x=d$x_target, y=d$y_target)
  d$z_target <- raster::extract(r_target, xy)


  # Compute z and xyz (3D) displacements
  d$z_disp <- d$z_target-d$z_source
  d$xyz_disp <- sqrt( (d$xy_disp)^2 + (d$z_disp)^2 )

  # Compute slope angle (degrees) of z displacement
  d$slope <- atan( d$z_disp/d$xy_disp )*180/pi


  ## Clean d for a nice output data.frame
  d_out <- data.frame(
    x_source=d$x_source,
    y_source=d$y_source,
    z_source=d$z_source,
    x_target=d$x_target,
    y_target=d$y_target,
    z_target=d$z_target,
    x_disp=d$x_disp,
    y_disp=d$y_disp,
    z_disp=d$z_disp,
    xy_disp = d$xy_disp, #aka total displacement for imcorr
    xyz_disp= d$xyz_disp,
    aspect=d$aspect,
    slope=d$slope
  )


  remove(d)
  return(d_out)

}

## Scale displacements #########################################

#' Scale displacements
#'
#' Applies a scale factor to adjust the magnitudes of displacement
#'     vectors
#'
#' @param d_tx a data frame containing the transformation vectors imported from bUnwarpJ
#' @param scale_factor the scale factor (\code{numeric}) to adjust the displacement magnitudes
#' @param xyz_disp  (\code{logical}) \code{TRUE} if x,y and z displacements
#'      are to be scaled, and \code{FALSE} if only x,y displacements are to be scaled.
#' @param return_sp  (\code{logical}) \code{TRUE} to return a
#'      SpatialPointsDataFrame based on the CRS projection of
#'      a reference \code{\link[raster]{raster}} object
#' @param r_ref a \code{\link[raster]{raster}} that has the CRS projection
#'      of the displacements
#' @return A data frame of the scaled displacements

displacements.scale <- function(d_tx, scale_factor, xyz_disp = TRUE, return_sp=FALSE,
                       r_ref = NA){
  # Computed scaled displacements
  d_scaled <- data.frame(x = rep(NA, len=nrow(d_tx)) )
  d_scaled$x <- d_tx$x_source + sin(d_tx$aspect*pi/180) * d_tx$xy_disp * scale_factor
  d_scaled$y <- d_tx$y_source + cos(d_tx$aspect*pi/180) * d_tx$xy_disp * scale_factor

  if(xyz_disp == TRUE){
    d_scaled$z <- d_tx$z_source + sin(d_tx$slope*pi/180) * d_tx$xyz_disp * scale_factor
  }

  if(return_sp == TRUE){
    if(xyz_disp == TRUE){
      d_scaled <- d_scaled[!is.na(d_scaled$z),]
    }
    raster::coordinates(d_scaled) <- ~x+y
    raster::crs(d_scaled) <- raster::crs(r_ref)
  }

  return(d_scaled)
}


## Transform raster using scaled displacements ###########################

#' Transform raster using scaled displacements
#'
#' The positions of values from a raster are transformed to another position
#'     based on a scale factor
#'
#' @param x a \code{\link[raster]{raster}} with values that should be
#'      transformed
#' @param d_tx a data frame containing the transformation vectors from
#'      obtained from the \code{bunwarpjDisplacementField()} function
#' @param scale_factor the scale factor (\code{numeric}) to adjust the displacement magnitudes
#' @param return_sp  (\code{logical}) \code{TRUE} to return a
#'      SpatialPointsDataFrame based on the CRS projection of the raster
#' @return A data frame of the scaled displacements


raster.transform <- function(x, d_tx, scale_factor, return_sp = FALSE){
  xy_crd <- data.frame(x=d_tx$x_source, y=d_tx$y_source)

  raster_values <- raster::extract(r_sd, xy_crd, method='simple')
  d_tx$raster_values <- raster_values

  # Calculate transformation of xy coordinates
  d_scaled <- data.frame(x = rep(NA, len=nrow(d_tx)) )
  d_scaled$x <- d_tx$x_source + sin(d_tx$aspect*pi/180) * d_tx$xy_disp * scale_factor
  d_scaled$y <- d_tx$y_source + cos(d_tx$aspect*pi/180) * d_tx$xy_disp * scale_factor
  d_scaled$sd <- d_tx$raster_values
  d_scaled <- d_scaled[!is.na(d_scaled$raster_values),]

  if(return_sp == TRUE){
    raster::coordinates(d_scaled) <- ~x+y
    raster::crs(d_scaled) <- raster::crs(r_ref)
  }

  return(d_scaled)
}


# This approach writes a tif format that imagej can read as an imput

## Functions ########################################################

#' GeoTIFF to TIFF format for ImageJ/Fiji
#'
#' ImageJ/Fiji does not always open GeoTIFF files. For example, written from
#'    ESRI products or the \code{\link[raster]{writeRaster}} function in
#'    \code{raster()}. This function gets around this issue by
#'    re-writing a problematic GeoTIFF using \code{rgdal()}.
#'
#' @param file_name File name of GeoTIFF raster to be used in ImageJ/Fiji.
#' @param write_name File name of new GeoTIFF raster to be used in ImageJ/Fiji.
#' @return A georeferenced image file (GeoTIFF) that can be read by ImageJ/Fiji.


geotiff.for.imagej <- function(file_name, write_name){
  #Convert a tif to a tif format that can be opened in imagej
  #require(rgdal)
  spgdf <- rgdal::readGDAL(file_name)
  rgdal::writeGDAL(spgdf, write_name)
}



