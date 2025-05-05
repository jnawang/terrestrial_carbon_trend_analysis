##' .. this function is designed to return land proportion given a lat and lon . ..
##'
##' @title get_land_proportion, this function will return land proportion.
##' @param input lat, lon
##' @return land proportion
##' @note In cycles of the ECMWF Integrated Forecasting System (IFS) from CY41R1 (introduced in May 2015) onwards, grid boxes where this parameter has a value above 0.5 can be comprised of a mixture of land and inland water but not ocean. Grid boxes with a value of 0.5 and below can only be comprised of a water surface. In the latter case, the lake cover is used to determine how much of the water surface is ocean or inland water.
##' @note In cycles of the IFS before CY41R1, grid boxes where this parameter has a value above 0.5 can only be comprised of land and those grid boxes with a value of 0.5 and below can only be comprised of ocean. In these older model cycles, there is no differentiation between ocean and inland water.
##' 

##' @author Boya ("Paul") Zhang
##' 


library(terra)
library(ecmwfr)

get_land_proportion <- function(nc_file, lat, lon) {
  
  if (lon < 0) {
    lon <- lon + 360
  }
  lat_nearest <- round(lat * 10) / 10
  lon_nearest <- round(lon * 10) / 10
  
  land_sea_mask <- rast(nc_file)
  if (is.null(land_sea_mask)) {
    stop("Error: Failed to read the NetCDF file.")
  }
  land_value <- extract(land_sea_mask, cbind(lon_nearest, lat_nearest))
  #print(land_value)
  return(land_value)
}

# Example usage
#nc_file <- "/Volumes/Malonelab/Research/ERA5_FLUX/Data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc"  # Your NetCDF file
#nc_file <- "z:/Research/ERA5_FLUX/Data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc"  
#latitude <- 25.76   # Example latitude
#longitude <- -80.19 # Example longitude (Miami, FL)
#longitude <- -61 ## somewhere in the ocean
#land_cover <- get_land_proportion(nc_file, latitude,longitude)
