# Soil Moisture from Terra Climate:
library(sf)
library(terra)
library(climateR)
library(tidyverse)


sites_long <- read.csv(file.path('data', 'long_term_ec_sites.csv'))

df <- data.frame(
  site_ID = sites_long$site_ID, 
  lon = sites_long$Long,
  lat = sites_long$Lat
)

p <- terra::vect(df, geom = c("lon", "lat"), crs = "EPSG:4326")


soil_moisture_data_all <-data.frame()

for( i in 1:length(p)) {
  print(i)
  soil_moisture_data <- getTerraClim(
    AOI = p[i,],
    varname = "soil",       # soil water content (mm)
    startDate = "1990-01-01",           # 1990-01-01
    endDate = "2025-12-31") %>% mutate( site_ID = p$site_ID[i] )
  
  soil_moisture_data_all <- rbind( soil_moisture_data_all, soil_moisture_data)
  print("done")
}

write.csv(soil_moisture_data_all, '/Users/junnawang/Downloads/terra_SWC_allsites.csv', row.names = F )


