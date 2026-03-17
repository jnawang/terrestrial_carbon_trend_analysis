# This script prepares data for projecting future NEE trends.
# Junna Wang, 03-17-2026
# 
# global monthly CO2 is obtained from https://zenodo.org/records/5021361
# citation: Cheng, W., Dan, L., Deng, X., Feng, J., Wang, Y., Peng, J., Tian, J., Qi, W., Liu, Z., Zheng, X., Zhou, D., Jiang, S., Zhao, H., & Wang, X. (2021). Global monthly distributions of atmospheric CO2 concentrations under the historical and future scenarios [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5021361
# for disturbance effects: use projections based various climate variables; be consistent. 


rm(list=ls())
# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'
#
library(librarian)
shelf(tidyverse, ggpubr, terra, ncdf4, future.apply, trend)

###########################################
# 1. get start and end of growing season
###########################################

data_grow_season <- read.csv(file.path('data', 'growing_season_yearly.csv'))
sites_long <- read.csv(file.path('data', 'long_term_ec_sites.csv'))
#
sos_eos <- data_grow_season %>% group_by(site_ID) %>% 
  summarise(sos=mean(sos_TRS+sos_DER, na.rm=T)/2, eos=mean(eos_TRS+eos_DER, na.rm=T)/2) %>% 
  left_join(sites_long[, c('site_ID', 'IGBP', 'Climate.class', 'category')]) %>% 
  mutate(group = case_when(category=='dry' ~ 'dry', category=='wet' ~ 'wet', 
                           category=='cold' & IGBP %in% c('ENF', 'EBF') ~ 'cold_ENF',
                           category=='cold' & IGBP %in% c('DBF', 'MF', 'DNF') ~ 'cold_DBF',
                           category=='cold' & !IGBP %in% c('ENF', 'EBF', 'DBF', 'MF', 'DNF') ~ 'cold_nonforest')) %>%
  arrange(group, Climate.class) 

sos_eos_group <- sos_eos %>% group_by(category, Climate.class) %>% 
  summarise(sos_group=round(mean(sos)), eos_group=round(mean(eos))) 
sos_eos_group$sos_group[sos_eos_group$category=='wet'] <- 1
sos_eos_group$eos_group[sos_eos_group$category=='wet'] <- 365

# adding three climate groups with year long growth
sos_eos_group <- rbind(sos_eos_group, data.frame(category = c("wet", "wet", "wet"),
                                                 Climate.class = c("Af", "Am", "Cwa"), 
                                                 sos_group = c(1, 1, 1), 
                                                 eos_group = c(365, 365, 365)))

# we added 10 climate classes whose start and end of growing season were given based on similar climates
sos_eos_group <- rbind(sos_eos_group, data.frame(category = c("cold", "cold", "cold", "cold", "cold", "cold", "cold", "cold", "cold", "cold"),
                                                 Climate.class = c("Cfc", "Csc", "Cwb", "Cwc", "Dsa", "Dsc", "Dsd", "Dwa", "Dwd", 'EF'), 
                       sos_group = sos_eos_group$sos_group[match(c('Cfb', 'Csb', 'Cfb', 'Cfb', "Dsb", "Dsb", "Dsb", "Dwb", "Dwc", 'ET'), sos_eos_group$Climate.class[1:10])], 
                       eos_group = sos_eos_group$eos_group[match(c('Cfb', 'Csb', 'Cfb', 'Cfb', "Dsb", "Dsb", "Dsb", "Dwb", "Dwc", 'ET'), sos_eos_group$Climate.class[1:10])]))
sos_eos_group <- sos_eos_group %>% arrange(category, Climate.class)
# we have 30 climates now. 


# match climate class and its coding ID
sos_eos_group <- sos_eos_group %>% mutate(climate_id = case_when(
  Climate.class == 'Dwa' ~ 21, Climate.class == 'Dwb' ~ 22,
  Climate.class == 'Dwc' ~ 23, Climate.class == 'Dwd' ~ 24,
  #
  Climate.class == 'Dfa' ~ 25, Climate.class == 'Dfb' ~ 26,
  Climate.class == 'Dfc' ~ 27, Climate.class == 'Dfd' ~ 28,
  Climate.class == 'ET' ~ 29, Climate.class == 'EF' ~ 30,
  #
  Climate.class == 'Csc' ~ 10, Climate.class == 'Csb' ~ 9, Climate.class == 'Csa' ~ 8, 
  Climate.class == 'BSk' ~ 7, Climate.class == 'BSh' ~ 6, 
  Climate.class == 'BWk' ~ 5, Climate.class == 'BWh' ~ 4,
  #
  Climate.class == 'Aw' ~ 3, Climate.class == 'Am' ~ 2, Climate.class == 'Af' ~ 1, 
  #
  Climate.class == 'Dsa' ~ 17, Climate.class == 'Dsb' ~ 18, 
  Climate.class == 'Dsc' ~ 19, Climate.class == 'Dsd' ~ 20, 
  #
  Climate.class == 'Cfa' ~ 14, Climate.class == 'Cfb' ~ 15, Climate.class == 'Cfc' ~ 16,
  #
  Climate.class == 'Cwa' ~ 11, Climate.class == 'Cwb' ~ 12, Climate.class == 'Cwc' ~ 13))


###########################################################
# 2. get global zonation based on climate and vegetation
###########################################################
koppen_climate <- terra::rast(file.path('data', 'koppen_geiger_0p00833333.tif'))
land_cover <- terra::rast(file.path('data', '968_Land_Cover_Class_0.25degree.nc'))
#
# plot(koppen_climate)
plot(land_cover[[2]]%in% c(4:5))
# large croplands on earth
# unify resolution: 0.25 degree
koppen_climate_0.25 <- terra::resample(koppen_climate, land_cover, method="near")
# plot(koppen_climate_0.25)
#
# a data frame to store zone info
df_zone <- data.frame(koppen_climate=as.integer(terra::values(koppen_climate_0.25)),
                      land_cover=as.integer(terra::values(land_cover[[2]])))
df_zone <- df_zone %>% mutate(group = case_when(
  koppen_climate %in% c(2, 14, 11:13) & !land_cover %in% c(12,13,15,16,0) ~ "wet",  
  koppen_climate %in% c(4:10) & !land_cover %in% c(12,13,15,16,0) ~ "dry",  # remove 0water, 12cropland, 13urban, 15ice, 16barren
  koppen_climate %in% c(15:30) & land_cover %in% c(6:9, 14) ~ "cold_nonforest", 
  koppen_climate %in% c(15:30) & land_cover %in% c(10) ~ "cold_nonforest", 
  koppen_climate %in% c(15:30) & land_cover %in% c(11) ~ "cold_nonforest", 
  koppen_climate %in% c(15:30) & land_cover %in% c(1:3) ~ "cold_ENF",
  koppen_climate %in% c(15:30) & land_cover %in% c(4:5) ~ "cold_DBF"))
# %>% mutate(grow_season_id = ifelse(is.na(koppen_climate) | is.na(group), NA, paste(group, koppen_climate, sep='-')))

zone_grow_season <- koppen_climate_0.25
terra::values(zone_grow_season)[!df_zone$koppen_climate %in% sos_eos_group$climate_id] <- NA   
plot(zone_grow_season)
# why there are so many missing points? these are water areas. 
# we have 20 climate zones, missing Aw, Am, Af, and Dsa, Dsb, Dsc, Dsd

# zone for climate-vegetation groups
zone_group <- koppen_climate_0.25
terra::values(zone_group) <- NA
terra::values(zone_group) <- df_zone$group
plot(zone_group)
writeRaster(zone_group, filename = 'data/zones_climate_vegetation_group.tif', overwrite=TRUE)

#------------------------------------
# calculate zone_grow_season masks
zone_masks <- rast()
for (i in 1:nrow(sos_eos_group)) {
  m <- zone_grow_season == sos_eos_group$climate_id[i]
  terra::values(m)[!terra::values(m)] <- NA
  names(m) <- paste0("zone_", sos_eos_group$climate_id[i])
  zone_masks <- c(zone_masks, m)
}

north_hemispher <- init(zone_masks[[1]], "y") > 0
south_hemispher <- init(zone_masks[[1]], "y") <= 0

zone_masks_north <- mask(zone_masks, north_hemispher, maskvalues=0)
zone_masks_south <- mask(zone_masks, south_hemispher, maskvalues=0)






# linked information about global land_cover_class maps
# https://thredds.daac.ornl.gov/thredds/ncss/grid/ornldaac/968/Land_Cover_Class_0.25degree.nc4/dataset.html
# resolution 0.25 degree
# type 12, croplands; type 13, urban areas; 
# type_1 = Evergreen Needleleaf Forests; type2: Evergreen Broadleaf Forests; type3: Deciduous Needleleaf Forests;  
# type_4 = Deciduous Broadleaf Forests; type_5 = Mixed Forests
# type_11 = Permanent Wetlands; 





