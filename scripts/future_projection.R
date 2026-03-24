# This script prepares data for projecting future NEE trends.
# Junna Wang, 03-17-2026
# 
# global monthly CO2 is obtained from https://zenodo.org/records/5021361
# citation: Cheng, W., Dan, L., Deng, X., Feng, J., Wang, Y., Peng, J., Tian, J., Qi, W., Liu, Z., Zheng, X., Zhou, D., Jiang, S., Zhao, H., & Wang, X. (2021). Global monthly distributions of atmospheric CO2 concentrations under the historical and future scenarios [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5021361
# for disturbance effects: use projections based various climate variables; be consistent. 


rm(list=ls())
# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Volumes/WZZ_disk'
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
land_cover <- terra::rast(file.path('data', '968_Land_Cover_Class_0.25degree_layer2.tif'))
# we only extract layer 2 of the 968_Land_Cover_Class_0.25degree.nc file. 
#
# linked information about global land_cover_class maps
# https://thredds.daac.ornl.gov/thredds/ncss/grid/ornldaac/968/Land_Cover_Class_0.25degree.nc4/dataset.html
# resolution 0.25 degree
# type 12, croplands; type 13, urban areas; 
# type_1 = Evergreen Needleleaf Forests; type2: Evergreen Broadleaf Forests; type3: Deciduous Needleleaf Forests;  
# type_4 = Deciduous Broadleaf Forests; type_5 = Mixed Forests
# type_11 = Permanent Wetlands; 
# 0water, 12cropland, 13urban, 15ice, 16barren

# plot(koppen_climate)
# plot(land_cover%in% c(4:5))
# large croplands on earth

# unify resolution: 0.25 degree
koppen_climate_0.25 <- terra::resample(koppen_climate, land_cover, method="near")
# plot(koppen_climate_0.25)

# aridity index
AI_0.25 <- terra::rast(file.path('data', 'Global-AI.tif'))
AI_0.25 <- terra::resample(AI_0.25, land_cover) / 10000       # values are enlarged by 10000
# Zomer, R.J., Xu, J. & Trabucco, A. Version 3 of the Global Aridity Index and Potential Evapotranspiration Database. Sci Data 9, 409 (2022). https://doi.org/10.1038/s41597-022-01493-1
#
# a data frame to store zone info
df_zone <- data.frame(koppen_climate=as.integer(terra::values(koppen_climate_0.25)),
                      land_cover=as.integer(terra::values(land_cover)),
                      AI = terra::values(AI_0.25))
names(df_zone)[3] <- "AI"

# the classification of dry and wet in tropics should based on aridity index. 
df_zone <- df_zone %>% mutate(group = case_when(
  koppen_climate %in% c(1:3, 11, 14) & !land_cover %in% c(12,13,15,16,0) & AI > 0.5 ~ "wet", 
  koppen_climate %in% c(1:3, 11, 14) & !land_cover %in% c(12,13,15,16,0) & AI < 0.5 ~ "dry", 
  koppen_climate %in% c(4:10) & !land_cover %in% c(12,13,15,16,0) ~ "dry",  # remove 0water, 12cropland, 13urban, 15ice, 16barren
  koppen_climate %in% c(12:13, 15:30) & land_cover %in% c(6:9, 14) ~ "cold_nonforest", 
  koppen_climate %in% c(12:13, 15:30) & land_cover %in% c(10) ~ "cold_nonforest", 
  koppen_climate %in% c(12:13, 15:30) & land_cover %in% c(11) ~ "cold_nonforest", 
  koppen_climate %in% c(12:13, 15:30) & land_cover %in% c(1:2) ~ "cold_ENF",
  koppen_climate %in% c(12:13, 15:30) & land_cover %in% c(3:5) ~ "cold_DBF"))


zone_grow_season <- koppen_climate_0.25

# zone for climate-vegetation groups
zone_group <- koppen_climate_0.25
terra::values(zone_group) <- NA
terra::values(zone_group) <- df_zone$group
plot(zone_group)
# why there are so many missing points? these are croplands, urban, ice, and barren.  

writeRaster(zone_group, filename = file.path('data', 'zones_climate_vegetation_group.tif'), overwrite=TRUE)

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

# Ideas about how to deal with northern and southern growing seasons.

###########################################################
# 3. a function to calculate growing season mean by zones
###########################################################

process_one_year <- function(r, zone_masks_north, zone_masks_south, sos_eos_group, gs=TRUE) {
  output <- terra::rast(r[[1]])
  terra::values(output) <- NA
  
  for (rid in 1:nrow(sos_eos_group)) {
    #    browser()
    days_north <- sos_eos_group$sos_group[rid]:sos_eos_group$eos_group[rid]
    days_south <- ifelse(days_north + 181 > 365, days_north + 181 - 365, days_north + 181)
    
    if (!gs) {
      days_north <- setdiff(1:365, days_north)
      days_south <- setdiff(1:365, days_south)
    }
    
    # Subset days for the region
    if (length(days_north) == 0) {next}
    r_north <- r[[days_north]]
    r_south <- r[[days_south]]
    
    # Compute mean across time
    # mean_north <- terra::app(r_north, mean, na.rm = TRUE, cores=8)
    # mean_south <- terra::app(r_south, mean, na.rm = TRUE, cores=8)
    # 
    # # mask results; this is faster than mask first
    # mean_north <- mask(mean_north, zone_masks_north[[rid]])
    # mean_south <- mask(mean_south, zone_masks_south[[rid]])
    
    # mean_north <- zonal(r_north, zone_masks_north[[rid]], fun = "mean", na.rm = TRUE)
    # mean_south <- zonal(r_south, zone_masks_south[[rid]], fun = "mean", na.rm = TRUE)
    cells_north <- which(!is.na(terra::values(zone_masks_north[[rid]])))
    cells_south <- which(!is.na(terra::values(zone_masks_south[[rid]])))
    #
    if (length(cells_north) > 0) {
      output[cells_north] <- rowMeans(r_north[cells_north], na.rm=T)
    }    
    #
    if (length(cells_south) > 0) {
      output[cells_south] <- rowMeans(r_south[cells_south], na.rm=T)
    }
    #
    # Merge them into final output
    # output <- cover(output, mean_north)
    # output <- cover(output, mean_south)
    # print(paste0('finish zone ', rid))
  }
  return(output)
}

################################################################################
# 4. read daily climate data and calculate growing season mean for each year
################################################################################
met_files <- list.files(file.path(dir_rawdata, '/ScaleUp_Global/met_day_2005_2023/'), pattern='.nc$', full.names = T)
#
sw_data <- list()
tam_data <- list()
tp_data <- list()
#
df_file_info <- data.frame(year=integer(), var=character(), file=character())
#
for (i in 1:length(met_files)) {
  nc <- nc_open(met_files[i])
  time_str <- ncatt_get(nc, "valid_time", "units")$value
  year <- as.numeric(sub(".*since (\\d{4})-.*", "\\1", time_str))
  varname <- names(nc$var)[1]
  df_file_info[i,1] <- year
  df_file_info[i,2] <- varname
  df_file_info[i,3] <- met_files[i]
  nc_close(nc)
}
df_file_info <- df_file_info %>% arrange(var, year)

# it takes 7-8 hours to rotate these rasters.
# work on t2m
# t2m_data <- list()
t2m_gs <- rast()
t2m_ngs <- rast()
df_file_info_subset <- subset(df_file_info, var=='t2m')    # nrow(df_file_info_subset)
for (i in 1:nrow(df_file_info_subset)) {
  print( paste(i, df_file_info_subset$var[i], df_file_info_subset$year[i], sep = '-') )
  met_data <- terra::rast(df_file_info_subset$file[i])
  met_data <- terra::rotate(met_data)
  # t2m_data <- c(t2m_data, met_data)
  if (i == 1) {
    # resample masks to be consistent with nc data
    zone_masks_north <- terra::resample(zone_masks_north, met_data, method = "near")
    zone_masks_south <- terra::resample(zone_masks_south, met_data, method = "near")
  }
  test <- process_one_year(met_data, zone_masks_north, zone_masks_south, sos_eos_group, gs=TRUE)
  t2m_gs <- c(t2m_gs, test)
  
  test <- process_one_year(met_data, zone_masks_north, zone_masks_south, sos_eos_group, gs=FALSE)
  t2m_ngs <- c(t2m_ngs, test)
  
  rm(met_data)
  gc()
  terra::tmpFiles(remove = TRUE)
}
writeRaster(t2m_gs, filename = file.path('data', 't2m_gs_2005_2023.tif'), overwrite=TRUE)
writeRaster(t2m_ngs, filename = file.path('data', 't2m_ngs_2005_2023.tif'), overwrite=TRUE)

#-----------------------------------
# tp_data <- list()
tp_gs <- rast()
tp_ngs <- rast()
df_file_info_subset <- subset(df_file_info, var=='tp')
for (i in 1:nrow(df_file_info_subset)) {
  print( paste(i, df_file_info_subset$var[i], df_file_info_subset$year[i], sep = '-') )
  met_data <- terra::rast(df_file_info_subset$file[i])
  met_data <- terra::rotate(met_data)
  # tp_data <- c(tp_data, met_data)
  test <- process_one_year(met_data, zone_masks_north, zone_masks_south, sos_eos_group, gs=TRUE)
  tp_gs <- c(tp_gs, test)
  
  test <- process_one_year(met_data, zone_masks_north, zone_masks_south, sos_eos_group, gs=FALSE)
  tp_ngs <- c(tp_ngs, test)
  
  rm(met_data)
  gc()
  terra::tmpFiles(remove = TRUE)
}
writeRaster(tp_gs, filename = file.path('data', 'tp_gs_2005_2023.tif'), overwrite=TRUE)
writeRaster(tp_ngs, filename = file.path('data', 'tp_ngs_2005_2023.tif'), overwrite=TRUE)

#-----------------------------------
# ssrd_data <- list()
ssrd_gs <- rast()
df_file_info_subset <- subset(df_file_info, var=='ssrd')
for (i in 1:nrow(df_file_info_subset)) {
  print( paste(i, df_file_info_subset$var[i], df_file_info_subset$year[i], sep = '-') )
  met_data <- terra::rast(df_file_info_subset$file[i])
  met_data <- terra::rotate(met_data)
  # ssrd_data <- c(ssrd_data, met_data)

  test <- process_one_year(met_data, zone_masks_north, zone_masks_south, sos_eos_group, gs=TRUE)
  ssrd_gs <- c(ssrd_gs, test)

  rm(met_data)
  gc()
  terra::tmpFiles(remove = TRUE)
}
writeRaster(ssrd_gs, filename = file.path('data', 'ssrd_gs_2005_2023.tif'), overwrite=TRUE)

###########################################################
# 5. process ET data from FLUXCOM in the past
###########################################################
ET_files <- list.files(file.path(dir_rawdata, '/ScaleUp_Global/et_day_2005_2021/'), pattern='.nc$', full.names = T)

et_gs <- rast()
for (i in 1:length(ET_files)) {
  print(i)
  met_data <- terra::rast(ET_files[i])
  # plot(met_data)
  if (i == 1) {
    # resample masks to be consistent with nc data
    zone_masks_north <- resample(zone_masks_north, met_data, method = "near")
    zone_masks_south <- resample(zone_masks_south, met_data, method = "near")
  }
  #
  test <- process_one_year(met_data, zone_masks_north, zone_masks_south, sos_eos_group)
  et_gs <- c(et_gs, test)
  
  rm(met_data)
  gc()
  terra::tmpFiles(remove = TRUE)
}
writeRaster(et_gs, filename = file.path('data', 'et_gs_2005_2023.tif'), overwrite=TRUE)

######################################################################
# 6. get future climate average during growing and non-growing season
######################################################################
# get the average of the three global circulation models: GFDL, MICRO6, and MPI-ESM
vars <- c('pr', 'rsds', 'tas')
for (var in vars) {
  folder <- paste0(var, '_day_2026_2050')
  #
  nc_files <- list.files(file.path(dir_rawdata, 'ScaleUp_Global', folder), pattern='.nc$', full.names = T)
  tp_data_fut3 <- list()
  for (i in 1:length(nc_files)) {
    met_data <- terra::rast(nc_files[i])
    tp_data_fut3 <- c(tp_data_fut3, met_data)
  }
  
  tp_data_fut <- list()
  nyears <- length(nc_files) / 3
  # we only used the data between 2026 and 2045
  for (i in 1:20) {
    # for (i in 1:nyears) {
    tp_data_fut <- c(tp_data_fut, (tp_data_fut3[[i]] + tp_data_fut3[[i+nyears]] + tp_data_fut3[[i+nyears*2]]) / 3)
    print(i)
  }
  
  # resample masks to be consistent with nc data
  test <- terra::rotate(tp_data_fut[[1]])
  zone_masks_north <- resample(zone_masks_north, test, method = "near")
  zone_masks_south <- resample(zone_masks_south, test, method = "near")
  
  # average tp of growing season and non-growing season
  tp_gs_fut <- rast()
  for (i in 1:length(tp_data_fut)) {
    tp_data_fut[[i]] <- terra::rotate(tp_data_fut[[i]])
    test <- process_one_year(tp_data_fut[[i]], zone_masks_north, zone_masks_south, sos_eos_group, gs=TRUE)
    tp_gs_fut <- c(tp_gs_fut, test)
    print(i)
  }
  writeRaster(tp_gs_fut, filename = paste0('data/', var, '_gs_2026_2050.tif'), overwrite=TRUE)
  
  # non-growing season
  if (var %in% c('pr', 'tas')) {
    tp_ngs_fut <- rast()
    for (i in 1:length(tp_data_fut)) {
      test <- process_one_year(tp_data_fut[[i]], zone_masks_north, zone_masks_south, sos_eos_group, gs=FALSE)
      tp_ngs_fut <- c(tp_ngs_fut, test)
      print(i)
    }
    writeRaster(tp_ngs_fut, filename = paste0('data/', var, '_ngs_2026_2050.tif'), overwrite=TRUE)      
  }
}

######################################################################
# 7. get change of ET from coarse global climate models
######################################################################
# I want to get the change of ET from 2015-2025 to 2026-2050. 
# assuming growing season is 5-9; 
# calculate changes in slopes using ET projections from 6 global climate models under ssp245;
# BCC-CSM2-MR, CAMS-CSM1-0, CAS-ESM2-0, FIO-ESM-2-0, MIROC6, CESM2_ssp245
folder <- 'et_day_2026_2050'
nc_files <- list.files(file.path(dir_rawdata, 'ScaleUp_Global', folder), pattern='.nc$', full.names = T)
et_data_fut3 <- list()
for (i in 1:length(nc_files)) {
  met_data <- terra::rast(nc_files[i])
  et_data_fut3 <- c(et_data_fut3, met_data)
}
et_data_fut3[[6]] <- c(et_data_fut3[[6]], et_data_fut3[[7]])

nmodels <- 6

et_slope_model_past <- list()
et_slope_model_fut <- list()
et_slope_model_past_fut_change <- list()
for (imodel in 1:nmodels) {
  et_data_year <- rast()
  data_north <- terra::mask(et_data_fut3[[imodel]], init(et_data_fut3[[imodel]], "y") > 0, maskvalues=0)
  data_south <- terra::mask(et_data_fut3[[imodel]], init(et_data_fut3[[imodel]], "y") <= 0, maskvalues=0)
  #
  for (iyear in 1:36) {
    r_north <- terra::app(data_north[[c(5:9) + (iyear-1)*12]], fun=mean, na.rm=T)
    r_south <- terra::app(data_south[[c(1:3, 11:12) + (iyear-1)*12]], fun=mean, na.rm=T)
    et_data_year <- c(et_data_year, terra::cover(r_north, r_south))
    print(iyear)
  }
  # 2015-2030
  et_slope_model_past[[imodel]] <- terra::app(et_data_year[[1:16]], function(x) {ifelse(any(is.na(x)), NA, sens.slope(x)$estimate / sd(x))})
  #  plot(et_slope_model_past[[imodel]])
  et_slope_model_past[[imodel]] <- terra::rotate(et_slope_model_past[[imodel]])
  # 2031-2046
  et_slope_model_fut[[imodel]] <- terra::app(et_data_year[[17:32]], function(x) {ifelse(any(is.na(x)), NA, sens.slope(x)$estimate / sd(x))})
  et_slope_model_fut[[imodel]] <- terra::rotate(et_slope_model_fut[[imodel]])
}


















# Compare the past ET trends from FLUXCOM with the current ET trends from climate models
et_gs_slope_FLUXCOM <- terra::rast('data/et_gs_slope_FLUXCOM_2005_2023.tif')
plot(et_gs_slope_FLUXCOM)
for (imodel in 1:nmodels) {
  et_slope_model_past[[imodel]] <- terra::resample(et_slope_model_past[[imodel]], et_gs_slope_FLUXCOM)
  et_slope_model_past[[imodel]] <- terra::mask(et_slope_model_past[[imodel]], et_gs_slope_FLUXCOM)
  plot(et_slope_model_past[[imodel]])
  print(cor.test(terra::values(et_gs_slope_FLUXCOM), terra::values(et_slope_model_past[[imodel]])))
  #
  et_slope_model_fut[[imodel]] <- terra::resample(et_slope_model_fut[[imodel]], et_gs_slope_FLUXCOM)
  et_slope_model_fut[[imodel]] <- terra::mask(et_slope_model_fut[[imodel]], et_gs_slope_FLUXCOM)
}

et_slope_model_past <- terra::app(terra::rast(et_slope_model_past), fun='mean', na.rm=T)
et_slope_model_fut  <- terra::app(terra::rast(et_slope_model_fut), fun='mean', na.rm=T)

print(cor.test(terra::values(et_gs_slope_FLUXCOM), terra::values(et_slope_model_past)))
# r = 0.20099 and significant


et_slope_model_past_fut_change <- et_slope_model_fut - et_slope_model_past
plot(et_slope_model_past_fut_change)
#
writeRaster(et_slope_model_past_fut_change, filename = file.path('data', 'et_slope_past_fut_change.tif'), overwrite=TRUE)

##################################
# 8. calculate past NEE trends
##################################




##################################
# 10. predict for future NEE trend
##################################







