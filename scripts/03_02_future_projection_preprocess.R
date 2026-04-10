# This script prepares data for projecting future NEE trends.
# Junna Wang, 03-17-2026
# 
# make the unit of different variables clear enough
# TA: C; P: mm/d; SW: W/m2; LE: W/m2; ET: mm/hr, and then should be converted to LE unit. 


# global monthly CO2 is obtained from https://zenodo.org/records/5021361
# citation: Cheng, W., Dan, L., Deng, X., Feng, J., Wang, Y., Peng, J., Tian, J., Qi, W., Liu, Z., Zheng, X., Zhou, D., Jiang, S., Zhao, H., & Wang, X. (2021). Global monthly distributions of atmospheric CO2 concentrations under the historical and future scenarios [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5021361
# for disturbance effects: use projections based various climate variables; be consistent. 


rm(list=ls())
# change this location where the raw data will be saved accordingly
dir_rawdata <- '/Volumes/WZZ_disk'
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'
#
library(librarian)
shelf(tidyverse, ggpubr, terra, ncdf4, future.apply, trend)

###########################################
# 1. get start and end of growing season
###########################################
source(file.path("scripts", "00_functions.R"))

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

sos_eos_group <- sos_eos_group %>% left_join(df_koppen_climate, by = 'Climate.class')


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
# 0water, 12cropland, 13urban, 14cropland_mosaic, 15 ice, 16barren
# for other types, see the table below. 

# plot(koppen_climate)
# plot(land_cover%in% c(4:5))
# large cropland areas on earth

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
  koppen_climate %in% c(1:3, 11, 14) & !land_cover %in% c(12,13,14,15,16,17,0) & AI > 0.5 ~ "wet", 
  koppen_climate %in% c(1:3, 11, 14) & !land_cover %in% c(12,13,14,15,16,17,0) & AI < 0.5 ~ "dry", 
  koppen_climate %in% c(4:10) & !land_cover %in% c(12,13,14,15,16,17,0) ~ "dry",  # remove 0water, 12cropland, 13urban, 15ice, 16barren
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
# work on t2m, unit: Kelvin
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
# work on tp, unit: m
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
# work on ssrd, unit: J m**-2
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
# original unit of ET data is mm/hr

et_gs_past <- rast()
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
  et_gs_past <- c(et_gs_past, test)
  
  rm(met_data)
  gc()
  terra::tmpFiles(remove = TRUE)
}
writeRaster(et_gs_past, filename = file.path('data', 'et_gs_2005_2023.tif'), overwrite=TRUE)

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

#########################################################################################################
# 7. get baseline (2005~2021) information about NEE, GPP, ET, MAT, MAP, climate and vegetation group 
#########################################################################################################
folder <- 'gpp_fluxcom_2005_2021'
nc_files <- list.files(file.path(dir_rawdata, 'ScaleUp_Global', folder), pattern='.nc$', full.names = T)
gpp_baseline <- terra::rast()
for (i in 1:length(nc_files)) {
  print(i)
  gpp_data <- terra::rast(nc_files[i])
  gpp_data <- terra::app(gpp_data, fun=sum, na.rm=T)
  gpp_baseline <- c(gpp_baseline, gpp_data)
}
gpp_baseline <- terra::app(gpp_baseline, fun=mean, na.rm=T)

writeRaster(gpp_baseline, filename = file.path('data', 'baseline_gpp_2005-2021.tif'), overwrite=TRUE)

#-----------
folder <- 'nee_fluxcom_2005_2021'
nee_baseline <- terra::rast()
for (i in 2005:2021) {
  print(i)
  nee_data <- terra::rast(file.path(dir_rawdata, 'ScaleUp_Global', folder, paste0('FLUXCOM-X-BASE net ecosystem exchange for ', i), paste0("NEE_", i, "_025_daily.nc")))
  nee_data <- terra::app(nee_data, fun=sum, na.rm=T)
  nee_baseline <- c(nee_baseline, nee_data) 
}
nee_baseline <- terra::app(nee_baseline, fun=mean, na.rm=T)

writeRaster(nee_baseline, filename = file.path('data', 'baseline_nee_2005-2021.tif'), overwrite=TRUE) 

#-----------
folder <- 'et_day_2005_2021'
nc_files <- list.files(file.path(dir_rawdata, 'ScaleUp_Global', folder), pattern='.nc$', full.names = T)
et_baseline <- terra::rast()
for (i in 1:length(nc_files)) {
  print(i)
  et_data <- terra::rast(nc_files[i])
  et_data <- terra::app(et_data, fun=mean, na.rm=T)
  et_baseline <- c(et_baseline, et_data)
}
et_baseline <- terra::app(et_baseline, fun=mean, na.rm=T)

writeRaster(et_baseline, filename = file.path('data', 'baseline_et_2005-2021.tif'), overwrite=TRUE) 

#--------MAT, MAP
# mat_baseline <- terra::rast(file.path('data', 'MAT.tif'))
# map_baseline <- terra::rast(file.path('data', 'MAP.tif'))
# 
# mat_baseline <- terra::resample(mat_baseline, gpp_baseline) 
# map_baseline <- terra::resample(map_baseline, gpp_baseline)
# 
# writeRaster(mat_baseline, filename = file.path('data', 'baseline_mat.tif'), overwrite=TRUE)
# writeRaster(map_baseline, filename = file.path('data', 'baseline_map.tif'), overwrite=TRUE)

# climate_baseline <-  koppen_climate
# igbp_baseline <-     land_cover

######################################################################
# 8. get change of ET from coarse global climate models and future ET
######################################################################
# I want to get the change of ET from 2015-2025 to 2026-2050. 
# assuming growing season is month 5-9; 
# calculate changes in slopes using ET projections from 6 global climate models under ssp245;
# BCC-CSM2-MR, CAMS-CSM1-0, CAS-ESM2-0, FIO-ESM-2-0, MIROC6, CESM2_ssp245
# the original data is monthly and the unit of evspsbl is mm/s; we convert this to mm/hr to be consistent with FLUXCOM unit
folder <- 'et_day_2026_2050'
nc_files <- list.files(file.path(dir_rawdata, 'ScaleUp_Global', folder), pattern='.nc$', full.names = T)
et_data_fut3 <- list()
for (i in 1:length(nc_files)) {
  # each file stores data from 2015 to 2100
  met_data <- terra::rast(nc_files[i])
  et_data_fut3 <- c(et_data_fut3, met_data)
}
et_data_fut3[[6]] <- c(et_data_fut3[[6]], et_data_fut3[[7]])

nmodels <- 6

# logics of conversion: 
# 1) get average et of each year; 2) calculate trend of 2015-2030 and that of 2031-2045; 
# 3) get difference of the two periods; 4) unit conversion: mm/s -> mm/hr 
et_data_year <- rast()
for (imodel in 1:nmodels) {
  data_north <- terra::mask(et_data_fut3[[imodel]], init(et_data_fut3[[imodel]], "y") > 0, maskvalues=0)
  data_south <- terra::mask(et_data_fut3[[imodel]], init(et_data_fut3[[imodel]], "y") <= 0, maskvalues=0)
  #
  for (iyear in 1:36) {
    r_north <- terra::app(data_north[[c(5:9) + (iyear-1)*12]], fun=mean, na.rm=T)
    r_south <- terra::app(data_south[[c(1:3, 11:12) + (iyear-1)*12]], fun=mean, na.rm=T)
    
    et_data_1year <- terra::cover(r_north, r_south)
    if (imodel > 1) {
      et_data_1year <- terra::resample(et_data_1year, et_data_year[[1]])
    }
    et_data_year <- c(et_data_year, et_data_1year)
    print(iyear)
  }
}

# average across different models
et_gs_2015_2045 <- terra::rast()
for (iyear in 1:36) {
  et_data_1year <- terra::app(et_data_year[[seq(iyear, 36*6, 36)]], fun=mean, na.rm=T)
  et_gs_2015_2045 <- c(et_gs_2015_2045, et_data_1year)
}
et_gs_2015_2045 <- rotate(et_gs_2015_2045)

et_at_model_past <- terra::app(et_gs_2015_2045[[1:15]], function(x) {ifelse(any(is.na(x)), NA, sens.slope(x)$estimate)})
et_at_model_fut <- terra::app(et_gs_2015_2045[[16:30]], function(x) {ifelse(any(is.na(x)), NA, sens.slope(x)$estimate)})

et_at_past_fut_change <- et_at_model_fut - et_at_model_past

et_at_past_fut_change <- et_at_past_fut_change * 3600   # convert unit from mm/s to mm/hr
et_at_past_fut_change <- terra::resample(et_at_past_fut_change, et_gs_past)

writeRaster(et_at_past_fut_change, filename = file.path('data', 'et_at_past_fut_change.tif'), overwrite=TRUE)

################################################################################
# 9. get current and future disturbance effects
################################################################################
trend_driver <- read.csv(file.path("data", "potential_trend_driver_sen.csv"))
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category")], by = "site_ID")
}

DIST_climate_IGBP_wet <- trend_driver %>% filter(category == 'wet') %>% 
  group_by(site_ID) %>% slice_max(year_mid) %>%
  filter(year_mid >= 2010) %>%
  arrange(Climate.class, IGBP) %>%
  group_by(Climate.class, IGBP) %>% 
  summarise(DIST_at = mean(DIST_at, na.rm = T), .groups = "drop")

DIST_climate_wet <- trend_driver %>% filter(category == 'wet') %>% 
  group_by(site_ID) %>% slice_max(year_mid) %>%
  group_by(Climate.class) %>% 
  summarise(DIST_at = mean(DIST_at, na.rm = T), .groups = "drop")


unique_climate_veg_groups <- levels(zone_group)[[1]]
cells <- which(terra::values(zone_group) == unique_climate_veg_groups$value[unique_climate_veg_groups$label=='wet'])
df_wet_DIST_at <- data.frame(icell = cells, 
                             climate_id = terra::values(koppen_climate_0.25)[cells], 
                             IGBP_id = terra::values(land_cover)[cells])
df_wet_DIST_at <- df_wet_DIST_at %>% left_join(df_koppen_climate, by = "climate_id") %>% 
  left_join(df_IGBP, by = "IGBP_id") 
  
Climate_IGBP_all <- df_wet_DIST_at[, c("Climate.class", "IGBP")] %>% distinct() %>% arrange(Climate.class, IGBP)
DIST_climate_IGBP_wet <- Climate_IGBP_all %>% left_join(DIST_climate_IGBP_wet, by = c("Climate.class", "IGBP")) %>% 
  left_join(DIST_climate_wet, by = "Climate.class", suffix = c("", "_fill")) %>%
  mutate(DIST_at = coalesce(DIST_at, DIST_at_fill)) %>%
  select(-DIST_at_fill)                                         

df_wet_DIST_at <- df_wet_DIST_at %>% left_join(DIST_climate_IGBP_wet, by = c("Climate.class", "IGBP"))

DIST_at_past <- zone_group
terra::values(DIST_at_past) <- NA
terra::values(DIST_at_past)[df_wet_DIST_at$icell] <- df_wet_DIST_at$DIST_at
plot(DIST_at_past)

writeRaster(DIST_at_past, filename = file.path('data', 'DIST_at_wet_past.tif'), overwrite=TRUE)

################################################################################
# 10. get current and future CO2 concentration
################################################################################
 
# https://doi.org/10.5281/zenodo.5021361   1-degree resolution
co2_1850_2013 <- terra::rast(file.path(dir_rawdata, 'ScaleUp_Global', 'co2', 'CO2_1deg_month_1850-2013.nc'))
co2_2015_2150 <- terra::rast(file.path(dir_rawdata, 'ScaleUp_Global', 'co2', 'CO2_SSP245_2015_2150.nc'))
terra::ext(co2_1850_2013) <- terra::ext(co2_2015_2150)
#
co2_baseline1 <- terra::app(co2_1850_2013[[1873:1968]], fun = mean, na.rm=T)  # 2005-2013
co2_baseline2 <- terra::app(co2_2015_2150[[1:96]], fun = mean, na.rm=T)       # 2015-2023

co2_baseline1 <- terra::resample(co2_baseline1, gpp_baseline)
co2_baseline2 <- terra::resample(co2_baseline2, gpp_baseline)

co2_baseline <- (co2_baseline1 + co2_baseline2) / 2
# plot(co2_baseline)

co2_fut <- terra::app(co2_2015_2150[[121:(120+12*25)]], fun = mean, na.rm=T)  # 2025-2050
co2_fut <- terra::resample(co2_fut, gpp_baseline)
# plot(co2_fut - co2_baseline)

writeRaster(co2_baseline, filename = file.path('data', 'co2_mean_2005-2023.tif'), overwrite=TRUE)
writeRaster(co2_fut, filename = file.path('data', 'co2_mean_2025-2050.tif'), overwrite=TRUE)

# I need to ensure that all the data used for projection are at the same resolution and the same extent. 
