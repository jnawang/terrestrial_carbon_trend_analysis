# This script prepares data for projecting future NEE trends.
# Junna Wang, 03-17-2026
# 
# make the unit of different variables clear enough
# TA: C; P: mm/d; SW: W/m2; LE: W/m2; ET: mm/hr, and then should be converted to LE unit. 


# global monthly CO2 is obtained from https://zenodo.org/records/5021361
# citation: Cheng, W., Dan, L., Deng, X., Feng, J., Wang, Y., Peng, J., Tian, J., Qi, W., Liu, Z., Zheng, X., Zhou, D., Jiang, S., Zhao, H., & Wang, X. (2021). Global monthly distributions of atmospheric CO2 concentrations under the historical and future scenarios [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5021361

# for disturbance effects: using calculated values based on global MODIS LAI, and we need SOS and EOS map for this purpose

rm(list=ls())
# change this location where the raw data will be saved accordingly
dir_rawdata <- '/Volumes/WZZ_disk'
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'
#
library(librarian)
shelf(tidyverse, ggpubr, terra, ncdf4, future.apply, trend, mgcv)

###########################################
# 1. get start and end of growing season
###########################################
source(file.path("scripts2026", "00_functions.R"))

data_grow_season <- read.csv(file.path("data2026", "growing_season_yearly.csv"))
sites_long <- read.csv(file.path("data2026", "long_term_ec_sites.csv"))

mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
year_round <- c("Af", "Am", "Cfa", "Cwa")

data_grow_season$sos_year <- rowMeans(
  data_grow_season[, c("sos_TRS", "sos_DER")], na.rm = TRUE
)
data_grow_season$eos_year <- rowMeans(
  data_grow_season[, c("eos_TRS", "eos_DER")], na.rm = TRUE
)
data_grow_season$sos_year[is.nan(data_grow_season$sos_year)] <- NA
data_grow_season$eos_year[is.nan(data_grow_season$eos_year)] <- NA

sos_eos_site <- data_grow_season %>%
  group_by(site_ID) %>%
  summarise(
    sos = mean_na(sos_year),
    eos = mean_na(eos_year),
    .groups = "drop"
  ) %>%
  left_join(
    sites_long %>%
      select(site_ID, Climate.class) %>%
      distinct(),
    by = "site_ID"
  ) %>%
  mutate(
    climate_parent = substr(Climate.class, 1, 2),
    climate_grandparent = substr(Climate.class, 1, 1)
  ) %>%
  filter(!is.na(Climate.class), is.finite(sos), is.finite(eos))

sos_eos_exact <- sos_eos_site %>%
  group_by(Climate.class) %>%
  summarise(
    nsite_exact = n_distinct(site_ID),
    sos_exact = mean(sos),
    eos_exact = mean(eos),
    .groups = "drop"
  )

sos_eos_parent <- sos_eos_site %>%
  group_by(climate_parent) %>%
  summarise(
    nsite_parent = n_distinct(site_ID),
    sos_parent = mean(sos),
    eos_parent = mean(eos),
    .groups = "drop"
  )

sos_eos_grandparent <- sos_eos_site %>%
  group_by(climate_grandparent) %>%
  summarise(
    nsite_grandparent = n_distinct(site_ID),
    sos_grandparent = mean(sos),
    eos_grandparent = mean(eos),
    .groups = "drop"
  )

sos_eos_group <- df_koppen_climate %>%
  mutate(
    climate_parent = substr(Climate.class, 1, 2),
    climate_grandparent = substr(Climate.class, 1, 1)
  ) %>%
  left_join(sos_eos_exact, by = "Climate.class") %>%
  left_join(sos_eos_parent, by = "climate_parent") %>%
  left_join(sos_eos_grandparent, by = "climate_grandparent") %>%
  mutate(
    assignment_source = case_when(
      Climate.class %in% year_round ~ "year_round",
      nsite_exact >= 2 ~ "exact_class",
      nchar(Climate.class) > 2 & !is.na(sos_parent) ~
        paste0("parent_", climate_parent),
      !is.na(sos_grandparent) ~
        paste0("grandparent_", climate_grandparent),
      TRUE ~ "missing"
    ),
    sos_group = case_when(
      Climate.class %in% year_round ~ 1,
      nsite_exact >= 2 ~ sos_exact,
      nchar(Climate.class) > 2 & !is.na(sos_parent) ~ sos_parent,
      !is.na(sos_grandparent) ~ sos_grandparent,
      TRUE ~ NA_real_
    ),
    eos_group = case_when(
      Climate.class %in% year_round ~ 366,
      nsite_exact >= 2 ~ eos_exact,
      nchar(Climate.class) > 2 & !is.na(eos_parent) ~ eos_parent,
      !is.na(eos_grandparent) ~ eos_grandparent,
      TRUE ~ NA_real_
    ),
    sos_group = pmin(pmax(round(sos_group), 1), 366),
    eos_group = pmin(pmax(round(eos_group), 1), 366),
    nsite_exact = replace_na(nsite_exact, 0L),
    nsite_parent = replace_na(nsite_parent, 0L),
    nsite_grandparent = replace_na(nsite_grandparent, 0L)
  ) %>%
  select(
    Climate.class, climate_id,
    # climate_parent, climate_grandparent,
    sos_group, eos_group, 
    # assignment_source, nsite_exact, nsite_parent, nsite_grandparent
  ) %>%
  arrange(climate_id)

if (anyNA(sos_eos_group$sos_group) ||
    anyNA(sos_eos_group$eos_group)) {
  stop("Some climate classes still lack SOS/EOS information.")
}

# print(sos_eos_group, n = Inf)


# -------------------- create SOS/EOS raster for Google Earth Engine --------------------
koppen_climate <- terra::rast(file.path('data2026', 'koppen_geiger_0p00833333.tif'))
land_cover <- terra::rast(file.path('data2026', '968_Land_Cover_Class_0.25degree_layer2.tif'))
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

SOS <- subst(
  koppen_climate_0.25,
  from = sos_eos_group$climate_id,
  to = sos_eos_group$sos_group,
  others = NA
)

EOS <- subst(
  koppen_climate_0.25,
  from = sos_eos_group$climate_id,
  to = sos_eos_group$eos_group,
  others = NA
)

names(SOS) <- "SOS"
names(EOS) <- "EOS"

SOS_EOS <- c(SOS, EOS)
# 
# writeRaster(
#   SOS_EOS,
#   file.path("data2026", "global_SOS_EOS_025degree.tif"),
#   datatype = "INT2S",
#   NAflag = -9999,
#   overwrite = TRUE
# )

plot(SOS, main = "Start of growing season")
plot(EOS, main = "End of growing season")

###########################################################
# 2. get global zonation based on climate and vegetation
###########################################################

# aridity index
AI_0.25 <- terra::rast(file.path('data2026', 'Global-AI.tif'))
AI_0.25 <- terra::resample(AI_0.25, land_cover) / 10000       # values are enlarged by 10000
# Zomer, R.J., Xu, J. & Trabucco, A. Version 3 of the Global Aridity Index and Potential Evapotranspiration Database. Sci Data 9, 409 (2022). https://doi.org/10.1038/s41597-022-01493-1

# mean annual temperature
MAT <- terra::rast(file.path('data2026', 'baseline_mat.tif'))
MAT <- terra::resample(MAT, land_cover)
#
# a data frame to store zone info
df_zone <- data.frame(koppen_climate=as.integer(terra::values(koppen_climate_0.25)),
                      land_cover=as.integer(terra::values(land_cover)),
                      AI = terra::values(AI_0.25),
                      MAT = terra::values(MAT))
names(df_zone)[3:4] <- c("AI", "MAT")

# using the same rules as we classify study sites
# Step 1: classify climate into: dry, wet-warm (wet), wet-cold (cold)
# rule 1: all arctic and subarctic climates ("Dfc", "Dfd", "Dsc", "Dsd", "Dwc", "Dwd", "ET") are cold
# rule 2: for the rest grid-cells, all grids with AI <= 0.5, as dry climate
# rule 3: all grids with AI > 0.5, as wet climate, and then wet climates are divided into cold and warm
# rule 4: climates in tropical and subtropical ("Am", "Af", "Aw", "Cfa", "Cwa", "Cwb", "Cwc", "Csa") are warm 
# rule 5: climates in continental climates ("Dfa", "Dfb", "Dsa", "Dsb", "Dwa", "Dwb") are cold
# rule 6: sites in temperate climate dependent vegetation type: evergreen broadleaf forest EBF is warm, and other vegetation types are cold
# rule 7: sites in dry climate (arid, semi-arid, mediterrean climate, and tropical dry climate): MAT > 10C, warm; MAT < 10C, cold. 

# Step 2: classify vegetation into ENF, DBF&MF (DBF), nonforest
df_zone <- df_zone %>% mutate(group_climate = case_when(
  koppen_climate %in% c(19:20, 23:24, 27:30) ~ "cold",
  AI < 0.5 & koppen_climate %in% (1:30) & land_cover != 11 ~ "dry",                   # ensure not in the ocean and not wetlands
  AI < 0.5 & koppen_climate %in% (1:30) & land_cover == 11 & MAT >= 10 ~ "wet",       # wetlands in warm zone are wet-warm ecosystems
  AI < 0.5 & koppen_climate %in% (1:30) & land_cover == 11 & MAT < 10 ~ "cold",       # wetlands in cold zone are wet-cold ecosystems
  AI >= 0.5 & koppen_climate %in% c(17:18, 21:22, 25:26) ~ "cold", 
  AI >= 0.5 & koppen_climate %in% c(1:2, 11:14) ~ "wet", 
  AI >= 0.5 & koppen_climate %in% c(15:16) & land_cover %in% c(2) ~ "wet", 
  AI >= 0.5 & koppen_climate %in% c(15:16) & !land_cover %in% c(2) ~ "cold", 
  AI >= 0.5 & koppen_climate %in% c(3:10) & MAT >= 10 ~ "wet", 
  AI >= 0.5 & koppen_climate %in% c(3:10) & MAT < 10 ~ "cold")) %>%
  mutate(group = case_when(group_climate == "wet" ~ "wet",
                           group_climate == "dry" ~ "dry",
                           group_climate == "cold" & land_cover %in% c(6:11, 14) ~ "cold_nonforest",
                           group_climate == "cold" & land_cover %in% c(1:2) ~ "cold_ENF",
                           group_climate == "cold" & land_cover %in% c(3:5) ~ "cold_DBF"))

zone_grow_season <- koppen_climate_0.25

# zone for climate-vegetation groups
zone_group <- koppen_climate_0.25
terra::values(zone_group) <- NA
terra::values(zone_group) <- df_zone$group
plot(zone_group)
# why there are so many missing points? these are croplands, urban, ice, and barren.  

writeRaster(zone_group, filename = file.path('data2026', 'zones_climate_vegetation_group.tif'), overwrite=TRUE)

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

}
writeRaster(t2m_gs, filename = file.path('data2026', 't2m_gs_2005_2023.tif'), overwrite=TRUE)
writeRaster(t2m_ngs, filename = file.path('data2026', 't2m_ngs_2005_2023.tif'), overwrite=TRUE)
terra::tmpFiles(remove = TRUE)
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
  
}
writeRaster(tp_gs, filename = file.path('data2026', 'tp_gs_2005_2023.tif'), overwrite=TRUE)
writeRaster(tp_ngs, filename = file.path('data2026', 'tp_ngs_2005_2023.tif'), overwrite=TRUE)
terra::tmpFiles(remove = TRUE)
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
  
}
writeRaster(ssrd_gs, filename = file.path('data2026', 'ssrd_gs_2005_2023.tif'), overwrite=TRUE)
terra::tmpFiles(remove = TRUE)
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
}
rm(met_data,test)
gc()
terra::tmpFiles(current=FALSE,orphan=TRUE,old=FALSE,remove=TRUE)
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
  writeRaster(tp_gs_fut, filename = paste0('data2026/', var, '_gs_2026_2050.tif'), overwrite=TRUE)
  
  
  # non-growing season
  if (var %in% c('pr', 'tas')) {
    tp_ngs_fut <- rast()
    for (i in 1:length(tp_data_fut)) {
      test <- process_one_year(tp_data_fut[[i]], zone_masks_north, zone_masks_south, sos_eos_group, gs=FALSE)
      tp_ngs_fut <- c(tp_ngs_fut, test)
      print(i)
    }
    writeRaster(tp_ngs_fut, filename = paste0('data2026/', var, '_ngs_2026_2050.tif'), overwrite=TRUE)
  }
  
  rm(tp_gs_fut,test)
  gc()
  terra::tmpFiles(current=FALSE,orphan=TRUE,old=FALSE,remove=TRUE)  
}

#########################################################################################################
# 7. get baseline (2005~2021) information about NEE, GPP, ET, MAT, MAP, climate and vegetation group 
#########################################################################################################
# folder <- 'gpp_fluxcom_2005_2021'
# nc_files <- list.files(file.path(dir_rawdata, 'ScaleUp_Global', folder), pattern='.nc$', full.names = T)
# gpp_baseline <- terra::rast()
# for (i in 1:length(nc_files)) {
#   print(i)
#   gpp_data <- terra::rast(nc_files[i])
#   gpp_data <- terra::app(gpp_data, fun=sum, na.rm=T)
#   gpp_baseline <- c(gpp_baseline, gpp_data)
# }
# gpp_baseline <- terra::app(gpp_baseline, fun=mean, na.rm=T)
# 
# writeRaster(gpp_baseline, filename = file.path('data2026', 'baseline_gpp_2005-2021.tif'), overwrite=TRUE)

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

writeRaster(nee_baseline, filename = file.path('data2026', 'baseline_nee_2005-2021.tif'), overwrite=TRUE) 

#-----------
# folder <- 'et_day_2005_2021'
# nc_files <- list.files(file.path(dir_rawdata, 'ScaleUp_Global', folder), pattern='.nc$', full.names = T)
# et_baseline <- terra::rast()
# for (i in 1:length(nc_files)) {
#   print(i)
#   et_data <- terra::rast(nc_files[i])
#   et_data <- terra::app(et_data, fun=mean, na.rm=T)
#   et_baseline <- c(et_baseline, et_data)
# }
# et_baseline <- terra::app(et_baseline, fun=mean, na.rm=T)
# 
# writeRaster(et_baseline, filename = file.path('data2026', 'baseline_et_2005-2021.tif'), overwrite=TRUE) 

#--------MAT, MAP
# mat_baseline <- terra::rast(file.path('data2026', 'MAT.tif'))
# map_baseline <- terra::rast(file.path('data2026', 'MAP.tif'))
# 
# mat_baseline <- terra::resample(mat_baseline, gpp_baseline) 
# map_baseline <- terra::resample(map_baseline, gpp_baseline)
# 
# writeRaster(mat_baseline, filename = file.path('data2026', 'baseline_mat.tif'), overwrite=TRUE)
# writeRaster(map_baseline, filename = file.path('data2026', 'baseline_map.tif'), overwrite=TRUE)

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

writeRaster(et_at_past_fut_change, filename = file.path('data2026', 'et_at_past_fut_change.tif'), overwrite=TRUE)

################################################################################
# 9. get current and future disturbance effects
################################################################################
data_dir <- file.path(dir_rawdata, "Carbon_trend/SiteData/", "LAI_4day_DBF")
zone_file <- file.path("data2026", "zones_climate_vegetation_group.tif")
sos_eos_file <- file.path("data2026", "global_SOS_EOS_025degree.tif")

files <- list.files(data_dir,pattern="^MODIS_LAI_QC_1000points_[0-9]{4}\\.csv$",full.names=TRUE)
if(!length(files)) stop("No GEE CSV files found.")
files <- files[1:19]       # get the first 19 lines, just use 2005-2023

read_GEE <- function(f) readr::read_csv(f,show_col_types=FALSE) %>%
  dplyr::select(sample_ID,sample_type,latitude,longitude,zone,
                tidyselect::matches("^\\d{4}_\\d{2}_\\d{2}_")) %>%
  tidyr::pivot_longer(
    tidyselect::matches("^\\d{4}_\\d{2}_\\d{2}_"),
    names_to=c("Date",".value"),
    names_pattern="^(\\d{4}_\\d{2}_\\d{2})_(.*)$"
  ) %>%
  mutate(Date=lubridate::ymd(str_replace_all(Date,"_","-")))

data_LAI <- map_dfr(files,read_GEE) %>%
  filter(sample_type=="random") %>%
  mutate(
    across(c(LAI,MODLAND,CloudState,Aerosol,Cloud_Shadow,Internal_CloudMask),as.numeric),
    LAI=if_else(LAI>20,NA_real_,LAI),
    year=lubridate::year(Date),
    doy=case_when(
      latitude < -10 & lubridate::leap_year(year) & lubridate::month(Date)<7 ~ lubridate::yday(Date)-182+366,
      latitude < -10 & lubridate::leap_year(year) & lubridate::month(Date)>=7 ~ lubridate::yday(Date)-182,
      latitude < -10 & !lubridate::leap_year(year) & lubridate::month(Date)<7 ~ lubridate::yday(Date)-181+365,
      latitude < -10 & !lubridate::leap_year(year) & lubridate::month(Date)>=7 ~ lubridate::yday(Date)-181,
      TRUE ~ lubridate::yday(Date)
    )
  ) %>%
  filter(is.finite(longitude),is.finite(latitude),
         between(longitude,-180,180),between(latitude,-90,90)) %>%
  arrange(sample_ID,Date)

sos_eos <- rast(sos_eos_file)
if(nlyr(sos_eos)<2) stop("SOS/EOS raster must contain two layers.")
names(sos_eos)[1:2] <- c("SOS","EOS")

point_info <- data_LAI %>%
  distinct(sample_ID,latitude,longitude)

point_vect <- vect(
  as.matrix(point_info[,c("longitude","latitude")]),
  type="points",
  crs="+proj=longlat +datum=WGS84"
)

sos_eos_values <- as.data.frame(
  terra::extract(sos_eos,point_vect,ID=FALSE)
)

point_info <- bind_cols(point_info,sos_eos_values) %>%
  transmute(sample_ID,latitude,longitude,SOS=round(SOS),EOS=round(EOS))

data_LAI <- data_LAI %>%
  left_join(point_info,by=c("sample_ID","latitude","longitude"))

calculate_deficit <- function(x) {
  x <- x %>% arrange(Date)
  sos <- first(na.omit(x$SOS))
  eos <- first(na.omit(x$EOS))
  
  if(!length(sos) || !length(eos) ||
     !is.finite(sos) || !is.finite(eos) || eos<sos) return(NULL)
  
  x <- x %>%
    mutate(
      strict_qc=CloudState==0 & Aerosol==0 &
        Cloud_Shadow==0 & Internal_CloudMask==0,
      weight=case_when(
        strict_qc & MODLAND==0 ~ 1,
        strict_qc & !is.na(MODLAND) & MODLAND!=0 ~ 0.5,
        TRUE ~ 0
      )
    )
  
  seasonal <- x %>%
    filter(strict_qc,is.finite(LAI)) %>%
    group_by(doy) %>%
    summarise(LAI=mean(LAI,na.rm=TRUE),.groups="drop")
  
  if(nrow(seasonal)<10 || n_distinct(seasonal$doy)<10) return(NULL)
  
  gam_model <- tryCatch(
    mgcv::gam(LAI~s(doy),data=seasonal),
    error=function(e) NULL
  )
  
  if(is.null(gam_model)) return(NULL)
  
  seasonal_smooth <- tibble(doy=sos:eos)
  pred <- predict(gam_model,newdata=seasonal_smooth,se.fit=TRUE)
  seasonal_smooth$LAI_gam <- as.numeric(pred$fit-1.96*pred$se.fit)
  
  x <- x %>%
    left_join(seasonal_smooth,by="doy")
  
  good_years <- x %>%
    group_by(year) %>%
    summarise(f=sum(weight==1)/n(),.groups="drop")
  
  mean_good <- mean(good_years$f,na.rm=TRUE)
  
  bad_years <- good_years %>%
    filter(is.finite(mean_good),f<mean_good*0.5) %>%
    pull(year)
  
  x <- x %>%
    filter(!year %in% bad_years)
  
  if(!nrow(x)) return(NULL)
  
  x$LAI_intp <- x$LAI
  x$LAI_intp[x$weight==0] <- NA_real_
  
  reject <- x$weight==0.5 &
    is.finite(x$LAI_intp) &
    is.finite(x$LAI_gam) &
    x$LAI_intp<0.5*x$LAI_gam
  
  x$LAI_intp[reject] <- NA_real_
  ok <- is.finite(x$LAI_intp) & !is.na(x$Date)
  
  if(sum(ok)<2) return(NULL)
  
  x$LAI_intp <- approx(
    x=as.numeric(x$Date[ok]),
    y=x$LAI_intp[ok],
    xout=as.numeric(x$Date),
    method="linear",
    rule=2
  )$y
  
  yearly_LAI <- x %>%
    filter(between(doy,sos,eos)) %>%
    group_by(year) %>%
    summarise(LAI_gs=mean(LAI_intp,na.rm=TRUE),.groups="drop") %>%
    filter(is.finite(LAI_gs))
  
  x$LAI_intp_detrend <- x$LAI_intp
  
  if(nrow(yearly_LAI)>=3 && n_distinct(yearly_LAI$year)>=3) {
    lm_model <- lm(LAI_gs~year,data=yearly_LAI)
    trend_p <- summary(lm_model)$coefficients["year","Pr(>|t|)"]
    
    if(is.finite(trend_p) && trend_p<0.05) {
      yearly_LAI$annual_trend <- fitted(lm_model)
      
      x <- x %>%
        left_join(
          yearly_LAI %>% dplyr::select(year,annual_trend),
          by="year"
        )
      
      x$LAI_intp_detrend <- x$LAI_intp-x$annual_trend+
        mean(x$annual_trend,na.rm=TRUE)
    }
  }
  
  LAI_gs <- mean(
    x$LAI_gam[between(x$doy,sos,eos)],
    na.rm=TRUE
  )*(eos-sos+1)
  
  if(!is.finite(LAI_gs) || LAI_gs<=0) return(NULL)
  
  x %>%
    filter(between(doy,sos,eos)) %>%
    mutate(deficit_obs=pmax(LAI_gam-LAI_intp_detrend,0)) %>%
    group_by(sample_ID,latitude,longitude,year) %>%
    summarise(
      deficit_obs=sum(deficit_obs,na.rm=TRUE)*4/LAI_gs,
      .groups="drop"
    )
}

data_split <- split(data_LAI,data_LAI$sample_ID)

annual_deficit <- map_dfr(
  data_split,
  ~tryCatch(calculate_deficit(.x),error=function(e) NULL)
)

if(!nrow(annual_deficit))
  stop("No annual LAI-deficit estimates were produced.")

source(file.path("scripts2026","00_functions.R"))

point_trends <- annual_deficit %>%
  group_by(sample_ID,latitude,longitude) %>%
  summarise(
    n_year=n_distinct(year[is.finite(deficit_obs)]),
    DIST_at=if(n_year>=8)
      cal_sens_slope(year,deficit_obs)
    else NA_real_,
    .groups="drop"
  ) %>%
  filter(n_year>=8,is.finite(DIST_at))

print(nrow(point_trends))
summary(point_trends$DIST_at)

if(nrow(point_trends)<20)
  stop("Only ",nrow(point_trends),
       " valid point trends were obtained; at least 20 are required.")

zone_raster <- rast(zone_file)
zone_mask <- zone_raster
values(zone_mask) <- ifelse(values(zone_raster,mat=FALSE)==1,1,NA_real_)

zone_df <- as.data.frame(zone_mask,xy=TRUE,na.rm=TRUE)
names(zone_df)[3] <- "zone"

trend_limit <- max(abs(point_trends$DIST_at),na.rm=TRUE)
if(!is.finite(trend_limit) || trend_limit==0) trend_limit <- 1

p_point_trends <- ggplot() +
  geom_raster(data=zone_df,aes(x=x,y=y),fill="grey85") +
  geom_point(
    data=point_trends,
    aes(x=longitude,y=latitude,color=DIST_at),
    size=1.8,
    alpha=0.9
  ) +
  scale_color_gradient2(
    low="#2166AC",
    mid="white",
    high="#B2182B",
    midpoint=0,
    limits=c(-trend_limit,trend_limit),
    name=expression("LAI deficit trend (yr"^{-1}*")")
  ) +
  coord_equal(
    xlim=c(-180,180),
    ylim=c(-60,85),
    expand=FALSE
  ) +
  labs(
    x="Longitude",
    y="Latitude",
    title=paste0(
      "LAI-deficit trends at ",
      nrow(point_trends),
      " sample points"
    )
  ) +
  theme_classic()

# print(p_point_trends)

lonlat_to_xyz <- function(longitude,latitude) {
  r <- 6371008.8
  lon <- as.numeric(longitude)*pi/180
  lat <- as.numeric(latitude)*pi/180
  
  m <- cbind(
    x=r*cos(lat)*cos(lon),
    y=r*cos(lat)*sin(lon),
    z=r*sin(lat)
  )
  
  storage.mode(m) <- "double"
  m
}

point_xyz <- lonlat_to_xyz(
  point_trends$longitude,
  point_trends$latitude
)

point_trends$x <- point_xyz[,1]
point_trends$y <- point_xyz[,2]
point_trends$z <- point_xyz[,3]

xyz_matrix <- function(x) {
  m <- cbind(
    x=as.numeric(x$x),
    y=as.numeric(x$y),
    z=as.numeric(x$z)
  )
  
  storage.mode(m) <- "double"
  m
}

idw_predict <- function(train_xyz,train_value,new_xyz,k,power,
                        max_distance=1500000) {
  train_xyz <- matrix(
    as.numeric(as.matrix(train_xyz)),
    ncol=3
  )
  
  new_xyz <- matrix(
    as.numeric(as.matrix(new_xyz)),
    ncol=3
  )
  
  train_value <- as.numeric(train_value)
  
  if(!nrow(train_xyz) || !nrow(new_xyz))
    stop("Empty coordinate matrix.")
  
  if(any(!is.finite(train_xyz)) || any(!is.finite(new_xyz)))
    stop("Coordinate matrices contain non-finite values.")
  
  if(length(train_value)!=nrow(train_xyz) ||
     any(!is.finite(train_value)))
    stop("Invalid IDW training values.")
  
  r <- 6371008.8
  k <- min(as.integer(k),nrow(train_xyz))
  
  prediction <- rep(NA_real_,nrow(new_xyz))
  nearest_distance_km <- rep(NA_real_,nrow(new_xyz))
  n_neighbors <- integer(nrow(new_xyz))
  
  for(i in seq_len(nrow(new_xyz))) {
    difference <- sweep(
      train_xyz,
      MARGIN=2,
      STATS=new_xyz[i,],
      FUN="-"
    )
    
    chord <- sqrt(rowSums(difference^2))
    distance <- 2*r*asin(pmin(chord/(2*r),1))
    index <- order(distance)[seq_len(k)]
    distance_k <- distance[index]
    valid <- distance_k<=max_distance
    
    nearest_distance_km[i] <- distance_k[1]/1000
    n_neighbors[i] <- sum(valid)
    
    if(any(valid)) {
      weights <- pmax(distance_k[valid],1)^(-power)
      
      prediction[i] <- weighted.mean(
        train_value[index[valid]],
        weights
      )
    }
  }
  
  tibble(
    DIST_at=prediction,
    nearest_distance_km=nearest_distance_km,
    n_neighbors=n_neighbors
  )
}

set.seed(123)

point_trends <- point_trends %>%
  mutate(
    block=paste(
      floor((longitude+180)/10),
      floor((latitude+90)/10),
      sep="_"
    )
  )

blocks <- sample(unique(point_trends$block))

if(length(blocks)<5)
  stop("Fewer than five spatial blocks are available.")

block_folds <- tibble(
  block=blocks,
  fold=rep(1:5,length.out=length(blocks))
)

point_trends <- point_trends %>%
  left_join(block_folds,by="block")

parameters <- expand_grid(
  k=c(5,10,15,20),
  power=c(1,1.5,2,3)
)

cv_results <- pmap_dfr(parameters,function(k,power) {
  fold_results <- map_dfr(1:5,function(f) {
    train <- point_trends %>% filter(fold!=f)
    test <- point_trends %>% filter(fold==f)
    
    if(!nrow(train) || !nrow(test)) return(NULL)
    
    pred <- idw_predict(
      xyz_matrix(train),
      train$DIST_at,
      xyz_matrix(test),
      k,
      power
    )
    
    tibble(
      obs=test$DIST_at,
      pred=pred$DIST_at
    ) %>%
      filter(is.finite(obs),is.finite(pred))
  })
  
  tibble(
    k=k,
    power=power,
    RMSE=if(nrow(fold_results))
      sqrt(mean((fold_results$obs-fold_results$pred)^2))
    else Inf
  )
})

best <- cv_results %>%
  filter(is.finite(RMSE)) %>%
  slice_min(RMSE,n=1,with_ties=FALSE)

if(!nrow(best))
  stop("No IDW parameter combination produced valid predictions.")

print(best)

zone_cells <- which(values(zone_raster,mat=FALSE)==1)
zone_xy <- xyFromCell(zone_raster,zone_cells)
zone_xyz <- lonlat_to_xyz(zone_xy[,1],zone_xy[,2])

interpolated <- idw_predict(
  xyz_matrix(point_trends),
  point_trends$DIST_at,
  zone_xyz,
  best$k[[1]],
  best$power[[1]]
)

output_raster <- zone_raster
values(output_raster) <- NA_real_
values(output_raster)[zone_cells] <- interpolated$DIST_at
names(output_raster) <- "DIST_at"

terra::writeRaster(
  output_raster,
  file.path(
    "data2026",
    "LAI_deficit_trends_DBF.tif"
  ),
  overwrite=TRUE,
  datatype="FLT4S",
  NAflag=-9999
)

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

co2_baseline1 <- terra::resample(co2_baseline1, nee_baseline)
co2_baseline2 <- terra::resample(co2_baseline2, nee_baseline)

co2_baseline <- (co2_baseline1 + co2_baseline2) / 2
# plot(co2_baseline)

co2_fut <- terra::app(co2_2015_2150[[121:(120+12*20)]], fun = mean, na.rm=T)  # 2025-2045
co2_fut <- terra::resample(co2_fut, nee_baseline)
# plot(co2_fut - co2_baseline)

writeRaster(co2_baseline, filename = file.path('data2026', 'co2_mean_2005-2023.tif'), overwrite=TRUE)
writeRaster(co2_fut, filename = file.path('data2026', 'co2_mean_2025-2045.tif'), overwrite=TRUE)

################################################################################
# 11. get global GPP trends in 2001-2020
################################################################################
files <- list.files(
  file.path(dir_rawdata, "CEDAR-GPP", "ST_CFE-ML_DT"),
  pattern = "\\.nc$|\\.tif$",
  full.names = TRUE
)

gpp_monthly <- rast(files, subds = "GPP_mean")

terra::ext(gpp_monthly) <- terra::ext(-180, 180, -90, 90)
terra::crs(gpp_monthly) <- "EPSG:4326"

##
start_year <- 2001

dates <- seq(
  from = as.Date(paste0(start_year, "-01-01")),
  by = "month",
  length.out = nlyr(gpp_monthly)
)

days_month <- lubridate::days_in_month(dates)
year_id <- lubridate::year(dates)

# convert daily GPP to monthly total GPP
gpp_monthly_total <- gpp_monthly * days_month

# If monthly GPP is already monthly total, use this:
gpp_annual <- terra::tapp(
  gpp_monthly_total,
  index = year_id,
  fun = sum,
  na.rm = TRUE
)

# compare with we calculated trends
# trend_driver <- read.csv(file.path("data2026", "potential_trend_driver.csv"))
# sites_long <- read.csv(file.path('data2026', 'long_term_ec_sites.csv'))
# 
# sites_long$map_GPP_at <- NA
# sites_xy <- terra::vect(
#   sites_long,
#   geom = c("Long", "Lat"),
#   crs = "EPSG:4326"
# )
# sites_long_gpp <- terra::extract(gpp_annual, sites_xy)
# for (i in 1:nrow(sites_long)) {
#   if (sites_long_gpp[i, 2] < 0) { next }
#   sites_long$map_GPP_at[i] <- sens.slope(t(sites_long_gpp[i, 2:21]))$estimate / 100
# }
# 
# trend_driver_one <- trend_driver %>%
#   dplyr::distinct(site_ID, .keep_all = TRUE) %>%
#   dplyr::select(site_ID, GPP_at_mean)
# 
# sites_long <- sites_long %>%
#   dplyr::left_join(trend_driver_one, by = "site_ID")
# 
# plot(sites_long$GPP_at_mean, sites_long$map_GPP_at)
# summary(lm(data = sites_long, map_GPP_at ~ GPP_at_mean))

gpp_annual[gpp_annual <= -9999] <- NA

# resolution from 0.05 to 0.25 degree
gpp_annual <- terra::resample(gpp_annual, zone_group)

# Sen's slope for annual GPP
sen_slope_fun <- function(y) {
  if (sum(!is.na(y)) < 3) return(NA_real_)
  
  out <- trend::sens.slope(y)
  
  as.numeric(out$estimates)
}

gpp_sen_slope <- terra::app(gpp_annual, fun = sen_slope_fun) / 100

names(gpp_sen_slope) <- "GPP_Sen_slope_per_year"

terra::writeRaster(
  gpp_sen_slope,
  file.path("data2026", "CEDAR_GPP_Sen_slope_ML.tif"),
  overwrite = TRUE
)