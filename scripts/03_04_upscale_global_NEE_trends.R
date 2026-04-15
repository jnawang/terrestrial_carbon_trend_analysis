# This script projects for future NEE trends.
# Junna Wang, 04-07-2026

# Be mindful that each data source has its own unit. 
# 

rm(list=ls())
library(librarian)
shelf(tidyverse, ggpubr, terra, future.apply, trend, mblm, ranger)

# use df_koppen_climate and df_IGBP
source(file.path("scripts", "00_functions.R"))

#--------------------------read in baseline information-----------------------------------
nee_baseline <- terra::rast(file.path("data", "baseline_nee_2005-2021.tif"))
gpp_baseline <- terra::rast(file.path("data", "baseline_gpp_2005-2021.tif"))
et_baseline <- terra::rast(file.path("data", "baseline_et_2005-2021.tif"))
# convert et (mm/hr) to le (W/m2)
le_baseline <- et_baseline / 3600 * 2.45 * 1000000
mat_baseline <- terra::rast(file.path("data", "baseline_mat.tif"))
map_baseline <- terra::rast(file.path("data", "baseline_map.tif"))

land_cover <- terra::rast(file.path('data', '968_Land_Cover_Class_0.25degree_layer2.tif'))
koppen_climate <- terra::rast(file.path('data', 'koppen_geiger_0p00833333.tif'))
aridity <- terra::rast(file.path('data', 'Global-AI.tif'))
zone_group <- terra::rast(file.path('data', 'zones_climate_vegetation_group.tif'))
#
koppen_climate <- terra::resample(koppen_climate, land_cover, method="near")
aridity <- terra::resample(aridity, land_cover) / 10000       # values are enlarged by 10000

co2_past <- terra::rast(file.path('data', 'co2_mean_2005-2023.tif'))
co2_fut <- terra::rast(file.path('data', 'co2_mean_2025-2050.tif'))

#---------------current and future climate trend--------------------------
ta_gs_past <- terra::rast(file.path('data', 't2m_gs_2005_2023.tif')) - 273.15      # convert Kelvin to Celcius
ta_ngs_past <- terra::rast(file.path('data', 't2m_ngs_2005_2023.tif')) - 273.15
ta_gs_fut <- terra::rast(file.path('data', 'tas_gs_2026_2050.tif')) - 273.15
ta_ngs_fut <- terra::rast(file.path('data', 'tas_ngs_2026_2050.tif')) - 273.15
#
# precipitation unit conversion
tp_gs_past <- terra::rast(file.path('data', 'tp_gs_2005_2023.tif')) * 1000 * 24      # convert unit from m/hr to mm/day
tp_ngs_past <- terra::rast(file.path('data', 'tp_ngs_2005_2023.tif')) * 1000 * 24    # convert unit from m/hr to mm/day 
tp_gs_fut <- terra::rast(file.path('data', 'pr_gs_2026_2050.tif')) * 86400           # convert mm/s to mm/day
tp_ngs_fut <- terra::rast(file.path('data', 'pr_ngs_2026_2050.tif')) * 86400         # convert mm/s to mm/day

sw_gs_past <- terra::rast(file.path('data', "ssrd_gs_2005_2023.tif")) / 3600       # convert J m⁻² to W m⁻² 
sw_gs_fut <- terra::rast(file.path('data', "rsds_gs_2026_2050.tif"))               # no need for conversion

et_gs_past <- terra::rast(file.path('data', "et_gs_2005_2023.tif"))                   # unit: mm/hr
et_at_past_fut_change <- terra::rast(file.path('data', "et_at_past_fut_change.tif"))   # unit: mm/hr

le_gs_past <- et_gs_past / 3600 * 2.45 * 1000000   # convert et to le by mm/hr -> mm/s -> W m⁻²
le_at_past_fut_change <- et_at_past_fut_change / 3600 * 2.45 * 1000000     # convert et to le by mm/hr -> mm/s -> W m⁻²

# calculate environmental trends
# more accurate methods, but take too much time
# ta_gs_at_past <- terra::app(ta_gs_past, function(x) {
#   time = 2005:2023
#   ifelse(any(is.na(x)), NA, summary(mblm(x ~ time))$coefficients[2, 1])
#   })
# mean(terra::values(ta_gs_at_past), na.rm=T)

ta_gs_at_past  <- terra::app(ta_gs_past, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(ta_gs_at_past), na.rm=T)
ta_gs_at_fut   <- terra::app(ta_gs_fut, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(ta_gs_at_fut), na.rm=T)
ta_ngs_at_past <- terra::app(ta_ngs_past, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(ta_ngs_at_past), na.rm=T)
ta_ngs_at_fut  <- terra::app(ta_ngs_fut, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(ta_ngs_at_fut), na.rm=T)
#
tp_gs_at_past <- terra::app(tp_gs_past, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(tp_gs_at_past), na.rm=T)
tp_ngs_at_past <- terra::app(tp_ngs_past, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(tp_ngs_at_past), na.rm=T)
tp_gs_at_fut   <- terra::app(tp_gs_fut, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(tp_gs_at_fut), na.rm=T)
tp_ngs_at_fut  <- terra::app(tp_ngs_fut, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(tp_ngs_at_fut), na.rm=T)
#
sw_gs_at_past <- terra::app(sw_gs_past, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(sw_gs_at_past), na.rm=T)
sw_gs_at_fut <- terra::app(sw_gs_fut, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(sw_gs_at_fut), na.rm=T)

le_gs_at_past <- terra::app(le_gs_past, function(x) { ifelse(any(is.na(x)), NA, sens.slope(x)$estimate) })
mean(terra::values(le_gs_at_past), na.rm=T)
le_gs_at_fut <- le_gs_at_past + le_at_past_fut_change    # LE trend is smaller. 
mean(terra::values(le_gs_at_fut), na.rm=T)

# get disturbance trends
DIST_at_wet_past <- terra::rast(file.path('data', "DIST_at_wet_past.tif"))
# assume future disturbance effects is at the same level with current disturbance effect. we can do sensitive analysis here, Junna. 
DIST_at_wet_fut <- DIST_at_wet_past

# ensure all the maps have the sampe extent and resolution
ta_gs_at_past <- terra::resample(ta_gs_at_past, land_cover)
ta_gs_at_fut <- terra::resample(ta_gs_at_fut, land_cover)
ta_ngs_at_past <- terra::resample(ta_ngs_at_past, land_cover)
ta_ngs_at_fut <- terra::resample(ta_ngs_at_fut, land_cover)

tp_gs_at_past <- terra::resample(tp_gs_at_past, land_cover)
tp_gs_at_fut <- terra::resample(tp_gs_at_fut, land_cover)
tp_ngs_at_past <- terra::resample(tp_ngs_at_past, land_cover)
tp_ngs_at_fut <- terra::resample(tp_ngs_at_fut, land_cover)

sw_gs_at_past <- terra::resample(sw_gs_at_past, land_cover)
sw_gs_at_fut <- terra::resample(sw_gs_at_fut, land_cover)

le_gs_at_past <- terra::resample(le_gs_at_past, land_cover)
le_gs_at_fut <- terra::resample(le_gs_at_fut, land_cover)

DIST_at_wet_past <- terra::resample(DIST_at_wet_past, land_cover)
DIST_at_wet_fut <- terra::resample(DIST_at_wet_fut, land_cover)

####################################projection starts here################################
# project growing-season length trends
trend_driver <- read.csv(file.path("data", "potential_trend_driver.csv"))
sites_long <- read.csv(file.path('data', 'long_term_ec_sites.csv'))

# adding IGBP, climate, and category into trend driver
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "AI")], by = "site_ID")
}

#
climate_vegs <- c('cold_ENF', 'cold_DBF', 'cold_nonforest', 'dry', 'wet')
key_vars <- c('CO2 + LE_gs_at + TA_gs_at + TA_ngs_at', 
              'SW_gs_at + Lgs_at + TA_ngs_at + CO2', 
              'P_gs_at',                          
              'LE_gs_at + TA_gs_at + SW_gs_at + P_ngs_at + SWC_gs_at', 
              'CO2 + TA_at + SW_at + DIST_at')
#
NEE_trend_past <- zone_group
NEE_trend_fut <- zone_group
terra::values(NEE_trend_past) <- NA
terra::values(NEE_trend_fut) <- NA
unique_climate_veg_groups <- levels(zone_group)[[1]]
names(unique_climate_veg_groups) <- c('value', 'label')

# parameters of random forests
m_try <- 4
min_node_size <- 5
# 
for (iclimate_veg in 1:length(climate_vegs)) {
  climate_veg <- climate_vegs[iclimate_veg]
  key_var <- key_vars[iclimate_veg]
  #
  print(climate_veg)
  #
  if (climate_veg == 'cold_ENF') {
    data_group <- trend_driver %>% 
      filter(IGBP %in% c('ENF', 'EBF') & category == "cold") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'cold_DBF') {
    data_group <- trend_driver %>% 
      filter(IGBP %in% c('DBF', 'DNF', 'MF') & category == "cold") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'cold_nonforest') {
    data_group <- trend_driver %>% 
      filter(!IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'dry') {
    data_group <- trend_driver %>% 
      filter(category == "dry") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'wet') {
    data_group <- trend_driver %>% 
      filter(category == "wet") %>% dplyr::filter(!is.na(LE_gs_at)) 
  }
  
  #--------------------project for Lgs_at for the decidulous forests------------------
  cells <- which(terra::values(zone_group) == unique_climate_veg_groups$value[unique_climate_veg_groups$label==climate_veg])
  df_project_past <- data.frame(cell_id = cells, CO2 = unname(co2_past[cells]), 
                                TA_gs_at = unname(ta_gs_at_past[cells]), TA_ngs_at = unname(ta_ngs_at_past[cells]), 
                                P_gs_at = unname(tp_gs_at_past[cells]), P_ngs_at = unname(tp_ngs_at_past[cells]), 
                                LE_gs_at = unname(le_gs_at_past[cells]), SW_gs_at = unname(sw_gs_at_past[cells]), 
                                DIST_at = unname(DIST_at_wet_past[cells]), 
                                MAP = unname(map_baseline[cells]), MAT = unname(mat_baseline[cells]), 
                                NEE = unname(nee_baseline[cells]), GPP = unname(gpp_baseline[cells]), 
                                LE = unname(le_baseline[cells]), AI = unname(aridity[cells]), 
                                source = "FLUXNET2020", 
                                climate_id = unname(koppen_climate[cells]), IGBP_id = unname(land_cover[cells]))
  df_project_past <- df_project_past %>% left_join(df_koppen_climate, by = 'climate_id') %>% 
    left_join(df_IGBP, by = 'IGBP_id')
  
  # Deal with climate classes never seen before. 
  new_climate <- setdiff(unique(df_project_past$Climate.class), unique(data_group$Climate.class))
  new_IGBP    <- setdiff(unique(df_project_past$IGBP), unique(data_group$IGBP))
  
  # Give the most common category
  if (length(new_climate) > 0) {
    df_project_past$Climate.class[df_project_past$Climate.class %in% new_climate] <- names(which.max(table(df_project_past$Climate.class)))
  }
  if (length(new_IGBP) > 0) {
    df_project_past$IGBP[df_project_past$IGBP %in% new_IGBP] <- names(which.max(table(df_project_past$IGBP)))
  }
  #
  df_project_fut <- df_project_past %>% mutate(CO2 = unname(co2_fut[cells]), 
                                               TA_gs_at = unname(ta_gs_at_fut[cells]), TA_ngs_at = unname(ta_ngs_at_fut[cells]), 
                                               P_gs_at = unname(tp_gs_at_fut[cells]), P_ngs_at = unname(tp_ngs_at_fut[cells]), 
                                               LE_gs_at = unname(le_gs_at_fut[cells]), SW_gs_at = unname(sw_gs_at_fut[cells]), 
                                               DIST_at = unname(DIST_at_wet_fut[cells]))
  if (climate_veg != 'wet') {
    df_project_past <- df_project_past %>% filter(if_all(-c(DIST_at), ~ !is.na(.)))
    df_project_fut <- df_project_fut %>% filter(if_all(-c(DIST_at), ~ !is.na(.)))      
  } else {
    df_project_past <- df_project_past %>% filter(if_all(-c(TA_ngs_at, P_ngs_at), ~ !is.na(.)))
    df_project_fut <- df_project_fut %>% filter(if_all(-c(TA_ngs_at, P_ngs_at), ~ !is.na(.))) 
  }
  #  
  if (climate_veg == 'cold_DBF') {
    # project for growing season length for this climate vegetation class
    rf_model <- ranger(
      formula = paste0("Lgs_at ~ TA_ngs_at + SW_gs_at + P_gs_at", " + MAP + MAT + NEE + GPP + LE + IGBP + Climate.class + source + AI"), 
      data = data_group, 
      mtry = m_try, 
      min.node.size = min_node_size
    )
    df_project_past$Lgs_at <- predict(rf_model, data = df_project_past)$predictions
    df_project_fut$Lgs_at <- predict(rf_model, data = df_project_fut)$predictions
  }
  #
  if (climate_veg == 'dry') {
    rf_model <- ranger(
      formula = paste0("SWC_gs_at ~ TA_gs_at + SW_gs_at + P_gs_at + LE_gs_at", " + MAP + MAT + NEE + GPP + LE + Climate.class + AI"), 
      data = data_group, 
      mtry = m_try, 
      min.node.size = min_node_size,
      importance = "permutation"
    )
    df_project_past$SWC_gs_at <- predict(rf_model, data = df_project_past)$predictions
    df_project_fut$SWC_gs_at <- predict(rf_model, data = df_project_fut)$predictions
  }
  
  #--------------project for NEE_at
  rf_model <- ranger(
    formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + NEE + GPP + LE + IGBP + Climate.class + source + AI"), 
    data = data_group,
    mtry = m_try, 
    min.node.size = min_node_size
  )
  
  if (climate_veg == 'wet') {
    # for wet-warm climate, we need to rename variable names; growing season is the whole year; 
    df_project_past <- df_project_past %>% rename(TA_at = TA_gs_at, SW_at = SW_gs_at) 
    df_project_fut  <- df_project_fut  %>% rename(TA_at = TA_gs_at, SW_at = SW_gs_at)
  }
  #
  df_project_past$NEE_at <- predict(rf_model, data = df_project_past)$predictions
  df_project_fut$NEE_at <- predict(rf_model, data = df_project_fut)$predictions
  
  # assign projections to the map
  terra::values(NEE_trend_past)[df_project_past$cell_id] <- df_project_past$NEE_at
  terra::values(NEE_trend_fut)[df_project_past$cell_id]  <- df_project_fut$NEE_at
  
}

terra::plot(NEE_trend_past)
terra::plot(NEE_trend_fut)
terra::plot(NEE_trend_fut - NEE_trend_past) 

hist(terra::values(NEE_trend_fut - NEE_trend_past))
density(terra::values(NEE_trend_fut - NEE_trend_past), na.rm=T)
quantile(terra::values(NEE_trend_fut - NEE_trend_past), probs = c(0.95, 0.75, 0.5, 0.25, 0.05), na.rm = T)

# writeRaster(NEE_trend_past, filename = file.path('data', 'NEE_trend_2005-2023.tif'), overwrite=TRUE)
# writeRaster(NEE_trend_fut, filename = file.path('data', 'NEE_trend_2026-2045.tif'), overwrite=TRUE)
writeRaster(NEE_trend_fut - NEE_trend_past, filename = file.path('data', 'NEE_trend_change.tif'), overwrite=TRUE)

