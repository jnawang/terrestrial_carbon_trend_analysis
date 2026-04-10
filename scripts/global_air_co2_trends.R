# This script calculate the amplitude of annual seasonal cycle of different climate classes. 
# Criteria of selecting a site-year:
# 1. at least one observation in a month. 
# 2. at least 15 years between 1995 and 2025. 
# 3. the annual amplitude is > 7 ppm. 

# Junna Wang, 4/4/2026

# New ideas for this analysis
# cluster analysis of these sites
# use relative magnitude
# remove the year with measurement < 2 days a year
# could we just look the amplitude of these valleys? the contribution of peaks and valleys on the amplitude. 
# should we look at the zero-crossing date?
# if we only use the sites above 50 degree N, could we see the increase in respiration activity ? or the sites between 50 and 70 degree. 
# I really think the shape of the curve is likely to change; larger peak and lower valley; 





library(librarian)
shelf(ggplot2, dplyr, tidyr, cowplot, ggpubr, trend, amerifluxr, lme4, ncdf4, lubridate, tibble, stringr, zoo, terra, rnaturalearthdata, tidyterra, lmerTest)


# look at atmospheric CO2
read_obspack_metadata <- function(file) {
  # file <- files[1]
  print(file)
  
  # read file as text
  txt <- readLines(file, warn = FALSE)
  
  # keep key : value lines
  meta_lines <- txt[grepl(":", txt)]
  
  # split into key / value
  kv <- strsplit(meta_lines, "\\s*:\\s*")
  meta_df <- data.frame(
    key = sapply(kv, `[`, 1),
    value = sapply(kv, `[`, 2),
    stringsAsFactors = FALSE
  )
  
  # extract fields
  dataset_start <- ymd_h(
    meta_df$value[meta_df$key == "dataset_start_date"],
    tz = "UTC"
  )
  dataset_stop <- ymd_h(
    meta_df$value[meta_df$key == "dataset_stop_date"],
    tz = "UTC"
  )
  
  data_duration_years <- time_length(
    interval(dataset_start, dataset_stop),
    "year"
  )
  
  data.frame(
    site_name = meta_df$value[meta_df$key == "site_name"],
    # site_country = meta_df$value[meta_df$key == "site_country"],
    site_latitude = as.numeric(meta_df$value[meta_df$key == "site_latitude"]),
    site_longitude = as.numeric(meta_df$value[meta_df$key == "site_longitude"]),
    site_elevation = as.numeric(meta_df$value[meta_df$key == "site_elevation"]),
    dataset_name = meta_df$value[meta_df$key == "dataset_name"],
    dataset_project = meta_df$value[meta_df$key == "dataset_project"],
    dataset_start_date = dataset_start,
    dataset_stop_date = dataset_stop,
    data_duration_years = round(data_duration_years, 2),
    stringsAsFactors = FALSE
  )
}

# files <- list.files(file.path(dir_rawdata, "obspack_co2", "metadata"), pattern = '.txt$', full.names = T)
files <- list.files('/Users/junnawang/YaleLab/data_CO2/obspack_co2_1_GLOBALVIEWplus_v11.0_2026-01-12_daily/metadata', pattern = '.txt', full.names = T)
all_sites_info <- do.call(rbind, lapply(files, read_obspack_metadata))

# ~ 300 stations have >10 years of data
# ~ 128 stations have >20 years of data
all_sites_info_10yr <- all_sites_info %>% filter(data_duration_years >= 20) %>% # filter(site_latitude > - 100) %>% 
  filter(dataset_project %in% c("surface-insitu", "shipboard-insitu", "surface-flask"))  # "tower-insitu", 
all_sites_info_10yr$magl <- as.numeric(str_extract(all_sites_info_10yr$dataset_name, "(?<=-)[0-9]+(?=magl)"))

# remove these duplicated datasets for one location, only save one dataset with longest data duration. 
all_sites_info_10yr_unique <- all_sites_info_10yr %>% 
  group_by(site_name, site_latitude, site_longitude) %>%
  slice_max(data_duration_years, n = 1) %>%  # keeps row with layer_height closest to 30
  ungroup() 

all_sites_info_10yr_unique$site_ID <- sapply(all_sites_info_10yr_unique$dataset_name, function(strname) { unlist(strsplit(strname, split = '_'))[2] })

# the spatial distribution of these sites
mundi <- vect(rnaturalearthdata::coastline110)
mundi <- crop(mundi,ext(c(-180,180,-60,90)))
#
p1 <- ggplot() +
  geom_spatvector(data = mundi, linewidth = .2) +
  geom_hline(yintercept=0, linetype = "dotted", colour = "grey") + 
  geom_point(data=all_sites_info_10yr_unique, aes(x=site_longitude, y=site_latitude), shape = 21, size = 2, stroke = 0.5, color='black') +
  scale_fill_manual(values=c('wet' = "#4393C3", 'dry' = "#E4003A", 'cold' = "#B790D4"), 
                    labels=c("wet"='Wet: temperate and tropical climates', 
                             "dry"="Dry: semi-arid climate",
                             "cold"="Cold: semi-arid climate")) + 
  labs(x='Longitude', y='Latitude', fill='Climate category') +
  coord_sf(expand=FALSE, xlim = c(-180, 180), ylim = c(-90, 90)) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black")) +
  theme(legend.position = c(0.25, 0.2),
        legend.text = element_text(size = 12), 
        legend.background = element_rect(fill = alpha('white', 0.3)))
p1

# Classification of these sites: dry, cold, warm
climate_koppen   <- terra::rast(file.path('data', 'koppen_geiger_0p00833333.tif'))
# plot(climate_koppen)
# aridicity
AI_0.25 <- terra::rast(file.path('data', 'Global-AI.tif'))
AI_0.25 <- AI_0.25 / 10000       # values are enlarged by 10000
plot(AI_0.25)
#
# get climate for the European and OZ fluxes based on coordinates
xy <- data.frame(x=all_sites_info_10yr_unique$site_longitude, 
                 y=all_sites_info_10yr_unique$site_latitude)
#
all_sites_info_10yr_unique$climate_code <- terra::extract(climate_koppen, xy)[,2]
all_sites_info_10yr_unique$ai <- terra::extract(AI_0.25, xy)[,2]
all_sites_info_10yr_unique <- all_sites_info_10yr_unique %>% mutate(koppen = case_when(climate_code %in% c(12:13, 15:30) ~ 'cold', 
                                                         climate_code %in% c(1:3, 11, 14) & ai > 0.5 ~ 'warm', 
                                                         climate_code %in% c(1:3, 11, 14) & ai < 0.5 ~ 'dry', 
                                                         climate_code %in% c(4:10) ~ 'dry', 
                                                         TRUE ~ 'marine'))

baseline_sites_info <- all_sites_info_10yr_unique 

# baseline_sites_info <- all_sites_info_10yr_unique %>% filter(dataset_name %in% c(
#   # "co2_brw_surface-insitu_1_allvalid", "co2_mlo_surface-insitu_1_allvalid",
#   # "co2_spo_surface-insitu_1_allvalid", "co2_syo_surface-insitu_8_allvalid",
#   "co2_alt_surface-insitu_6_allvalid", "co2_nwr_surface-flask_1_representative",
#   "co2_cgo_surface-flask_1_representative", "co2_uta_surface-flask_1_representative",
#   "co2_lef_surface-flask_1_representative", "co2_sgp_surface-insitu_64_allvalid-60magl",
#   "co2_hei_surface-insitu_22_allvalid"
#   # , 
#   # "co2_hei_surface-insitu_22_allvalid", 
#   # "co2_mnm_surface-insitu_19_representative", 	
#   # "co2_yon_surface-insitu_19_representative", "co2_pal_surface-insitu_30_allvalid", 
#   # "co2_cpt_surface-insitu_36_marine", 
#   # "co2_kas_surface-insitu_53_allvalid"
# ))

# take a look at how CO2 concentration change over elevation; the magnitude of CO2 variation should decline with elevation
df_obspack_co2 <- data.frame()
dir <- "/Users/junnawang/YaleLab/data_CO2/obspack_co2_1_GLOBALVIEWplus_v11.0_2026-01-12_nc/data/nc/"
# dir <- file.path(dir_rawdata, "obspack_co2", "obspack_co2_1_GLOBALVIEWplus_v11.0_2026-01-12", "data", "nc")


for (dataset in baseline_sites_info$dataset_name) {
  file <- paste0(dataset, '.nc')
  nc_file <- nc_open(file.path(dir, file))
  df <- data.frame(dataset_name = dataset,
                   datetime = ymd_hms(ncvar_get(nc_file, "datetime")), 
                   co2 = ncvar_get(nc_file, "value") * 1000000)
  
  plot(df$datetime, df$co2, main = dataset)
  
  df_obspack_co2 <- bind_rows(df_obspack_co2, df)
}

df_obspack_co2_daily <- df_obspack_co2 %>% mutate(date = as.Date(datetime)) %>% 
  group_by(dataset_name, date) %>% summarise(co2 = mean(co2, na.rm=T), .groups = "drop") %>%
  mutate(year = year(date), month = month(date)) %>%
  left_join(all_sites_info_10yr_unique[, c('dataset_name', 'site_ID')], by = 'dataset_name')

# remove site years without enough measurements; at least one observation per month.
site_years_all <- df_obspack_co2_daily %>% filter(!is.na(co2)) %>%
  distinct(site_ID, year)  %>% arrange(site_ID, year)

site_years_bad <- df_obspack_co2_daily %>% filter(!is.na(co2)) %>%
  group_by(site_ID, year) %>% summarise(nmonth = n_distinct(month), .groups = "drop") %>%
  filter(nmonth < 12) %>% select(site_ID, year)

sites_good <- site_years_all %>%
  anti_join(site_years_bad, by = c("site_ID", "year")) %>% 
  filter(between(year, 1995, 2025)) %>%
  group_by(site_ID) %>% summarise(nyear = n_distinct(year), .groups = "drop") %>%
  filter(nyear >= 15) %>% 
  pull(site_ID) %>% as.character()

sites_to_remove <- setdiff(unique(site_years_all$site_ID), sites_good)
  
df_obspack_co2_daily <- df_obspack_co2_daily %>% filter(!site_ID %in% sites_to_remove)
site_years_bad <- site_years_bad %>% filter(!site_ID %in% sites_to_remove)

# get site-specific characteristics
for (dataset in unique(df_obspack_co2_daily$site_ID)) {
  # 3-day, 5-day, 7-day average
  df_brw <- df_obspack_co2_daily %>% filter(site_ID == dataset) 
  df_brw$co2_m3 <- rollmean(df_brw$co2, k = 3, align = 'right', fill = NA) 
  df_brw$co2_m5 <- rollmean(df_brw$co2, k = 5, align = 'right', fill = NA) 
  df_brw$co2_m7 <- rollmean(df_brw$co2, k = 7, align = 'right', fill = NA) 
  df_brw$co2_m10 <- rollmean(df_brw$co2, k = 10, align = 'right', fill = NA) 
  
  df_brw$decimal_date = decimal_date(df_brw$date)
  
#  df_brw <- df_brw %>% filter(decimal_date <= 2025.002739726)
  # if the last row is close to the end of a year, adding one row to indicate the end of a year.
  nrows <- nrow(df_brw)
  if (df_brw$decimal_date[nrows] - df_brw$year[nrows] > 0.95) {
    df_brw <- bind_rows(df_brw, data.frame(decimal_date = df_brw$year[nrows] + 1, 
                                           co2 = df_brw$co2[nrows] + 0.1))
  }
  #
  p2 <- df_brw %>% 
    filter(year >= 1990) %>% 
    ggplot() +
    geom_line(aes(x=date, y=co2), col = 'black') +
    geom_line(aes(x=date, y=co2_m10), col = 'purple') +
    labs(title = dataset) +
    theme_classic2()
  print(p2)
  
  # #
  # df_brw_year <- df_brw %>% group_by(year) %>% summarise(co2_min = min(co2_m10), co2_max = max(co2_m10), 
  #                                                        co2_mean = mean(co2_m10), mag = co2_max - co2_min, 
  #                                                        mag_min = co2_mean - co2_min, mag_max = co2_max - co2_mean)
  # 
  # p3 <- ggplot(df_brw_year) +
  #   geom_line(aes(x = year, y = co2_min), col = 'blue') + 
  #   geom_line(aes(x = year, y = co2_max), col = 'red') + 
  #   geom_line(aes(x = year, y = co2_mean), col = 'black') + 
  #   labs(title = dataset) + 
  #   theme_classic2()
  # print(p3)
  # 
  # p4 <- ggplot(df_brw_year) +
  #   geom_line(aes(x = year, y = mag_min), col = 'blue') +
  #   geom_line(aes(x = year, y = mag), col = 'black') +
  #   geom_line(aes(x = year, y = mag_max), col = 'red') +
  #   labs(title = dataset) +
  #   theme_classic2()
  # print(p4)
  # 
  # summary(lm(data = df_brw_year, mag ~ year))
  
  write.table(df_brw[, c("decimal_date", "co2")], file = paste0('ccg_filter/data/', dataset, '.txt'), sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#-------------------------------------------------------------------------------------------
# Call the python program to calculate amplitude. Here is the command line to call the py scripts. 
# for f in data/*.txt; do
# base=$(basename "$f")
# python3 ccgcrv.py --file "output/$base" --amp "$f"
# done

#------------------
# I did not see what I expected to see. peak trends are increasing not faster than valley trends. 
# except the Harwii and germany sites, all other sites co2 magnitude is increasing. 
# I need to detrend and save the result. 
# autumn peak or summer peak. duration of net carbon uptake. 

# What I need to show is that correlation of AMP and max peak is strengthening; contribution of maximum peak is increasing.
# The AMP increasing trends are declining. I also need to calculate trends. the minimum AMP trends are declining, and the maximum AMP trends are increasing. 

# which stations are close to global drylands. 

# I need to know amp and two timings of zero crossing
# Results I need: year, amp, peak, valley, peak date, valley date, zero_up, zero_dn

# https://gml.noaa.gov/ccgg/mbl/crvfit/index.html

# co2_amp.csv is the result from the NOAA package. 
co2_amp <- data.frame()
outfiles <- list.files(file.path("ccg_filter", "output"), full.names = F)
for (outfile in outfiles) {
  df <- read.table(file.path("ccg_filter", "output", outfile), sep = '', skip = 3, header = FALSE)
  names(df) <- c('year', 'amp', 'amp_peak', 'month_max', 'day_max', 
                                'amp_valley', 'month_valley', 'day_valley')
  # remove bad years
  site <- substring(outfile, 1, 3)
  if (sum(site_years_bad$site_ID == site) > 0) {
    df <- df %>% filter(!year %in% site_years_bad$year[site_years_bad$site_ID == site])
  }
  co2_amp <- bind_rows(co2_amp, data.frame(site_ID = site, df))
}

co2_amp_feature <- co2_amp %>% group_by(site_ID) %>% summarise(amp_mean = mean(amp, na.rm = T)) %>%
  left_join(all_sites_info_10yr_unique[, c('site_ID', 'climate_code', 'koppen', 'site_name', "site_latitude", "site_longitude")], by = 'site_ID')

# mannually assign climate_code to sites not identified in the global climate map.
co2_amp_feature$koppen[co2_amp_feature$site_ID %in% c("coi", "cya", "hba", "ice", "mqa", "psa", "syo")] <- 'cold'
co2_amp_feature$koppen[co2_amp_feature$site_ID %in% c("hat", "bmw", "key", "kum", "mnm", "tap", "mid", "mlo")] <- 'warm'
co2_amp_feature$koppen[co2_amp_feature$site_ID %in% c("nwr")] <- 'dry'

if (!'koppen' %in% names(co2_amp)) {
  co2_amp <- co2_amp %>% left_join(co2_amp_feature[, c('site_ID', 'koppen', "site_latitude", "site_longitude", "amp_mean")], by = 'site_ID')
}

###########
sites <- unique(co2_amp$site_ID)
window_size = 15

source(file.path("scripts", "00_functions.R"))

# Characteristics of annual seasonal cycle for each site
df_amp_at <- data.frame(site_ID = character(), year_mid = integer(), 
                        amp_at = double(), amp_at_peak = double(), amp_at_valley = double())
id = 1
for (site in sites) {
  co2_amp_site <- co2_amp[co2_amp$site_ID == site, ]
  year_start <- min(co2_amp_site$year)
  year_end <- max(co2_amp_site$year)
  nwindow = max(year_end - year_start + 1 - window_size + 1, 1)
  for (iwindow in 1:nwindow) {
    years_window <- (year_start + iwindow - 1):(year_start + iwindow + window_size - 2)
    df <- co2_amp_site %>% filter(year %in% years_window) %>% arrange(year)
    df_amp_at[id, "site_ID"] <- site
    df_amp_at[id, "year_mid"] <- mean(years_window)
    df_amp_at[id, "amp_at"] <- cal_sens_slope(df$year, df$amp)                # sens.slope(df$amp) 
    df_amp_at[id, "amp_at_peak"] <- cal_sens_slope(df$year, df$amp_peak)      # 
    df_amp_at[id, "amp_at_valley"] <- cal_sens_slope(df$year, df$amp_valley)  # 
    #
    df_amp_at[id, 6] <- cal_sens_slope_scale(df$year, df$amp, df$amp)[2]
    df_amp_at[id, 7] <- cal_sens_slope_scale(df$year, df$amp_peak, df$amp_peak)[2]
    df_amp_at[id, 8] <- cal_sens_slope_scale(df$year, df$amp_valley, df$amp_valley)[2]
    id = id + 1
  }
}

site = "mlo" # 'brw'
co2_amp %>% filter(site_ID == site) %>% 
  ggplot(aes(x = year)) +  # , col = site
  geom_line(aes(y = amp), linetype = 1) +
  geom_line(aes(y = amp_peak), linetype = 2) + 
  geom_line(aes(y = amp_valley), linetype = 3) + 
  geom_hline(yintercept = 0, linetype = 3) +
  labs(title = site) + 
  theme_classic2()

df_amp_at %>% filter(site_ID == site) %>%      # "mlo", 'brw', 'alt', 'syo'
  ggplot(aes(x = year_mid)) +  # , col = site
  geom_line(aes(y = amp_at), linetype = 1) +
  geom_line(aes(y = amp_at_peak), linetype = 2) +      # dashed.
  geom_line(aes(y = amp_at_valley), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(title = site) + 
  theme_classic2()

###########overall trend of all climates
df_amp_at_group <- data.frame()
# using all sites or sites within certain category. 
climate_vegs <- c('cold', 'dry', 'warm')

# remove the sites with too small amp. What should be the criteria?   # amp >= either 1 or 4 ppm
for (climate_veg in climate_vegs) {
  NEEtrend <- data.frame(year_mid=integer(), nsite = integer(),
                         amp=double(), Tlmer=double(), Tlmerlow=double(), Tlmerup=double(),
                         amp_peak=double(), Tlmer_peak=double(), Tlmerlow_peak=double(), Tlmerup_peak=double(),
                         amp_valley=double(), Tlmer_valley=double(), Tlmerlow_valley=double(), Tlmerup_valley=double()
                         )  
  #
  if (climate_veg == 'cold') {
    data <- co2_amp %>% filter(koppen == "cold") %>% filter(amp_mean > 5) # filter(between(site_latitude, 45, 90))
  } else if (climate_veg == 'warm') {
    data <- co2_amp %>% filter(koppen == "warm") %>% filter(amp_mean > 5)
  } else if (climate_veg == 'dry') {
    data <- co2_amp %>% filter(koppen == "dry") %>% filter(amp_mean > 5)
  }
  # data <- co2_amp 
  
  #
  year_start <- max(min(data$year), 1990)
  year_end <- min(max(data$year), 2024)
  #
  for (iwindow in year_start:(year_end - window_size + 1)) {
    print(iwindow)
    subdata <- data %>% filter(between(year, iwindow, iwindow + window_size - 1))   
    mod_lmer <- lmer(data=subdata, amp ~ year + (1 | site_ID))
    mod_lmer_peak <- lmer(data=subdata, amp_peak ~ year + (1 | site_ID))
    mod_lmer_valley <- lmer(data=subdata, amp_valley ~ year + (1 | site_ID))
    #
    NEEtrend[iwindow-year_start, 1] <- iwindow + (window_size - 1) / 2
    NEEtrend[iwindow-year_start, 2] <- length(unique(subdata$site_ID))
    NEEtrend[iwindow-year_start, 3] <- mean(subdata$amp[subdata$year == NEEtrend[iwindow-year_start, 1]], na.rm = T)
    NEEtrend[iwindow-year_start, 4] <- summary(mod_lmer)$coefficients[2,c(1)]
    NEEtrend[iwindow-year_start, 5:6] <- confint(mod_lmer, parm = "beta_", method = "profile")[2,1:2]
    NEEtrend[iwindow-year_start, 7] <- mean(subdata$amp_peak[subdata$year == NEEtrend[iwindow-year_start, 1]], na.rm = T)
    NEEtrend[iwindow-year_start, 8] <- summary(mod_lmer_peak)$coefficients[2,c(1)]
    NEEtrend[iwindow-year_start, 9:10] <- confint(mod_lmer_peak, parm = "beta_", method = "profile")[2,1:2]
    NEEtrend[iwindow-year_start, 11] <- mean(subdata$amp_valley[subdata$year == NEEtrend[iwindow-year_start, 1]], na.rm = T)
    NEEtrend[iwindow-year_start, 12] <- summary(mod_lmer_valley)$coefficients[2,c(1)]
    NEEtrend[iwindow-year_start, 13:14] <- confint(mod_lmer_valley, parm = "beta_", method = "profile")[2,1:2]
  }
  
  #-----------
  p0 <- ggplot(data = NEEtrend, aes(x=year_mid)) +
    geom_line(aes(y = amp), col = 'black') +
    geom_line(aes(y = amp_peak), col = 'red') +
    geom_line(aes(y = amp_valley), col = 'blue') +
    geom_hline(yintercept = 0) +
    lims(x=c(1990, 2020)) +
    labs(title = climate_veg) +
    theme_classic()
  print(p0)
  
  p1 <- ggplot(data = NEEtrend, aes(x=year_mid, y = Tlmer)) +
    geom_line(col = 'green') +
    geom_point(col = 'green') +
    geom_errorbar(aes(ymin=Tlmerlow, ymax=Tlmerup)) +
    geom_hline(yintercept = 0) +
    lims(x=c(1990, 2020)) +
    labs(title = climate_veg) +
    theme_classic()
  print(p1)

  p2 <- ggplot(data = NEEtrend, aes(x=year_mid, y = Tlmer_peak)) +
    geom_line(col = 'green') +
    geom_point(col = 'green') +
    geom_errorbar(aes(ymin=Tlmerlow_peak, ymax=Tlmerup_peak)) +
    geom_hline(yintercept = 0) +
    lims(x=c(1990, 2020)) +
    labs(title = climate_veg) +
    theme_classic()
  print(p2)  
  
  p3 <- ggplot(data = NEEtrend, aes(x=year_mid, y = Tlmer_valley)) +
    geom_line(col = 'green') +
    geom_point(col = 'green') +
    geom_errorbar(aes(ymin=Tlmerlow_valley, ymax=Tlmerup_valley)) +
    geom_hline(yintercept = 0) +
    lims(x=c(1990, 2020)) +
    labs(title = climate_veg) +
    theme_classic()
  print(p3)    
  
  # merge them together
  df_amp_at_group <- bind_rows(df_amp_at_group, data.frame(climate_veg = climate_veg, NEEtrend))
}

write.csv(df_amp_at_group, file.path("data", "air_co2_amp_trend.csv"), row.names = T)
# Lingering questions for this analysis. 
# How to classify these sites?
# a better way to represent the intensity of net uptake from the terrestrial ecosystems. 
# the threshold of classification
# warm: 8 sites
# dry: 13 sites
# cold: 29 sites

