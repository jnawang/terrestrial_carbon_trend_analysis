# A function to compare two SWC products (measured vs modelled)
# 1/28/2026, Junna Wang

rm(list=ls())

# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'
#
library(librarian)
shelf(tidyverse, cowplot, ggpubr, trend)

# function to calculate sens' slope
cal_sens_slope <- function(x) {
  x <- x[!is.na(x)]
  # x <- scale(x)
  return(as.numeric(sens.slope(x)[c(1,3)]))
}

cal_sens_slope_scale <- function(x) {
  x <- x[!is.na(x)]
  x <- scale(x)
  return(as.numeric(sens.slope(x)[c(1,3)]))
}

# data I needed: site, yearStart, yearEnd, trend_obsv, trend_anls
SWC_trend <- data.frame(site_ID = character(), yearStart = integer(), yearEnd = integer(), 
                        trend_obsv = double(), trend_obsvp = double(), 
                        trend_obsv_scale = double(), trend_obsv_scalep = double(), 
                        trend_obsv_Y = double(), trend_obsvp_Y = double(), 
                        trend_obsv_Y_scale = double(), trend_obsv_Y_scalep = double(), 
                        trend_mold = double(), trend_moldp = double(), 
                        trend_mold_scale = double(), trend_mold_scalep = double())


# Find all the sites with FLUXNET products with good measured SWC values
sites_long    <- read.csv(file.path('data', 'long_term_ec_sites.csv'))
# How many sub time series do we have?
rm_years   <- read.csv(file.path('data', 'years_to_remove.csv'))

source(file.path("scripts", "00_functions.R"))
df_file_structure <- get_rawflux_file_structure(dir_rawdata)

# read analyzed SWC
SWC_ERA5 <- read.csv(file.path('/Volumes/MaloneLab/Research/Stability_Project', 'Thermal_Acclimation', "ERA5_daily_swc_1990_2024_allsites.csv"))
SWC_ERA5$date <- ymd(SWC_ERA5$date)
SWC_ERA5$year <- year(SWC_ERA5$date) 
sites_SWC_ERA5 <- unique(SWC_ERA5$name)

for (i in 1:nrow(sites_long)) {
  if (! sites_long$source[i] %in% c('FLUXNET2020', 'AmeriFlux_FLUXNET', 'FLUXNET2020_ICOS')) next
  # i = 15
  site_ID <- sites_long$site_ID[i]
  if (!site_ID %in% sites_SWC_ERA5) next
  
  print(i)
  print(site_ID)
  folder <- unlist(df_file_structure$folders[df_file_structure$sources == sites_long$source[i]])
  DDfile <- df_file_structure$DDfiles[df_file_structure$sources == sites_long$source[i]]
  YYfile <- df_file_structure$YYfiles[df_file_structure$sources == sites_long$source[i]]
  a_DD <- read.csv(folder[grepl(pattern=paste0(site_ID, DDfile), folder)])
  a_YY <- read.csv(folder[grepl(pattern=paste0(site_ID, YYfile), folder)])
  a_DD[a_DD==-9999] <- NA
  a_YY[a_YY==-9999] <- NA
  #
  a_DD <- a_DD %>% mutate(TIMESTAMP = ymd(a_DD$TIMESTAMP), YEAR = year(TIMESTAMP))
  a_DD_good <- a_DD %>% filter(!YEAR %in% rm_years$years[rm_years$site_ID==site_ID])
  a_YY_good <- a_YY %>% filter(!TIMESTAMP %in% rm_years$years[rm_years$site_ID==site_ID])
  #
  # 
  if ("SWC_F_MDS_1" %in% names(a_DD_good) & sum(!is.na(a_YY_good$SWC_F_MDS_1)) > 5) {
    print("SWC_F_MDS_1")
    plot(a_DD_good$TIMESTAMP, a_DD_good$SWC_F_MDS_1, main = site_ID)
  } else if ("SWC_F_MDS" %in% names(a_DD_good) & sum(!is.na(a_YY_good$SWC_F_MDS)) > 5) {
    print("SWC_F_MDS")
  } else {
    next
  }
  #
  SWC_trend[i, 1] <- site_ID
  SWC_trend[i, 2] <- min(a_DD_good$YEAR, na.rm = T)
  SWC_trend[i, 3] <- max(a_DD_good$YEAR, na.rm = T)
  SWC_trend[i, 4:5] <- cal_sens_slope(a_DD_good$SWC_F_MDS_1)
  SWC_trend[i, 6:7] <- cal_sens_slope_scale(a_DD_good$SWC_F_MDS_1)
  SWC_trend[i, 8:9] <- cal_sens_slope(a_YY_good$SWC_F_MDS_1)
  SWC_trend[i, 10:11] <- cal_sens_slope_scale(a_YY_good$SWC_F_MDS_1)
  # compare these sites with climate reanalysis products
  SWC_ERA5_site <- SWC_ERA5 %>% filter(name == site_ID) %>% filter(between(year, SWC_trend[i, 2], SWC_trend[i, 3]))
  SWC_trend[i, 12:13] <- cal_sens_slope(SWC_ERA5_site$SWC)
  SWC_trend[i, 14:15] <- cal_sens_slope_scale(SWC_ERA5_site$SWC)
}


# the last thing is to ask models to obtain more results. 
SWC_trend %>% filter(!site_ID %in% c("CA-Mer", "US-BZS", "US-Ho1", "US-ICs", "UC-ICt",
                                     "UC-MOz", "US-NC4", "US-Syv",
                                     "US-UMB", "US-Wkg",  "BE-Bra", 
                                     "BE-Dor", "BE-Vie", "CZ-BK1", "CZ-Stn", "CZ-wet", 
                                     "RU-Fyo", "IL-Yat", "SE-Deg", "IT-Ren", "FR-Fon", 
                                     "FR-FBn", "DK-Sor", "DE-Tha",  "DE-Obe", "US-Ton")) %>% 
  # filter(abs(trend_obsv) < 0.004) %>%
  # ggplot(aes(x=trend_obsv, y=trend_mold)) +
  ggplot(aes(x=trend_obsv_scale, y=trend_mold_scale)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  # stat_regline_equation(
  #   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))
  # ) +
  theme_classic()




