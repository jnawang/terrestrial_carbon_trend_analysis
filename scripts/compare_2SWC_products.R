# A function to compare two SWC products (measured vs modelled)
# 1/28/2026, Junna Wang

# important: SWC_ERA5 does not include data from IT-Noe. 
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
                        trend_mold_scale = double(), trend_mold_scalep = double(),
                        trend_mold_Y = double(), trend_moldp_Y = double(), 
                        trend_mold_Y_scale = double(), trend_mold_Y_scalep = double())


# Find all the sites with FLUXNET products with good measured SWC values
sites_long    <- read.csv(file.path('data', 'long_term_ec_sites.csv'))
# How many sub time series do we have?
rm_years   <- read.csv(file.path('data', 'years_to_remove.csv'))

source(file.path("scripts", "00_functions.R"))
df_file_structure <- get_rawflux_file_structure(dir_rawdata)

# read analyzed SWC
SWC_ERA5 <- read.csv(file.path(dir_rawdata, "ERA5L_daily_SWC_trend_sites_1990_2025.csv"))
names(SWC_ERA5)[3] <- 'SWC'

# terra soil water data
# SWC_ERA5 <- read.csv(file.path(dir_rawdata, "terra_SWC_allsites.csv"))
# names(SWC_ERA5) <- c("date", "SWC", "name")

#
SWC_ERA5$date <- ymd(SWC_ERA5$date)
SWC_ERA5$year <- year(SWC_ERA5$date) 
sites_SWC_ERA5 <- unique(SWC_ERA5$name)

sites_SWC_remove <- c("RU-Fyo", "DE-Obe", "CZ-Stn", "SE-Htm", "SE-Deg", "FR-Bil", "DK-Sor",
                      "CZ-wet", "CZ-BK1", "CH-Dav", "BE-Vie", "BE-Dor", "US-Wkg", "BE-Bra", 
                      "US-Syv", "US-Ro4", "US-Ho1", "US-Esm", "US-EML", "US-BZS", "FR-FBn",
                      "IL-Yat", "IT-Ren", "FR-Fon", "IT-SRo", "IT-Ro2", "IT-Noe", "NL-Loo")
txt <- "
site_ID, goodyear_start, goodyear_end
ES_LJu, 2007, 2025
IT-Ren, 2000, 2009
IT-Cp2, 2013, 2025
DE-Tha, 1997, 2018
US-UMd, 2016, 2025
US-UMB, 2008, 2025
US-Ton, 2003, 2020
US-NC4, 2014, 2020
US-MOz, 2012, 2025
US-Kon, 2004, 2018
UC-ICt, 2013, 2025
US-ICs, 2014, 2025
CA-Mer, 2010, 2020
AU-Stp, 2009, 2016
AU-Cow, 2012, 2025
JP-Tef, 2006, 2021
JP-Ta2, 2007, 2016
"
sites_SWC_trunck <- read.csv(text = txt)

for (i in nrow(sites_long)) {
  # if (! sites_long$source[i] %in% c('FLUXNET2020', 'AmeriFlux_FLUXNET', 'FLUXNET2020_ICOS')) next
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
  # if SWC_F_MDS in a_DD and a_YY, but SWC_F_MDS_1 not, add one more column SWC_F_MDS_1
  if ("SWC_F_MDS" %in% names(a_DD) & ! "SWC_F_MDS_1" %in% names(a_DD) ) {
    a_DD$SWC_F_MDS_1 <- a_DD$SWC_F_MDS
    a_YY$SWC_F_MDS_1 <- a_YY$SWC_F_MDS
  }

  #
  a_DD <- a_DD %>% mutate(TIMESTAMP = ymd(a_DD$TIMESTAMP), YEAR = year(TIMESTAMP))
  a_DD_good <- a_DD %>% filter(!YEAR %in% rm_years$years[rm_years$site_ID==site_ID])
  a_YY_good <- a_YY %>% filter(!TIMESTAMP %in% rm_years$years[rm_years$site_ID==site_ID])
  #
  year_start <- min(a_YY_good$TIMESTAMP)
  year_end <- max(a_YY_good$TIMESTAMP)
  
  if ("SWC_F_MDS_1" %in% names(a_DD_good) & sum(!is.na(a_YY_good$SWC_F_MDS_1)) > 5) {
    # the unit of SWC is inconsistent; make them consistent
    if (mean(a_YY_good$SWC_F_MDS_1, na.rm=T) < 1) {
      a_DD_good$SWC_F_MDS_1 <- a_DD_good$SWC_F_MDS_1 * 100
      a_YY_good$SWC_F_MDS_1 <- a_YY_good$SWC_F_MDS_1* 100   
    }
    #
    print("SWC_F_MDS_1")
    plot(a_DD_good$TIMESTAMP, a_DD_good$SWC_F_MDS_1, main = site_ID)
    #
    if (site_ID %in% sites_SWC_remove) {
      next
    } else if (site_ID %in% sites_SWC_trunck$site_ID) {
      year_start <- sites_SWC_trunck$goodyear_start[sites_SWC_trunck$site_ID == site_ID]
      year_end <- sites_SWC_trunck$goodyear_end[sites_SWC_trunck$site_ID == site_ID]
    }
  } else {
    next
  }
  #
  a_DD_good <- a_DD_good %>% filter(between(YEAR, year_start, year_end))
  a_YY_good <- a_YY_good %>% filter(between(TIMESTAMP, year_start, year_end))
  # 
  SWC_trend[i, 1] <- site_ID
  SWC_trend[i, 2] <- min(a_DD_good$YEAR, na.rm = T)
  SWC_trend[i, 3] <- max(a_DD_good$YEAR, na.rm = T)
  SWC_trend[i, 4:5] <- cal_sens_slope(a_DD_good$SWC_F_MDS_1)
  SWC_trend[i, 6:7] <- cal_sens_slope_scale(a_DD_good$SWC_F_MDS_1)
  SWC_trend[i, 8:9] <- cal_sens_slope(a_YY_good$SWC_F_MDS_1)
  SWC_trend[i, 10:11] <- cal_sens_slope_scale(a_YY_good$SWC_F_MDS_1)
  # compare these sites with climate reanalysis products
  SWC_ERA5_site <- SWC_ERA5 %>% filter(name == site_ID) %>% filter(between(year, year_start, year_end))
  SWC_trend[i, 12:13] <- cal_sens_slope(SWC_ERA5_site$SWC * 100)
  SWC_trend[i, 14:15] <- cal_sens_slope_scale(SWC_ERA5_site$SWC * 100)
  SWC_ERA5_site_Y <- SWC_ERA5_site %>% group_by(year) %>% summarise(SWC = mean(SWC)) 
  SWC_trend[i, 16:17] <- cal_sens_slope(SWC_ERA5_site_Y$SWC * 100)
  SWC_trend[i, 18:19] <- cal_sens_slope_scale(SWC_ERA5_site_Y$SWC * 100)
}

# the last thing is to ask models to obtain more results. 
SWC_trend %>% 
  filter(abs(trend_obsv) < 0.005) %>%
  # ggplot(aes(x=trend_obsv, y=trend_mold)) +
  # ggplot(aes(x=trend_obsv_Y, y=trend_mold_Y)) +
  # ggplot(aes(x=trend_obsv_Y_scale, y=trend_mold_Y_scale)) +
  ggplot(aes(x=trend_obsv_scale, y=trend_mold_scale)) +  # trend_mold
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  # stat_regline_equation(
  #   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))
  # ) +
  theme_classic() #

#  we have 95 site or this comparison
# one site IT-Noe, their modeled SWC IS missing, we will use measured data. 
SWC <- SWC_ERA5 %>% group_by(name)  %>% summarise(SWC = mean(SWC, na.rm = T))
 






