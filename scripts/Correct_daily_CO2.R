# This script is used to correct sites with abnormal CO2 concentrations. 
# It produces a CSV file of daily CO2 for all modeled sites. 
# It produces a CSV file of yearly CO2 for all modeled sites. (optional)

rm(list=ls())

# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'

#
library(librarian)

shelf(ggplot2, dplyr, tidyr, cowplot, ggpubr, trend, amerifluxr, lme4)

######-------------------------------------------------------------------------------
# Look at globally CO2 dataset first. 
sites_long <- read.csv(file.path('data', 'long_term_ec_sites.csv'))
rm_years   <- read.csv(file.path('data', 'years_to_remove.csv'))

source(file.path("scripts", "00_functions.R"))
df_file_structure <- get_rawflux_file_structure(dir_rawdata)

a_YY_all <- data.frame(site_ID=character(), TIMESTAMP=integer(), NEE_VUT_REF=double(), 
                       GPP_VUT_REF=double(), RECO_VUT_REF=double(), CO2_F_MDS=double())
#
sites_long$nyear <- 0
#
for (i in 1:nrow(sites_long)) {
  # i = 20
  print(i)
  site_ID <- sites_long$site_ID[i]
  print(site_ID)
  folder <- unlist(df_file_structure$folders[df_file_structure$sources == sites_long$source[i]])
  YYfile <- df_file_structure$YYfiles[df_file_structure$sources == sites_long$source[i]]
  a_YY <- read.csv(folder[grepl(pattern=paste0(site_ID, YYfile), folder)])
  a_YY[a_YY==-9999] <- NA
  #
  a_YY_good <- a_YY %>% filter(!TIMESTAMP %in% rm_years$years[rm_years$site_ID==site_ID])
  sites_long$nyear[i] <- sum(!is.na(a_YY_good$NEE_VUT_REF))
  #
  # plot(a_YY_good$TIMESTAMP, a_YY_good$NEE_VUT_REF, main=site_ID)
  # plot(a_YY_good$TIMESTAMP, a_YY_good$GPP_NT_VUT_REF, main=site_ID)
  # plot(a_YY_good$TIMESTAMP, a_YY_good$RECO_NT_VUT_REF, main=site_ID)
  
  if (sites_long$NEEpartition[i] == 'DT' & "CO2_F_MDS" %in% names(a_YY_good)) {
    df <- a_YY_good[, c("TIMESTAMP", "NEE_VUT_REF", "GPP_DT_VUT_REF", "RECO_DT_VUT_REF", "CO2_F_MDS")]
  } else if (sites_long$NEEpartition[i] == 'NT' & "CO2_F_MDS" %in% names(a_YY_good)) {
    df <- a_YY_good[, c("TIMESTAMP", "NEE_VUT_REF", "GPP_NT_VUT_REF", "RECO_NT_VUT_REF", "CO2_F_MDS")]
  }
  names(df)[3:4] <- c("GPP_VUT_REF", "RECO_VUT_REF")
  # remove unreasonable GPP and RECO. 
  df[df$GPP_VUT_REF < 10 | is.na(df$GPP_VUT_REF), 3:4] <- NA
  #
  a_YY_all <- rbind(a_YY_all, data.frame(site_ID=site_ID, df))
}

# Adding IGBP and categories
if (!'category' %in% names(a_YY_all)) {
  a_YY_all <- a_YY_all %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category")], by = "site_ID")
}

a_YY_all_decade <- a_YY_all %>% filter(between(TIMESTAMP,1990, 2000))
summary(lmer(data = a_YY_all_decade, CO2_F_MDS ~ TIMESTAMP + (1 | site_ID)))

a_YY_all %>% # filter(site_ID %in% c("BE-Vie", "CH-Dav", "DE-Tha", "FI-Hyy", "US-Ha1")) %>% 
  ggplot(aes(x=TIMESTAMP, y=CO2_F_MDS, color=site_ID)) +
  geom_point() +
  geom_line() +
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme_classic() + 
  theme(legend.position = "none")


a_YY_all %>% # filter(site_ID %in% c("BE-Vie", "CH-Dav", "DE-Tha", "FI-Hyy", "US-Ha1")) %>% 
  ggplot(aes(x=TIMESTAMP, y=CO2_F_MDS)) +
  geom_point() +
  geom_line() +
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme(legend.position = "none")

# Read some AmeriFlux data
# US-Cwt: we only have abnormal CO2 concentration; suggested range ~410 – 440 ppm on average; diurnal variation (~±5–10 ppm); 
# US-SP3: no measured data for this site, but we can use US-SP1 (2000-2020), and US-SP2 (1999-2008). 
name_site <-     
data_BASE <- amf_read_base('/Users/junnawang/Downloads/drive-download-20260204T075928Z-3-001/AMF_US-SP3_BASE-BADM_3-1.zip', parse_timestamp=TRUE, unzip = T)
plot(data_BASE$TIMESTAMP, data_BASE$CO2_2)

