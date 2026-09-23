# An Rscript to compile modeled surface soil water content for each site
# 1/28/2026, Junna Wang

# important: SWC_ERA5 does not include data from IT-Noe, used measured data for this site.
# it takes ~40 min to complete. 
rm(list=ls())

# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'
#
dir_rawdata <- '/Volumes/WZZ_disk/Carbon_trend/SiteData'

# Find all the sites with FLUXNET products with good measured SWC values
sites_long    <- read.csv(file.path('data2026', 'long_term_ec_sites.csv'))
# How many sub time series do we have?
rm_years   <- read.csv(file.path('data2026', 'years_to_remove.csv'))

source(file.path("scripts2026", "00_functions.R"))
df_file_structure <- get_rawflux_file_structure(dir_rawdata)

# read analyzed SWC
SWC_ERA5 <- read.csv(file.path(dir_rawdata, "ERA5L_daily_SWC_trend_allsites_1990_2025_download.csv"))
names(SWC_ERA5)[3] <- 'SWC'

# terra soil water data
# SWC_ERA5 <- read.csv(file.path(dir_rawdata, "terra_SWC_allsites.csv"))
# names(SWC_ERA5) <- c("date", "SWC", "name")

SWC_ERA5$date <- ymd(SWC_ERA5$date)
SWC_ERA5$year <- year(SWC_ERA5$date) 
sites_SWC_ERA5 <- unique(SWC_ERA5$name)

for (i in 1:nrow(sites_long)) {
  site_ID <- sites_long$site_ID[i]
  if (!site_ID %in% "IT-Noe") next
  #
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
  year_start <- min(a_YY_good$TIMESTAMP, na.rm = T)
  year_end <- max(a_YY_good$TIMESTAMP, na.rm = T)
  measured_SWC_availabe <- TRUE
  
  if ("SWC_F_MDS_1" %in% names(a_DD_good) & sum(!is.na(a_YY_good$SWC_F_MDS_1)) > 5) {
    # the unit of SWC is inconsistent; make them consistent
    if (mean(a_YY_good$SWC_F_MDS_1, na.rm=T) < 1) {
      a_DD_good$SWC_F_MDS_1 <- a_DD_good$SWC_F_MDS_1 * 100
      a_YY_good$SWC_F_MDS_1 <- a_YY_good$SWC_F_MDS_1 * 100
    }
    #
    print("SWC_F_MDS_1")
    plot(a_DD_good$TIMESTAMP, a_DD_good$SWC_F_MDS_1, main = site_ID)
    if (site_ID == "IT-Noe") {
      df = data.frame(name = site_ID, date = a_DD_good$TIMESTAMP, SWC = a_DD_good$SWC_F_MDS_1)
    }
  }
}

# One site IT-Noe, their modeled SWC is missing, we will use measured data. 
SWC_ERA5 <- SWC_ERA5 %>% 
  mutate(SWC = SWC * 100) %>% 
  select("date", "name", "SWC") %>% 
  filter(!name == "IT-Noe") %>% 
  bind_rows(df)

# note: remotely sendsed SWC at GL-HaZ is constant, using measured data in script 02_05. 

# add CH-Oe3 using CH-Oe1, and add ES-LM1 using ES-LM2; because these pairs are really close (<10km)
tmp1 <- SWC_ERA5 %>% filter(name == "CH-Oe1") %>% mutate(name = "CH-Oe3")
tmp2 <- SWC_ERA5 %>% filter(name == "ES-LM2") %>% mutate(name = "ES-LM1")
SWC_ERA5 <- bind_rows(SWC_ERA5, tmp1, tmp2)

write.csv(SWC_ERA5, file.path(dir_rawdata, "ERA5L_daily_SWC_trend_allsites_1990_2025_correct.csv"), row.names = F)

