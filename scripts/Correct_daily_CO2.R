# This script is used to correct sites with abnormal CO2 concentrations.
# It produces a CSV file of yearly CO2 for all modeled sites. (We will do this)
# Junna Wang, Feb 10, 2026

rm(list=ls())

# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'

##
library(librarian)
shelf(ggplot2, dplyr, tidyr, cowplot, ggpubr, trend, amerifluxr, lme4, ncdf4, lubridate, tibble, stringr)

######Step 1: Look at CO2 data in our flux data-set-------------------------------------------------------------------------------
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

a_YY_all %>% 
  ggplot(aes(x=TIMESTAMP, y=CO2_F_MDS, color=site_ID)) +
  geom_point() +
  geom_line() +
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme_classic() + 
  theme(legend.position = "none")

write.csv(a_YY_all[, 1:6], file = file.path(data, "a_YY_all.csv"), row.names = F)
######end of Step 1-------------------------------------------------------------------------------

######Step 2: pick up site-years with good CO2 measurements from our flux datasets---------------- 

# Read some AmeriFlux data to confirm problematic CO2 data. 
# US-Cwt: we only have abnormal CO2 concentration; suggested range ~410 – 440 ppm on average; diurnal variation (~±5–10 ppm); 
# US-SP3: no measured data for this site, but we can use US-SP1 (2000-2020), and US-SP2 (1999-2008). 
# data_BASE <- amf_read_base('/Users/junnawang/Downloads/drive-download-20260204T075928Z-3-001/AMF_US-SP3_BASE-BADM_3-1.zip', parse_timestamp=TRUE, unzip = T)
# plot(data_BASE$TIMESTAMP, data_BASE$CO2_2)

a_YY_all <- read.csv(file = file.path(data, "a_YY_all.csv"))

# 29 sites with good CO2
flux_sites_good_CO2 <- c('CA-Ca1', 'CA-Oas', 'CA-Obs', "CH-Fru", "DE-Obe", "DE-RuR", "DE-Tha", "FI-Hyy",
"FI-Let", "FI-Sod", "FR-FBn", "IT-Ro2", "JP-Api", "JP-Tak", "US-GLE", "US-KL3", "US-KL2", "US-KM4",
"US-MOz", "US-PFa", "US-Rws", "US-Seg", "US-SRG", "US-SRM", "US-Tw1", "US-UMB", "US-UMd", "US-Whs", "US-Wjs")
flux_sites_partial_good_CO2 <- tribble(
  ~site_ID, ~years2remove,
  "BR-Sa1", c(2002:2005, 2009:2011, 2016:2019), 
  "CA-Cbo", 1996, 
  "CA-Gro", 2008, 
  "CA-KP1", 2008, 
  "CA-Man", 2007, 
  "CA-Mer", 1999, 
  "CA-Ojp", c(2002, 2003), 
  "CA-TP4", 2003, 
  "CA-TPD", 2012, 
  "CZ-wet", c(2008:2010), 
  "DE-Gri", 2004,
  "DE-Hai", c(2013:2025),
  "ES-LJu", c(2017), 
  "GF-Guy", c(2017:2025),
  "IL-Yat", 2002,
  "JP-Fhk", c(2008, 2009, 2020), 
  "JP-Fjy", 2020, 
  "JP-Spp", c(2002, 2017), 
  "JP-Ta2", 2010,
  "JP-Tef", c(2021, 2022, 2023),  
  "RU-Fyo", 2017,
  "US-ChR", c(2015:2025),
  "US-Ha1", c(2006, 2007),
  "US-Ho1", c(2008, 2009, 2010), 
  "US-ICh", c(2021:2025), 
  "US-ICs", c(2013:2015, 2021), 
  "US-ICt", 2021, 
  "US-Jo2", 2012, 
  "US-KM3", 2012, 
  "US-Kon", 2007, 
  "US-LL1", c(2012, 2015), 
  "US-Los", 2025, 
  "US-MMS", c(1999, 2002), 
  "US-Mpj", c(2014:2016, 2018:2020), 
  "US-MtB", c(2012, 2013, 2021), 
  "US-NC2", 2009,
  "US-NR1", c(2004, 2005), 
  "US-Oho", c(2004, 2005), 
  "US-Ses", c(2011, 2016:2018), 
  "US-Syv", c(2000:2010, 2019), 
  "US-Ton", c(2004, 2005), 
  "US-Var", 2004,
  "US-Vcm", 2017, 
  "US-Vcp", c(2016, 2017), 
  "US-WCr", c(2021, 2025), 
  "US-Wjs", 2011, 
  "US-Wkg", 2007,
  "AU-Cpr", c(2014, 2018:2020), 
  "AU-Dry", c(2010:2016),
  "AU-Stp", c(2011, 2018), 
  "AU-Tum", c(2005, 2008)
)

# In the end, every site-year in a_YY_all has to be filled with values. 

a_YY_co2_good <- a_YY_all[, c("site_ID", "TIMESTAMP", "CO2_F_MDS")] %>% 
  filter(site_ID %in% flux_sites_good_CO2) %>% 
  filter(!is.na(CO2_F_MDS)) 

a_YY_all[, c("site_ID", "TIMESTAMP", "CO2_F_MDS")]
for (i in 1: nrow(flux_sites_partial_good_CO2)) {
  #
  print(i)
  df <- a_YY_all[, c("site_ID", "TIMESTAMP", "CO2_F_MDS")] %>% filter(site_ID == flux_sites_partial_good_CO2$site_ID[i]) %>%
    filter(! TIMESTAMP %in% flux_sites_partial_good_CO2$years2remove[i]) %>%
    filter(!is.na(CO2_F_MDS))
  a_YY_co2_good <- bind_rows(a_YY_co2_good, df)
}

a_YY_co2_good <- a_YY_co2_good %>% left_join(sites_long[, c("site_ID", "Lat", "Long")], by = "site_ID") 
#----------------------------------------------
# create a table for ObsPack metadata 
# site, country, lat, long, start, end, measurement

read_obspack_metadata <- function(file) {
  # file <- files[1]
  print(file)
  # file = "/Users/junnawang/YaleLab/data_CO2/obspack_co2_1_GLOBALVIEWplus_v11.0_2026-01-12/metadata/co2_afo_aircraft-insitu_1_allvalid.txt"
  
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

files <- list.files('/Users/junnawang/YaleLab/data_CO2/obspack_co2_1_GLOBALVIEWplus_v11.0_2026-01-12_daily/metadata', pattern = '.txt', full.names = T)
all_sites_info <- do.call(rbind, lapply(files, read_obspack_metadata))

# ~ 300 stations have >10 years of data
# ~ 128 stations have >20 years of data
all_sites_info_10yr <- all_sites_info %>% filter(data_duration_years > 10) %>% filter(site_latitude > - 100) %>% 
  filter(dataset_project %in% c("surface-insitu", "tower-insitu", "shipboard-insitu"))
all_sites_info_10yr$magl <- as.numeric(str_extract(all_sites_info_10yr$dataset_name, "(?<=-)[0-9]+(?=magl)"))

# at one site, only choose the one at around canopy ~ 30 magl
# remove these duplicated datasets for one location, only save one dataset closest to 30 magl
all_sites_info_10yr_unique <- all_sites_info_10yr %>% 
  group_by(site_name, site_latitude, site_longitude) %>%
  slice_min(abs(magl - 30), n = 1) %>%  # keeps row with layer_height closest to 30
  # arrange(abs(magl - 30)) %>%
  # slice(1) %>%
  ungroup()

# take a look at how CO2 concentration change over elevation; the magnitude of CO2 variation should decline with elevation
# the magnitude of CO2 concentration can be so high during winter time. This is some information we can use. 
df_obspack_co2_yearly <- data.frame(dataset_name = character(), year = integer(), co2 = double())
dir <- "/Users/junnawang/YaleLab/data_CO2/obspack_co2_1_GLOBALVIEWplus_v11.0_2026-01-12_nc/data/nc/"

for (dataset in all_sites_info_10yr_unique$dataset_name) {
  # the two sites co2 is higher than other sites. both are in Utah state. not sure why. 
  # co2_utrpk_tower-insitu_432_allvalid and co2_utsug_tower-insitu_432_allvalid
  # if (!dataset %in% c("co2_utrpk_tower-insitu_432_allvalid", "co2_utsug_tower-insitu_432_allvalid")) next
  file <- paste0(dataset, '.nc') 
  nc_file <- nc_open(file.path(dir, file))
  # plot(ymd_hms(ncvar_get(nc_file, "datetime")), ncvar_get(nc_file, "value") * 1000000, main = dataset)
  df <- data.frame(datetime = ymd_hms(ncvar_get(nc_file, "datetime")), value = ncvar_get(nc_file, "value") * 1000000)
  df <- df %>% mutate(year = year(datetime)) %>% group_by(year) %>%
    summarise(co2 = mean(value)) %>% ungroup %>%
    mutate(dataset_name = dataset)
  #  
  df_obspack_co2_yearly <- bind_rows(df_obspack_co2_yearly, df)
}

# data of this site year is exceedingly high: co2_cvo_surface-insitu_456_allvalid (2016, 444.7239); 
df_obspack_co2_yearly <- df_obspack_co2_yearly %>% filter(! (dataset_name == 'co2_cvo_surface-insitu_456_allvalid' & year == 2016))

# take a look at the dataset
df_obspack_co2_yearly %>% 
  ggplot(aes(x=year, y=co2, color=dataset_name)) +
  geom_point() +
  geom_line() +
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme_classic() + 
  theme(legend.position = "none")

# the quality of these data is high. 
# Sable Island, Nova Scotia; co2_wsa_surface; year 2022 may be problematic, but those are quality-checked-data.  

# I need to create an annual database it stores data_source, lat, long, year, co2 values, 
# plot all the figures together, and see the similarity of these patterns. 
# 
