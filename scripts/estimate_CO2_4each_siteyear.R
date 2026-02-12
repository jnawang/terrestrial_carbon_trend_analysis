# This script is used to correct sites with abnormal CO2 concentrations.
# It produces a CSV file of yearly CO2 for all modeled sites. (We will do this)
# Junna Wang, Feb 10, 2026
# This takes ~10 mins

rm(list=ls())

# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'

##
library(librarian)
shelf(ggplot2, dplyr, tidyr, cowplot, ggpubr, trend, amerifluxr, lme4, ncdf4, lubridate, tibble, stringr, zoo, geosphere, rstan)

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

write.csv(a_YY_all[, 1:6], file = file.path("data", "a_YY_all.csv"), row.names = F)
######end of Step 1-------------------------------------------------------------------------------

######Step 2: pick up site-years with good CO2 measurements from our flux datasets---------------- 

# Read some AmeriFlux data to confirm problematic CO2 data. 
# US-Cwt: we only have abnormal CO2 concentration; suggested range ~410 – 440 ppm on average; diurnal variation (~±5–10 ppm); 
# US-SP3: no measured data for this site, but we can use US-SP1 (2000-2020), and US-SP2 (1999-2008). 
# data_BASE <- amf_read_base('/Users/junnawang/Downloads/drive-download-20260204T075928Z-3-001/AMF_US-SP3_BASE-BADM_3-1.zip', parse_timestamp=TRUE, unzip = T)
# plot(data_BASE$TIMESTAMP, data_BASE$CO2_2)

a_YY_all <- read.csv(file = file.path("data", "a_YY_all.csv"))

# 29 sites with good CO2
flux_sites_good_CO2 <- c('CA-Ca1', 'CA-Oas', 'CA-Obs', "CH-Fru", "DE-Obe", "DE-RuR", "DE-Tha", "FI-Hyy",
"FI-Let", "FI-Sod", "FR-FBn", "IT-Ro2", "JP-Api", "JP-Tak", "US-GLE", "US-KL3", "US-KL2", "US-KM4",
"US-MOz", "US-PFa", "US-Rws", "US-Seg", "US-SRG", "US-SRM", "US-Tw1", "US-UMB", "US-UMd", "US-Whs", "US-Wjs")
flux_sites_partial_good_CO2 <- tribble(
  ~site_ID, ~years2remove, 
  "BR-Sa1", c(2020, 2021), 
  "CA-Cbo", 1996, 
  "CA-Gro", 2008, 
  "CA-KP1", 2008, 
  "CA-Man", 2007, 
  "CA-Mer", 1999, 
  "CA-Ojp", c(2002, 2003), 
  "CA-TP4", 2003, 
  "CA-TPD", 2012, 
  # "CZ-wet", c(2008:2010),     # CO2 is super high and not sure if it is realistic. 
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

for (i in 1: nrow(flux_sites_partial_good_CO2)) {
  #
  print(i)
  df <- a_YY_all[, c("site_ID", "TIMESTAMP", "CO2_F_MDS")] %>% filter(site_ID == flux_sites_partial_good_CO2$site_ID[i]) %>%
    filter(! TIMESTAMP %in% unlist(flux_sites_partial_good_CO2$years2remove[i])) %>%
    filter(!is.na(CO2_F_MDS))
  a_YY_co2_good <- bind_rows(a_YY_co2_good, df)
}

a_YY_co2_good <- a_YY_co2_good %>% left_join(sites_long[, c("site_ID", "Lat", "Long")], by = "site_ID") 

# take a look at the data
a_YY_co2_good %>% 
  ggplot(aes(x=TIMESTAMP, y=CO2_F_MDS, color=site_ID)) +
  geom_point() +
  geom_line() +
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme_classic() + 
  theme(legend.position = "none")

a_YY_co2_good %>% filter(CO2_F_MDS < 350)

######end of Step 2------------------------------------------------------------------------ 

######Step 3: get yearly CO2 data from ObsPack---------------------------------------------

# create a table for ObsPack metadata including site, country, lat, long, start, end, measurement
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

files <- list.files('/Users/junnawang/YaleLab/data_CO2/obspack_co2_1_GLOBALVIEWplus_v11.0_2026-01-12_daily/metadata', pattern = '.txt', full.names = T)
all_sites_info <- do.call(rbind, lapply(files, read_obspack_metadata))

# ~ 300 stations have >10 years of data
# ~ 128 stations have >20 years of data
all_sites_info_10yr <- all_sites_info %>% filter(data_duration_years >= 10) %>% filter(site_latitude > - 100) %>% 
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

# take a look at the dataset: the quality of these data is high. 
df_obspack_co2_yearly %>% 
  ggplot(aes(x=year, y=co2, color=dataset_name)) +
  geom_point() +
  geom_line() +
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme_classic() + 
  theme(legend.position = "none")

# Sable Island, Nova Scotia; co2_wsa_surface; year 2022 may be problematic, but those are quality-checked-data.  
######end of Step 3---------------------------------------------------------------------------

######Step 4: combine ObsPack and a_YY_co2_good---------------------------------------------
# check if there are duplicated sites between the two datasets
dup_pairs <- sites_long %>% filter(site_ID %in% c(flux_sites_good_CO2, flux_sites_partial_good_CO2$site_ID)) %>%
  select(c("site_ID", "Lat", "Long")) %>%
  mutate(id1 = row_number()) %>%
  inner_join(
    all_sites_info_10yr_unique[, c("site_latitude", "site_longitude")] %>% mutate(id2 = row_number()),
    by = character() ) %>%
  filter(abs(Lat - site_latitude) <= 1e-3, abs(Long - site_longitude) <= 1e-3) %>%
  pull(site_ID)

# US-PFa is the duplicated site; remove it from the good flux data
a_YY_co2_good <- a_YY_co2_good %>% filter(site_ID != dup_pairs)

# combine the two data
df_obspack_ec_co2_yearly <- df_obspack_co2_yearly %>% left_join(all_sites_info_10yr_unique[, c('dataset_name', 'site_latitude', 'site_longitude')], by = "dataset_name") %>%
  select(-dataset_name) %>% bind_rows(a_YY_co2_good[, c("TIMESTAMP", "CO2_F_MDS", "Lat", "Long")] %>% 
  rename("year"="TIMESTAMP", "co2"="CO2_F_MDS", "site_latitude"='Lat', "site_longitude"="Long")) %>%
  mutate(site = paste0(round(site_latitude, 4), '-', round(site_longitude, 4)))

# plot all the data together to see the similarity of these patterns. 
df_obspack_ec_co2_yearly %>% 
  ggplot(aes(x=year, y=co2, color = site)) +
  geom_point() +
  geom_line() + 
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme_classic() + 
  theme(legend.position = "none")

######look at CO2 trends in different regions and periods
subdata <- df_obspack_ec_co2_yearly %>% 
  filter(year < 2026 & year >= 2020) %>% filter(site_latitude > 35 & site_latitude < 50)
  
# site_latitude > -90 & site_latitude < -23.5
# site_latitude < 23.5 & site_latitude > -23.5
# site_latitude > 23.5 & site_latitude < 35
# site_latitude > 35 & site_latitude < 50
# site_latitude > 50 & site_latitude < 90

summary(lmer(data = subdata, co2 ~ year + (1 | site)))
######end of Step 4-------------------------------------------------------------------------

######Step 5: estimating CO2 for every site-year--------------------------------------------
baseline_lat_long <- df_obspack_ec_co2_yearly %>%
  count(site, site_latitude, site_longitude, name = "nyear")
# 8 baseline sites with >= 35-year data
id_baseline_longterm <- which(baseline_lat_long$nyear >= 35)

# One function to calculate distance between two locations. 
dist_km <- function(lat0, lon0, lat, lon) {
  distHaversine(
    matrix(c(lon0, lat0), ncol = 2),
    cbind(lon, lat)
  ) / 1000
}

# Rstan model to estimate CO2 at missing site-years
stan_code <- "
data {
  int<lower=1> T;
  int<lower=1> P;
  int<lower=1> N_obs;

  int<lower=1,upper=T> obs_t[N_obs];
  vector[N_obs] y_obs;

  int<lower=1,upper=P> period[T];

  vector[P] beta_mu;
  vector[P] beta_se;
}
parameters {
  real y0;                 // initial absolute level
  vector[P] beta;          
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[T] y;

  y[1] = y0;

  for (t in 2:T)
    y[t] = y[t-1] + beta[period[t]];
}
model {
  // anchor absolute level with prior
  y0 ~ normal(355, 8.4);

  beta ~ normal(beta_mu, beta_se);
  sigma_obs ~ exponential(1);

  y_obs ~ normal(y[obs_t], sigma_obs);
}
"

years <- 1990:2025
T <- length(years)

period <- cut(
  years,
  breaks = c(1990, 2000, 2010, 2020, 2025),
  labels = 1:4,
  # labels = paste0('p', 1:4),
  right = TRUE, 
  include.lowest = TRUE
) |> as.integer()

# End of information for Rstan model. 

# A data frame to store gap-filled CO2 data
a_YY_co2_gap_filled <- a_YY_all[, c("site_ID", "TIMESTAMP", "CO2_F_MDS")]
a_YY_co2_gap_filled$CO2_F_MDS[!a_YY_co2_gap_filled$site_ID %in% c(flux_sites_good_CO2, flux_sites_partial_good_CO2$site_ID)] <- NA
for (i in 1: nrow(flux_sites_partial_good_CO2)) {
  # print(i)
  a_YY_co2_gap_filled$CO2_F_MDS[a_YY_co2_gap_filled$site_ID == flux_sites_partial_good_CO2$site_ID[i] & a_YY_co2_gap_filled$TIMESTAMP %in% unlist(flux_sites_partial_good_CO2$years2remove[i]) ] <- NA
}

sites_long$dist_co2_proxy <- 0.0
# gap-fill missed annual CO2
for (i in 1:nrow(sites_long)) {
  fluxsite <- sites_long$site_ID[i]
  a_YY_co2_site <- a_YY_co2_gap_filled %>% filter(site_ID == fluxsite)
  if (fluxsite %in% c(flux_sites_good_CO2, flux_sites_partial_good_CO2$site_ID)) {
    print(plot(a_YY_co2_site$TIMESTAMP, a_YY_co2_site$CO2_F_MDS, main = paste0('before--', fluxsite)))
  }
  
  if (fluxsite %in% flux_sites_good_CO2) {
    # for sites with good measurements, just use interpolation to gap-fill NA
    if ( any(is.na(a_YY_co2_site$CO2_F_MDS)) ) {
      a_YY_co2_site$CO2_F_MDS <- na.approx(a_YY_co2_site$CO2_F_MDS, x = a_YY_co2_site$TIMESTAMP, na.rm = FALSE)
      if (is.na(a_YY_co2_site$CO2_F_MDS[1])) {
        # use the first 10 values
        mod <- lm(data = a_YY_co2_site[1:10, ], CO2_F_MDS ~ TIMESTAMP)
        first_non_na <- match(TRUE, !is.na(a_YY_co2_site$CO2_F_MDS))
        a_YY_co2_site$CO2_F_MDS[1:(first_non_na - 1)] <- predict(mod, newdata = data.frame(TIMESTAMP = a_YY_co2_site$TIMESTAMP[1:(first_non_na - 1)]))
      }
      if (is.na(a_YY_co2_site$CO2_F_MDS[nrow(a_YY_co2_site)])) {
        # use the last 10 values
        id_last <- nrow(a_YY_co2_site)
        mod <- lm(data = a_YY_co2_site[(id_last-9):id_last, ], CO2_F_MDS ~ TIMESTAMP)
        last_non_na <- id_last - match(TRUE, rev(!is.na(a_YY_co2_site$CO2_F_MDS))) + 1
        a_YY_co2_site$CO2_F_MDS[(last_non_na + 1):id_last] <- predict(mod, newdata = data.frame(TIMESTAMP = a_YY_co2_site$TIMESTAMP[(last_non_na + 1):id_last]))
      }
    }
    
  } else {
    # find the nearest site
    dist_km_point <- dist_km(
      lat0 = sites_long$Lat[i], lon0 = sites_long$Long[i],
      lat  = baseline_lat_long$site_latitude,
      lon  = baseline_lat_long$site_longitude
    )
    sites_long$dist_co2_proxy[i] <- min(dist_km_point)
    idx <- which.min(dist_km_point)
    # print(paste0(fluxsite, '--', dist_km_point[idx], '--', baseline_lat_long$site_latitude[idx], '--', baseline_lat_long$site_longitude[idx]) )

    # get observation data
    obs <- df_obspack_ec_co2_yearly %>% filter(site == baseline_lat_long$site[idx]) %>% 
      select(c(year, co2)) %>% rename('value' = 'co2') %>% filter(year >= 1990)
    obs_t <- match(obs$year, years)
    
    # get trend information
    # the default trend
    beta_mu <- c(1.81, 2.04, 2.36, 2.85)
    beta_se <- c(0.5, 0.5, 0.6, 0.6)    
    # update trend using the period trends of nearby 20 sites
    idx <- order(dist_km_point)[1:20]
    # it has to include one site with > 35 years of data to cover the trend from 1990 to 2025
    if (max(baseline_lat_long$nyear[idx]) < 35) {
      idx <- c(idx, id_baseline_longterm[which.min(dist_km_point[id_baseline_longterm])])
    }
    
    subdata <- df_obspack_ec_co2_yearly %>% filter(site %in% baseline_lat_long$site[idx])
    for (j in unique(period)) {
      subdata_period <- subdata %>% filter(year %in% years[period == j])
      if (length(unique(subdata_period$site)) > 1) {
        mod <- lmer(data = subdata_period, co2 ~ year + (1 | site))
        beta_mu[j] = fixef(mod)[2]        
      } else {
        mod <- lm(data = subdata_period, co2 ~ year)
        beta_mu[j] = coef(mod)[2]        
      }
    }
    
    # get stan data
    stan_data <- list(
      T = T,
      P = 4,
      N_obs = nrow(obs),
      obs_t = obs_t,
      y_obs = obs$value,
      period = period,
      beta_mu = beta_mu,
      beta_se = beta_se
    )
    
    # run the model
    fit <- stan(
      model_code = stan_code,
      data = stan_data,
      chains = 4,
      iter = 2000,
      control = list(
        adapt_delta = 0.99,
        max_treedepth = 15
      )
    )
    # summary(fit)
    
    # get y values
    post <- rstan::extract(fit)
    y_mean <- apply(post$y, 2, mean)

    # gap fill 
    id_na <- which(is.na(a_YY_co2_site$CO2_F_MDS))
    a_YY_co2_site$CO2_F_MDS[id_na] <- y_mean[match(a_YY_co2_site$TIMESTAMP[id_na], years)]
  }
  print(plot(a_YY_co2_site$TIMESTAMP, a_YY_co2_site$CO2_F_MDS, main = paste0('after--', fluxsite)))
  
  a_YY_co2_gap_filled <- a_YY_co2_gap_filled %>%
    rows_update(a_YY_co2_site, by = c("site_ID", "TIMESTAMP"))
}

# ~ 10 sites, we have to use data about 2000 km away. 
hist(sites_long$dist_co2_proxy)
# if use >10-year data from obspack, 2 sites use >2000km (max: 2472); 8 sites use >1000km; 15 sites use > 500km; 29 sites use > 200km; 47 sites use > 100km;
# if use >8-year data from obspack, 2 sites use >2000km (max: 2472); 6 sites use >1000km; 13 sites use > 500km; 25 sites use > 200km; 46 sites use > 100km;
# small differences between them, so use 10-year as a criteria. 


# take a look at the data before we finally save it.
a_YY_co2_gap_filled %>% 
  ggplot(aes(x=TIMESTAMP, y=CO2_F_MDS, color = site_ID)) +
  geom_point() +
  geom_line() + 
  labs(x='Year', y='CO2 concentration (ppm)') + 
  theme_classic() + 
  theme(legend.position = "none")

# save the new dataset
write.csv(a_YY_co2_gap_filled, file = file.path("data", "a_YY_co2_gap_filled.csv"), row.names = F)

######end of Step 5----------------------------------------------------------------------------------