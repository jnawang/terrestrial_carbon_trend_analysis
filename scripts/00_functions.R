
# This function organizes different sources of raw flux data; updated on Jan, 2026
# It takes a few mins to run this function
get_rawflux_file_structure <- function(dir_rawdata){
  
  files.unzip.Ameri <- list.files(file.path(dir_rawdata, "Ameriflux_FLUXNET", "unzip"), pattern =  "_FLUXNET_FULLSET_(HR|HH|DD|YY)", full.names = T)
  files.unzip.FLUXNET2020 <- list.files(file.path(dir_rawdata, "FLUXNET2020", "unzip"), pattern = "_FLUXNET2015_FULLSET_(HR|HH|DD|YY)", full.names = T)
  files.unzip.FLUXNET2015 <- list.files(file.path(dir_rawdata, "FLUXNET2015_10yr_only"), pattern = "", full.names = T)
  files.unzip.FLUXNET2020_ICOS <- list.files(file.path(dir_rawdata, "ICOS_FLUXNET2025-1", "unzip_connect2020"), full.names = T)
  files.unzip.FLUXNET2015_ICOS <- list.files(file.path(dir_rawdata, "ICOS_FLUXNET2025-1", "unzip_connect2015"), full.names = T)
  files.Ameri.ReddyProcGapFill <- list.files(file.path(dir_rawdata, "Ameriflux_BASE", "ReddyProcGapFill"), full.names=T)
  files.unzip.JapanFlux2024 <- list.files(file.path(dir_rawdata, "JapanFlux2024", "unzip"), full.names=T)
  files.csv.OzFlux2025 <- list.files(file.path(dir_rawdata, "OzFlux2025", "csv"), full.names=T)
  files.csv.HRBFlux2022 <- list.files(file.path(dir_rawdata, "HRBFlux2022", "csv"), full.names=T)
  #
  df_file_structure <- data.frame(sources = c("AmeriFlux_FLUXNET", "AmeriFlux_BASE", 
                                              "FLUXNET2020_ICOS", "FLUXNET2015_ICOS", "FLUXNET2020", "FLUXNET2015", 
                                              "JapanFlux2024", "OzFlux2025", "HRBFlux2022"), 
                                  folders = I(list(list(files.unzip.Ameri), list(files.Ameri.ReddyProcGapFill), 
                                                   list(files.unzip.FLUXNET2020_ICOS), list(files.unzip.FLUXNET2015_ICOS),
                                                   list(files.unzip.FLUXNET2020), list(files.unzip.FLUXNET2015), 
                                                   list(files.unzip.JapanFlux2024), list(files.csv.OzFlux2025), 
                                                   list(files.csv.HRBFlux2022))), 
                                  HHfiles = c("_FLUXNET_FULLSET_(HR|HH)", "_EDDYPROC_FULLSET_HH",
                                              "_FLUXNET_HH", "_FLUXNET_HH", 
                                              "_FLUXNET2015_FULLSET_(HR|HH)", "_FLUXNET2015_FULLSET_(HR|HH)", 
                                              "_JapanFLUX2024_ALLVARS_HH", "_PyFluxPro_FULLSET_HH",
                                              "_ReddyProcRF_FULLSET_HH"), 
                                  DDfiles = c("_FLUXNET_FULLSET_DD", "_EDDYPROC_FULLSET_DD",
                                              "_FLUXNET_DD", "_FLUXNET_DD", 
                                              "_FLUXNET2015_FULLSET_DD", "_FLUXNET2015_FULLSET_DD", 
                                              "_JapanFLUX2024_ALLVARS_DD", "_PyFluxPro_FULLSET_DD",
                                              "_ReddyProcRF_FULLSET_DD"),
                                  YYfiles = c("_FLUXNET_FULLSET_YY", "_EDDYPROC_FULLSET_YY", 
                                              "_FLUXNET_YY", "_FLUXNET_YY", 
                                              "_FLUXNET2015_FULLSET_YY", "_FLUXNET2015_FULLSET_YY",
                                              "_JapanFLUX2024_ALLVARS_YY", "_PyFluxPro_FULLSET_YY",
                                              "_ReddyProcRF_FULLSET_YY"))
  return(df_file_structure)
}

#---------------functions to calculate absolute or relative slopes--------------
# these are old methods ignoring irregular time series
# cal_sens_slope <- function(x) {
#   x <- x[!is.na(x)]
#   return(as.numeric(sens.slope(x)[c(1)]))
# }
# 
# cal_sens_slope_scale <- function(x) {
#   x <- x[!is.na(x)]
#   x <- scale(x)
#   return(as.numeric(sens.slope(x)[c(1,3)]))
# }

# use new methods covering irregular time series

library(mblm)
cal_sens_slope <- function(time, y) {
  id.nonna <- !is.na(y) & !is.na(time) 
  y <- y[id.nonna]
  time <- time[id.nonna]
  #
  model <- mblm(y ~ time)
  # summary(model)
  return(summary(model)$coefficients[2, 1])
}

cal_sens_slope_scale <- function(time, y, yall) {
  id.nonna <- !is.na(y) & !is.na(time) 
  y <- y[id.nonna]
  time <- time[id.nonna]
  #
  if (is.null(yall)) {
    ymean <- mean(y)
    ystd <- sd(y)
  } else {
    ymean <- mean(yall, na.rm=T)
    ystd <- sd(yall, na.rm=T)
  }
  # scale y
  y <- (y - ymean) / ystd
  model <- mblm(y ~ time)
  # summary(model)
  return(as.numeric(summary(model)$coefficients[2,c(1, 4)]))
}


sites_grow_yearlong <- c('BR-Ma2', 'BR-Sa1', 'US-Elm', 'US-Esm', 'US-LL1', 'US-SP1', 'US-SP3', 
                         'GF-Guy', 'ID-PaD', 'AU-How', 'AU-Rob', "AU-Tum", "AU-Cum", "AU-Cow")

# match climate class and its coding ID
df_koppen_climate <- tribble(
  ~Climate.class, ~climate_id,
  "Af", 1,  "Am", 2,  "Aw", 3,
  "BWh", 4, "BWk", 5, "BSh", 6, "BSk", 7,
  "Csa", 8, "Csb", 9, "Csc", 10,
  "Cwa", 11, "Cwb", 12, "Cwc", 13,
  "Cfa", 14, "Cfb", 15, "Cfc", 16,
  "Dsa", 17, "Dsb", 18, "Dsc", 19, "Dsd", 20,
  "Dwa", 21, "Dwb", 22, "Dwc", 23, "Dwd", 24,
  "Dfa", 25, "Dfb", 26, "Dfc", 27, "Dfd", 28,
  "ET", 29, "EF", 30
)


df_IGBP <- tribble(
  ~IGBP, ~IGBP_id,
  "ENF", 1, "EBF", 2, "DNF", 3, 
  "DBF", 4, "MF", 5,  "CSH", 6,
  "OSH", 7, "WSA", 8, "SAV", 9, 
  "GRA", 10, "WET", 11, 
  "Croplands", 12, "Urban", 13, 
  "CropMosaic", 14, "Ice", 15, 
  "Barren", 16, "Water", 17, 
  "Water", 0)

