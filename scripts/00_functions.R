
# This function organizes different sources of raw flux data; updated on Jan, 2026
# It takes a few min to run this function
get_rawflux_file_structure <- function(dir_rawdata){
  
  files.unzip.Ameri <- list.files(file.path(dir_rawdata, "Ameriflux_FLUXNET", "unzip"), pattern =  "_FLUXNET_FULLSET_(HR|HH|DD|YY)", full.names = T)
  files.unzip.FLUXNET2020 <- list.files(file.path(dir_rawdata, "FLUXNET2020", "unzip"), pattern = "_FLUXNET2015_FULLSET_(HR|HH|DD|YY)", full.names = T)
  files.unzip.FLUXNET2015 <- list.files(file.path(dir_rawdata, "FLUXNET2015_10yr_only"), pattern = "", full.names = T)
  files.unzip.FLUXNET2020_ICOS <- list.files(file.path(dir_rawdata, "ICOS_FLUXNET2025-1", "unzip_connect2020"), full.names = T)
  files.unzip.FLUXNET2015_ICOS <- list.files(file.path(dir_rawdata, "ICOS_FLUXNET2025-1", "unzip_connect2015"), full.names = T)
  files.Ameri.ReddyProcGapFill <- list.files(file.path(dir_rawdata, "Ameriflux_BASE", "ReddyProcGapFill"), full.names=T)
  files.unzip.JapanFlux2024 <- list.files(file.path(dir_rawdata, "JapanFlux2024", "unzip"), full.names=T)
  files.csv.OzFlux2025 <- list.files(file.path(dir_rawdata, "OzFlux2025", "csv"), full.names=T)
  #
  df_file_structure <- data.frame(sources = c("AmeriFlux_FLUXNET", "AmeriFlux_BASE", 
                                              "FLUXNET2020_ICOS", "FLUXNET2015_ICOS", "FLUXNET2020", "FLUXNET2015", 
                                              "JapanFlux2024", "OzFlux2025"), 
                                  folders = I(list(list(files.unzip.Ameri), list(files.Ameri.ReddyProcGapFill), 
                                                   list(files.unzip.FLUXNET2020_ICOS), list(files.unzip.FLUXNET2015_ICOS),
                                                   list(files.unzip.FLUXNET2020), list(files.unzip.FLUXNET2015), 
                                                   list(files.unzip.JapanFlux2024), list(files.csv.OzFlux2025))), 
                                  HHfiles = c("_FLUXNET_FULLSET_(HR|HH)", "_EDDYPROC_FULLSET_HH",
                                              "_FLUXNET_HH", "_FLUXNET_HH", 
                                              "_FLUXNET2015_FULLSET_(HR|HH)", "_FLUXNET2015_FULLSET_(HR|HH)", 
                                              "_JapanFLUX2024_ALLVARS_HH", "_PyFluxPro_FULLSET_HH"), 
                                  DDfiles = c("_FLUXNET_FULLSET_DD", "_EDDYPROC_FULLSET_DD",
                                              "_FLUXNET_DD", "_FLUXNET_DD", 
                                              "_FLUXNET2015_FULLSET_DD", "_FLUXNET2015_FULLSET_DD", 
                                              "_JapanFLUX2024_ALLVARS_DD", "_PyFluxPro_FULLSET_DD"),
                                  YYfiles = c("_FLUXNET_FULLSET_YY", "_EDDYPROC_FULLSET_YY", 
                                              "_FLUXNET_YY", "_FLUXNET_YY", 
                                              "_FLUXNET2015_FULLSET_YY", "_FLUXNET2015_FULLSET_YY",
                                              "_JapanFLUX2024_ALLVARS_YY", "_PyFluxPro_FULLSET_YY"))
  return(df_file_structure)
}