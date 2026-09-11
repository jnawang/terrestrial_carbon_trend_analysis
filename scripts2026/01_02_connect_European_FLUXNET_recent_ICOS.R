# title: "Connect FLUXNET Shuttle, FLUXNET2015, and FLUXNET2025"
# output: pdf_document
# date: "2026-07-18"
# author: "Junna Wang"
#--------------------------------------------------------------------------------
# it will take about 1 hour to run this script. 

rm(list=ls())
# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'
dir_rawdata_origin <- '/Volumes/WZZ_disk/Thermal_Acclimation/SiteData'
dir_rawdata <- '/Volumes/WZZ_disk/Carbon_trend/SiteData'

library(librarian)
shelf(tidyverse, lubridate)

source(file.path("scripts2026", "00_functions.R"))

AmeriFlux <- read.csv(file.path('data', 'AmeriFlux-site-search-results-202512092232.csv'))
FLUXNET_shuttle <- read.csv(file.path('data2026', 'fluxnet_shuttle_snapshot_20260717T143751_selectedsites.csv'))
# excluding croplands
FLUXNET_shuttle <- FLUXNET_shuttle %>% filter(!igbp %in% c("CRO", "CVM", "BSV", "SNO"))


files_FLUXNET2015 <- list.files(path = file.path(dir_rawdata, "FLUXNET2015"), pattern = "\\.zip$")
FLUXNET2015_metadata <- data.frame(site_ID = substring(files_FLUXNET2015, 5, 10), 
                                   year_start = as.numeric(substring(files_FLUXNET2015, 32, 35)), 
                                   year_end = as.numeric(substring(files_FLUXNET2015, 37, 40)), stringsAsFactors = FALSE)

files_FLUXNET2025 <- list.files(path =  file.path(dir_rawdata_origin, "FLUXNET07202025"), pattern = "\\.zip$")
FLUXNET2025_metadata <- data.frame(site_ID = substring(files_FLUXNET2025, 5, 10),
                                   year_start = as.numeric(substring(files_FLUXNET2025, 35, 38)), 
                                   year_end = as.numeric(substring(files_FLUXNET2025, 40, 43)),
                                   stringsAsFactors = FALSE)

FLUXNET_all <- data.frame(site_ID = unique(c(FLUXNET_shuttle$site_id, FLUXNET2015_metadata$site_ID, FLUXNET2025_metadata$site_ID)),
                          year_start = NA, year_end = NA)
# only 4 sites in FLUXNET2015, but not in FLUXNET-Shuttle
setdiff(FLUXNET_all$site_ID, FLUXNET_shuttle$site_id)

for (i in 1:nrow(FLUXNET_all)) {
  site = FLUXNET_all$site_ID[i]
  year_start = 100000
  year_end = 0
  if (site %in% FLUXNET_shuttle$site_id) {
    year_start <- min(year_start, FLUXNET_shuttle$first_year[FLUXNET_shuttle$site_id == site])
    year_end <- max(year_end, FLUXNET_shuttle$last_year[FLUXNET_shuttle$site_id == site])
  }
  #
  if (site %in% FLUXNET2025_metadata$site_ID) {
    year_start <- min(year_start, FLUXNET2025_metadata$year_start[FLUXNET2025_metadata$site_ID == site])
    year_end <- max(year_end, FLUXNET2025_metadata$year_end[FLUXNET2025_metadata$site_ID == site])
  }
  #
  if (site %in% FLUXNET2015_metadata$site_ID) {
    year_start <- min(year_start, FLUXNET2015_metadata$year_start[FLUXNET2015_metadata$site_ID == site])
    year_end <- max(year_end, FLUXNET2015_metadata$year_end[FLUXNET2015_metadata$site_ID == site])
  }
  
  FLUXNET_all$year_start[i] <- year_start
  FLUXNET_all$year_end[i] <- year_end
}

# get sites longer than 10 years
FLUXNET_all_10yr <- FLUXNET_all %>% filter(year_end - year_start + 1 >= data_length_thresh) %>%
  mutate(year_start_shuttle = NA, year_end_shuttle = NA, 
         year_start_2025 = NA, year_end_2025 = NA,
         year_start_2015 = NA, year_end_2015 = NA)

# combine the three datasets
for (i in 1:nrow(FLUXNET_all_10yr)) {
  site = FLUXNET_all_10yr$site_ID[i]
  #
  missing_years <- FLUXNET_all_10yr$year_start[i]:FLUXNET_all_10yr$year_end[i]
  #
  if (site %in% FLUXNET_shuttle$site_id) {
    year_start <- FLUXNET_shuttle$first_year[FLUXNET_shuttle$site_id == site]
    year_end <- max(year_end, FLUXNET_shuttle$last_year[FLUXNET_shuttle$site_id == site])
    #
    FLUXNET_all_10yr$year_start_shuttle[i] <- min(intersect(missing_years, year_start:year_end))
    FLUXNET_all_10yr$year_end_shuttle[i] <- max(intersect(missing_years, year_start:year_end))
    #
    missing_years <- setdiff(missing_years, year_start:year_end)
    if (length(missing_years) == 0) { next }
  }
  #
  if (site %in% FLUXNET2025_metadata$site_ID) {
    year_start <- FLUXNET2025_metadata$year_start[FLUXNET2025_metadata$site_ID == site]
    year_end <- FLUXNET2025_metadata$year_end[FLUXNET2025_metadata$site_ID == site]
    #
    FLUXNET_all_10yr$year_start_2025[i] <- min(intersect(missing_years, year_start:year_end))
    FLUXNET_all_10yr$year_end_2025[i] <- max(intersect(missing_years, year_start:year_end))    
    #
    missing_years <- setdiff(missing_years, year_start:year_end)
    if (length(missing_years) == 0) { next }
  }
  #
  if (site %in% FLUXNET2015_metadata$site_ID) {
    year_start <- FLUXNET2015_metadata$year_start[FLUXNET2015_metadata$site_ID == site]
    year_end <- FLUXNET2015_metadata$year_end[FLUXNET2015_metadata$site_ID == site]
    #
    FLUXNET_all_10yr$year_start_2015[i] <- min(intersect(missing_years, year_start:year_end))
    FLUXNET_all_10yr$year_end_2015[i] <- max(intersect(missing_years, year_start:year_end))    
    #
    missing_years <- setdiff(missing_years, year_start:year_end)
    # print(missing_years)
  }
}

################
# create folder
if (sum(!is.na(FLUXNET_all_10yr$year_start_shuttle)) > 0) {
  if (!dir.exists(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_connect2015"))) {
    dir.create(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_connect2015"), recursive = TRUE)
  }
  #
  if (!dir.exists(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_connect2025"))) {
    dir.create(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_connect2025"), recursive = TRUE)
  } 
  #
  if (!dir.exists(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_select"))) {
    dir.create(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_select"), recursive = TRUE)
  } 
}

if (sum(!is.na(FLUXNET_all_10yr$year_start_2015) & 
        is.na(FLUXNET_all_10yr$year_start_shuttle) & 
        is.na(FLUXNET_all_10yr$year_start_2025)) > 0) {
  if (!dir.exists(file.path(dir_rawdata, "FLUXNET2015", "unzip_select"))) {
    dir.create(file.path(dir_rawdata, "FLUXNET2015", "unzip_select"), recursive = TRUE)
  } 
}

if (sum(!is.na(FLUXNET_all_10yr$year_start_2025) & 
        is.na(FLUXNET_all_10yr$year_start_shuttle) ) > 0) {
  if (!dir.exists(file.path(dir_rawdata, "FLUXNET2025", "unzip_select"))) {
    dir.create(file.path(dir_rawdata, "FLUXNET2025", "unzip_select"), recursive = TRUE)
  } 
}


#### get file names
unzip_files_FLUXNET_shuttle_HH <- list.files(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip"), pattern = "_FLUXNET_FLUXMET_(HH|HR)", full.names = TRUE)
unzip_files_FLUXNET_shuttle_DD <- list.files(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip"), pattern = "_FLUXNET_FLUXMET_(DD)", full.names = TRUE)  
unzip_files_FLUXNET_shuttle_YY <- list.files(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip"), pattern = "_FLUXNET_FLUXMET_(YY)", full.names = TRUE)
unzip_files_FLUXNET_shuttle <- c(unzip_files_FLUXNET_shuttle_HH, unzip_files_FLUXNET_shuttle_DD, unzip_files_FLUXNET_shuttle_YY)

unzip_files_FLUXNET2015_HH <- list.files(file.path(dir_rawdata, "FLUXNET2015", "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)", full.names = TRUE)
unzip_files_FLUXNET2015_DD <- list.files(file.path(dir_rawdata, "FLUXNET2015", "unzip"), pattern = "_FLUXNET2015_FULLSET_(DD)", full.names = TRUE)  
unzip_files_FLUXNET2015_YY <- list.files(file.path(dir_rawdata, "FLUXNET2015", "unzip"), pattern = "_FLUXNET2015_FULLSET_(YY)", full.names = TRUE)
unzip_files_FLUXNET2015 <- c(unzip_files_FLUXNET2015_HH, unzip_files_FLUXNET2015_DD, unzip_files_FLUXNET2015_YY)
#
unzip_files_FLUXNET2025_HH <- list.files(file.path(dir_rawdata_origin, "FLUXNET07202025", "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)", full.names = TRUE)
# unzip_files_FLUXNET2025_DD <- list.files(file.path(dir_rawdata_origin, "FLUXNET07202025", "unzip"), pattern = "_FLUXNET2015_FULLSET_DD", full.names = TRUE)  
# unzip_files_FLUXNET2025_YY <- list.files(file.path(dir_rawdata_origin, "FLUXNET07202025", "unzip"), pattern = "_FLUXNET2015_FULLSET_YY", full.names = TRUE)
# unzip_files_FLUXNET2025 <- c(unzip_files_FLUXNET2025_HH, unzip_files_FLUXNET2025_DD, unzip_files_FLUXNET2025_YY)
# we only have HH data for FLUXNET2025
unzip_files_FLUXNET2025 <- c(unzip_files_FLUXNET2025_HH)


# do the connection
for (i in 1:nrow(FLUXNET_all_10yr)) {
  site_ID <- FLUXNET_all_10yr$site_ID[i]
  
  print(i)
  # only if FLUXNET shuttle, copy three files into unzip_select
  # I have to standardize these variable names. 
  if (!is.na(FLUXNET_all_10yr$year_start_shuttle[i]) &
      is.na(FLUXNET_all_10yr$year_start_2025[i]) &
      is.na(FLUXNET_all_10yr$year_start_2015[i])) {
    #
    a_YY <- read.csv(unzip_files_FLUXNET_shuttle_YY[grepl(site_ID, unzip_files_FLUXNET_shuttle_YY)])
    a_YY[a_YY == -9999] <- NA
    #
    # only copy files when it has partitioned GPP and ER
    if (!"NEE_CUT_REF" %in% names(a_YY) && !"NEE_VUT_REF" %in% names(a_YY)) { next }
    if (!"GPP_DT_VUT_REF" %in% names(a_YY) && !"GPP_NT_VUT_REF" %in% names(a_YY) && 
        !"GPP_DT_CUT_REF" %in% names(a_YY) && !"GPP_NT_CUT_REF" %in% names(a_YY)) { next }
    
    if (!"NEE_VUT_REF" %in% names(a_YY) && "NEE_CUT_REF" %in% names(a_YY)) {
      a_YY$NEE_VUT_REF <- a_YY$NEE_CUT_REF
    }
    
    # deal with the five sites with complete YY FLUXNET data
    if (site_ID == "US-Rls") {
      a_YY_FLUXNET_old <- read.csv(file.path(dir_rawdata, "publication", "AMF_US-Rls_FLUXNET_FULLSET_YY_2014-2023_4-6.csv"))
    } else if (site_ID == "US-Rms") {
      a_YY_FLUXNET_old <- read.csv(file.path(dir_rawdata, "publication", "AMF_US-Rms_FLUXNET_FULLSET_YY_2014-2023_4-6.csv"))
    } else if (site_ID == 'US-Rws') {
      a_YY_FLUXNET_old <- read.csv(file.path(dir_rawdata, "publication", "AMF_US-Rws_FLUXNET_FULLSET_YY_2014-2023_4-6.csv"))
    } 
    if (site_ID %in% c("US-Rls", "US-Rms", 'US-Rws')) {
      a_YY_FLUXNET_old <- left_join(data.frame(TIMESTAMP = a_YY$TIMESTAMP), a_YY_FLUXNET_old, by = "TIMESTAMP")
      a_YY$NEE_VUT_REF[is.na(a_YY$NEE_VUT_REF)] <- a_YY_FLUXNET_old$NEE_VUT_REF[is.na(a_YY$NEE_VUT_REF)]
      a_YY[a_YY == -9999] <- NA
    }

    idx <- which(!is.na(a_YY$NEE_VUT_REF))
    
    if (length(idx) == 0 || max(idx) - min(idx) + 1 < 5) { next }

    # re-read the three files and adding NEE_VUT_REF, GPP_DT_VUT_REF, RECO_DT_VUT_REF, GPP_NT_VUT_REF, RECO_NT_VUT_REF
    a_YY <- read.csv(unzip_files_FLUXNET_shuttle_YY[grepl(site_ID, unzip_files_FLUXNET_shuttle_YY)])
    a_DD <- read.csv(unzip_files_FLUXNET_shuttle_DD[grepl(site_ID, unzip_files_FLUXNET_shuttle_DD)])
    a_HH <- read.csv(unzip_files_FLUXNET_shuttle_HH[grepl(site_ID, unzip_files_FLUXNET_shuttle_HH)])
    #
    if (!is_col_exist(a_YY, "NEE_VUT_REF") & is_col_exist(a_YY, "NEE_CUT_REF")) {
      a_YY$NEE_VUT_REF <- a_YY$NEE_CUT_REF
    }
    if (!is_col_exist(a_YY, "GPP_DT_VUT_REF") & is_col_exist(a_YY, "GPP_DT_CUT_REF")) {
      a_YY$GPP_DT_VUT_REF <- a_YY$GPP_DT_CUT_REF
      a_YY$RECO_DT_VUT_REF <- a_YY$RECO_DT_CUT_REF
    }
    if (!is_col_exist(a_YY, "GPP_NT_VUT_REF") & is_col_exist(a_YY, "GPP_NT_CUT_REF")) {
      a_YY$GPP_NT_VUT_REF <- a_YY$GPP_NT_CUT_REF
      a_YY$RECO_NT_VUT_REF <- a_YY$RECO_NT_CUT_REF      
    }
    
    #### update the three sites with earlier version of FLUXNET data
    if (site_ID %in% c("US-Rls", "US-Rms", 'US-Rws')) {
      a_YY$NEE_VUT_REF[a_YY$NEE_VUT_REF == -9999] <- a_YY_FLUXNET_old$NEE_VUT_REF[a_YY$NEE_VUT_REF == -9999]
      a_YY$NEE_VUT_REF_QC[a_YY$NEE_VUT_REF_QC == -9999] <- a_YY_FLUXNET_old$NEE_VUT_REF_QC[a_YY$NEE_VUT_REF_QC == -9999]
      if (is_col_exist(a_YY, "GPP_DT_VUT_REF") & is_col_exist(a_YY_FLUXNET_old, "GPP_DT_VUT_REF")) {
        a_YY$GPP_DT_VUT_REF[a_YY$GPP_DT_VUT_REF == -9999] <- a_YY_FLUXNET_old$GPP_DT_VUT_REF[a_YY$GPP_DT_VUT_REF == -9999]
        a_YY$RECO_DT_VUT_REF[a_YY$RECO_DT_VUT_REF == -9999] <- a_YY_FLUXNET_old$RECO_DT_VUT_REF[a_YY$RECO_DT_VUT_REF == -9999]
      }
      if (is_col_exist(a_YY, "GPP_NT_VUT_REF") & is_col_exist(a_YY_FLUXNET_old, "GPP_NT_VUT_REF")) {
        a_YY$GPP_NT_VUT_REF[a_YY$GPP_NT_VUT_REF == -9999] <- a_YY_FLUXNET_old$GPP_NT_VUT_REF[a_YY$GPP_NT_VUT_REF == -9999]
        a_YY$RECO_NT_VUT_REF[a_YY$RECO_NT_VUT_REF == -9999] <- a_YY_FLUXNET_old$RECO_NT_VUT_REF[a_YY$RECO_NT_VUT_REF == -9999]
      }
    }
    #####

    if (!is_col_exist(a_DD, "NEE_VUT_REF") & is_col_exist(a_DD, "NEE_CUT_REF")) {
      a_DD$NEE_VUT_REF <- a_DD$NEE_CUT_REF
    }    
    if (!is_col_exist(a_DD, "GPP_DT_VUT_REF") & is_col_exist(a_DD, "GPP_DT_CUT_REF")) {
      a_DD$GPP_DT_VUT_REF <- a_DD$GPP_DT_CUT_REF
      a_DD$RECO_DT_VUT_REF <- a_DD$RECO_DT_CUT_REF
    }   
    if (!is_col_exist(a_DD, "GPP_NT_VUT_REF") & is_col_exist(a_DD, "GPP_NT_CUT_REF")) {
      a_DD$GPP_NT_VUT_REF <- a_DD$GPP_NT_CUT_REF
      a_DD$RECO_NT_VUT_REF <- a_DD$RECO_NT_CUT_REF
    }    
    #
    if (!is_col_exist(a_HH, "NEE_VUT_REF") & is_col_exist(a_HH, "NEE_CUT_REF")) {
      a_HH$NEE_VUT_REF <- a_HH$NEE_CUT_REF
      a_HH$NEE_VUT_REF_QC <- a_HH$NEE_CUT_REF_QC
    }    
    if (!is_col_exist(a_HH, "GPP_DT_VUT_REF") & is_col_exist(a_HH, "GPP_DT_CUT_REF")) {
      a_HH$GPP_DT_VUT_REF <- a_HH$GPP_DT_CUT_REF
      a_HH$RECO_DT_VUT_REF <- a_HH$RECO_DT_CUT_REF
    }   
    if (!is_col_exist(a_HH, "GPP_NT_VUT_REF") & is_col_exist(a_HH, "GPP_NT_CUT_REF")) {
      a_HH$GPP_NT_VUT_REF <- a_HH$GPP_NT_CUT_REF
      a_HH$RECO_NT_VUT_REF <- a_HH$RECO_NT_CUT_REF
    }
    #
    file_name <- paste0('FLX_', site_ID, '_FLUXNET_FLUXMET_HH_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    var_name <- c("TIMESTAMP_START", "TIMESTAMP_END", "TA_F_MDS", "TA_F_MDS_QC", "TA_ERA", "TA_F", "TS_F_MDS_1", "TS_F_MDS_1_QC",	
                  "SWC_F_MDS_1", "SWC_F_MDS_1_QC", "VPD_F_MDS", "VPD_F_MDS_QC", "VPD_ERA", "VPD_F", "RH_F_MDS", "RH_F_MDS_QC", "RH", 
                  "SW_IN_F_MDS", "SW_IN_F_MDS_QC", "SW_IN_F", "SW_IN_ERA", "LE_F_MDS", "LE_F_MDS_QC", "P_ERA", "WS_F", 
                  "NEE_VUT_REF", "NEE_VUT_REF_QC", "GPP_NT_VUT_REF", "RECO_NT_VUT_REF", "GPP_DT_VUT_REF", "RECO_DT_VUT_REF", 
                  "CO2_F_MDS", "CO2_F_MDS_QC")
    write.csv(a_HH %>% select(any_of(var_name)), file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_select", file_name), row.names = FALSE)
    #
    file_name <- paste0('FLX_', site_ID, '_FLUXNET_FLUXMET_DD_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    var_name <- c("TIMESTAMP", "TA_F_MDS", "TA_F_MDS_QC", "TA_ERA", "TA_F", "TS_F_MDS_1", "TS_F_MDS_1_QC", 
                  "SWC_F_MDS_1", "SWC_F_MDS_1_QC", "VPD_F_MDS", "VPD_F_MDS_QC", "VPD_ERA", "VPD_F", "RH_F_MDS", "RH_F_MDS_QC", "RH", 
                  "SW_IN_F_MDS", "SW_IN_F_MDS_QC", "SW_IN_F", "SW_IN_ERA", "LE_F_MDS", "LE_F_MDS_QC", "P_ERA", 
                  "NEE_VUT_REF", "GPP_NT_VUT_REF", "RECO_NT_VUT_REF", "GPP_DT_VUT_REF", "RECO_DT_VUT_REF", "CO2_F_MDS")
    write.csv(a_DD %>% select(any_of(var_name)), file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_select", file_name), row.names = FALSE)
    #
    file_name <- paste0('FLX_', site_ID, '_FLUXNET_FLUXMET_YY_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    write.csv(a_YY %>% select(any_of(var_name)), file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_select", file_name), row.names = FALSE)
  }
  
  # only FLUXNET2025 into FLUXNET2025 folder, unzip_select
  if (is.na(FLUXNET_all_10yr$year_start_shuttle[i]) & 
      !is.na(FLUXNET_all_10yr$year_start_2025[i]) & 
      is.na(FLUXNET_all_10yr$year_start_2015[i])) {
    files2copy <- unzip_files_FLUXNET2025[grepl(site_ID, unzip_files_FLUXNET2025)]
    file.copy(from = files2copy, to = file.path(dir_rawdata, "FLUXNET2025", "unzip_select"))
  }  
  
  # only FLUXNET2015 into FLUXNET2015 folder, unzip_select
  if (is.na(FLUXNET_all_10yr$year_start_shuttle[i]) & 
      is.na(FLUXNET_all_10yr$year_start_2025[i]) & 
      !is.na(FLUXNET_all_10yr$year_start_2015[i])) {
    files2copy <- unzip_files_FLUXNET2015[grepl(site_ID, unzip_files_FLUXNET2015)]
    file.copy(from = files2copy, to = file.path(dir_rawdata, "FLUXNET2015", "unzip_select"))
  }  
  
  # FLUXNET2015 + FLUXNET2025 into FLUXNET2025 folder, unzip_select
  if (is.na(FLUXNET_all_10yr$year_start_shuttle[i]) & 
      !is.na(FLUXNET_all_10yr$year_start_2025[i]) & 
      !is.na(FLUXNET_all_10yr$year_start_2015[i])) {
    # we do not have this case. 
  }  
  
  # FLUXNET shuttle + FLUXNET2015 into FLUXNET shuttle, unzip_select
  if (!is.na(FLUXNET_all_10yr$year_start_shuttle[i]) & 
      is.na(FLUXNET_all_10yr$year_start_2025[i]) & 
      !is.na(FLUXNET_all_10yr$year_start_2015[i])) {
    #
    a_HH1 <- read.csv(unzip_files_FLUXNET2015_HH[grepl(site_ID, unzip_files_FLUXNET2015_HH)])
    a_HH2 <- read.csv(unzip_files_FLUXNET_shuttle_HH[grepl(site_ID, unzip_files_FLUXNET_shuttle_HH)])
    # filter proper years from a_HH1, 
    # mutate mean time, mutate year, filter year, and delete new columns
    dt <- ymd_hm(a_HH1$TIMESTAMP_END[1]) - ymd_hm(a_HH1$TIMESTAMP_START[1])
    a_HH1 <- a_HH1 %>% mutate(YEAR=year( ymd_hm(TIMESTAMP_START) + dt / 2) ) %>%
      filter(between(YEAR, FLUXNET_all_10yr$year_start_2015[i], FLUXNET_all_10yr$year_end_2015[i])) %>%
      select(-YEAR)
    #
    a_HH <- bind_rows(a_HH1, a_HH2)
    
    # calculate a_DD
    a_DD1 <- read.csv(unzip_files_FLUXNET2015_DD[grepl(site_ID, unzip_files_FLUXNET2015_DD)])
    a_DD2 <- read.csv(unzip_files_FLUXNET_shuttle_DD[grepl(site_ID, unzip_files_FLUXNET_shuttle_DD)])
    # filter years for a_DD1
    a_DD1 <- a_DD1 %>% mutate(YEAR = year(ymd(TIMESTAMP))) %>%
      filter(between(YEAR, FLUXNET_all_10yr$year_start_2015[i], FLUXNET_all_10yr$year_end_2015[i])) %>%
      select(-YEAR)
    a_DD <- bind_rows(a_DD1, a_DD2)
    
    a_YY1 <- read.csv(unzip_files_FLUXNET2015_YY[grepl(site_ID, unzip_files_FLUXNET2015_YY)])
    a_YY2 <- read.csv(unzip_files_FLUXNET_shuttle_YY[grepl(site_ID, unzip_files_FLUXNET_shuttle_YY)])
    # filter years for a_YY1
    a_YY1 <- a_YY1 %>% filter(between(TIMESTAMP, FLUXNET_all_10yr$year_start_2015[i], FLUXNET_all_10yr$year_end_2015[i]))
    a_YY <- bind_rows(a_YY1, a_YY2)
    
    # write the data out
    file_name <- paste0('ICOS_', site_ID, '_FLUXNET_FLUXMET_HH_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    write.csv(a_HH, file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_connect2015", file_name), row.names = FALSE)
    file_name <- paste0('ICOS_', site_ID, '_FLUXNET_FLUXMET_DD_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    write.csv(a_DD, file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_connect2015", file_name), row.names = FALSE)
    file_name <- paste0('ICOS_', site_ID, '_FLUXNET_FLUXMET_YY_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    write.csv(a_YY, file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_connect2015", file_name), row.names = FALSE)
  }
  
  #####connect FLUXNET2025 with FLUXNET shuttle
  if (!is.na(FLUXNET_all_10yr$year_start_shuttle[i]) & 
      !is.na(FLUXNET_all_10yr$year_start_2025[i]) & 
      is.na(FLUXNET_all_10yr$year_start_2015[i])) {
    #
    a_HH1 <- read.csv(unzip_files_FLUXNET2025_HH[grepl(site_ID, unzip_files_FLUXNET2025_HH)])
    a_HH2 <- read.csv(unzip_files_FLUXNET_shuttle_HH[grepl(site_ID, unzip_files_FLUXNET_shuttle_HH)])
    # filter proper years from a_HH1, 
    # mutate mean time, mutate year, filter year, and delete new columns
    dt <- ymd_hm(a_HH1$TIMESTAMP_END[1]) - ymd_hm(a_HH1$TIMESTAMP_START[1])
    a_HH1 <- a_HH1 %>% mutate(YEAR=year( ymd_hm(TIMESTAMP_START) + dt / 2) ) %>%
      filter(between(YEAR, FLUXNET_all_10yr$year_start_2025[i], FLUXNET_all_10yr$year_end_2025[i])) %>%
      select(-YEAR)
    #
    a_HH <- bind_rows(a_HH1, a_HH2)
    
    # estimate a_DD and a_YY data
    dt <- as.numeric(difftime(ymd_hm(a_HH1$TIMESTAMP_END[2]), ymd_hm(a_HH1$TIMESTAMP_END[1]), units=c('days')))
    a_HH1[a_HH1 == -9999] <- NA
    if (!"SWC_F_MDS_1" %in% names(a_HH1)) {
      a_HH1$SWC_F_MDS_1 <- NA
    }
    if (!"TS_F_MDS_1" %in% names(a_HH1)) {
      a_HH1$TS_F_MDS_1 <- NA
    }    
    
    a_DD1 <- a_HH1 %>% mutate(YEAR=year(ymd_hm(TIMESTAMP_START)), MONTH=month(ymd_hm(TIMESTAMP_START)),     
                              DAY=day(ymd_hm(TIMESTAMP_START)), TIMESTAMP=sprintf("%04d%02d%02d", YEAR, MONTH, DAY)) %>% 
      group_by(TIMESTAMP) %>% summarise(TA_F_MDS=mean(TA_F_MDS), VPD_F_MDS=mean(VPD_F_MDS), 
                                        SW_IN_F_MDS=mean(SW_IN_F_MDS), LE_F_MDS=mean(LE_F_MDS), 
                                        CO2_F_MDS=mean(CO2_F_MDS), P_ERA=sum(P_ERA),
                                        TA_ERA=mean(TA_ERA), VPD_ERA=mean(VPD_ERA), 
                                        TS_F_MDS_1=mean(TS_F_MDS_1), SWC_F_MDS_1=mean(SWC_F_MDS_1), 
                                        NEE_VUT_REF=sum(NEE_VUT_REF)*dt*0.0864*12,
                                        GPP_DT_VUT_REF=sum(GPP_DT_VUT_REF)*dt*0.0864*12,
                                        RECO_DT_VUT_REF=sum(RECO_DT_VUT_REF)*dt*0.0864*12, 
                                        GPP_NT_VUT_REF=sum(GPP_NT_VUT_REF)*dt*0.0864*12,
                                        RECO_NT_VUT_REF=sum(RECO_NT_VUT_REF)*dt*0.0864*12)
    # Yearly NEE, GPP, and RECO: gC m-2 year-1
    a_YY1 <- a_HH1 %>% mutate(TIMESTAMP=year(ymd_hm(TIMESTAMP_START))) %>% group_by(TIMESTAMP) %>% 
      summarise(TA_F_MDS=mean(TA_F_MDS), VPD_F_MDS=mean(VPD_F_MDS), 
                SW_IN_F_MDS=mean(SW_IN_F_MDS), LE_F_MDS=mean(LE_F_MDS), 
                CO2_F_MDS=mean(CO2_F_MDS), P_ERA=sum(P_ERA),
                TA_ERA=mean(TA_ERA), VPD_ERA=mean(VPD_ERA), 
                TS_F_MDS_1=mean(TS_F_MDS_1), SWC_F_MDS_1=mean(SWC_F_MDS_1), 
                NEE_VUT_REF=sum(NEE_VUT_REF)*dt*0.0864*12,
                GPP_DT_VUT_REF=sum(GPP_DT_VUT_REF)*dt*0.0864*12, 
                RECO_DT_VUT_REF=sum(RECO_DT_VUT_REF)*dt*0.0864*12,
                GPP_NT_VUT_REF=sum(GPP_NT_VUT_REF)*dt*0.0864*12,
                RECO_NT_VUT_REF=sum(RECO_NT_VUT_REF)*dt*0.0864*12)
    # Daytime and Nighttime NEE
    NEE_DAY <- a_HH1 %>% filter(NIGHT==0) %>% mutate(TIMESTAMP=year(ymd_hm(TIMESTAMP_START))) %>% group_by(TIMESTAMP) %>% 
      summarise(NEE_VUT_REF_DAY=sum(NEE_VUT_REF)*dt*0.0864*12)
    NEE_NIGHT <- a_HH1 %>% filter(NIGHT==1) %>% mutate(TIMESTAMP=year(ymd_hm(TIMESTAMP_START))) %>% group_by(TIMESTAMP) %>% 
      summarise(NEE_VUT_REF_NIGHT=sum(NEE_VUT_REF)*dt*0.0864*12)
    a_YY1 <- a_YY1 %>% left_join(NEE_DAY, by="TIMESTAMP") %>% left_join(NEE_NIGHT, by="TIMESTAMP")
    # plot the daily and yearly time series to ensure it makes sense. 
    plot(ymd_hm(a_HH1$TIMESTAMP_START), a_HH1$NEE_VUT_REF, main=site_ID)
    plot(ymd(a_DD1$TIMESTAMP), a_DD1$NEE_VUT_REF, main=site_ID)
    plot(a_YY1$TIMESTAMP, a_YY1$NEE_VUT_REF, main=site_ID)    
    
    #
    a_DD2 <- read.csv(unzip_files_FLUXNET_shuttle_DD[grepl(site_ID, unzip_files_FLUXNET_shuttle_DD)])
    a_DD2 <- a_DD2 %>% mutate(TIMESTAMP = as.character(TIMESTAMP))
    a_DD <- bind_rows(a_DD1, a_DD2)
    #
    a_YY2 <- read.csv(unzip_files_FLUXNET_shuttle_YY[grepl(site_ID, unzip_files_FLUXNET_shuttle_YY)])
    a_YY <- bind_rows(a_YY1, a_YY2)
    
    # write out results
    file_name <- paste0('ICOS_', site_ID, '_FLUXNET_FLUXMET_HH_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    write.csv(a_HH, file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_connect2025", file_name), row.names = FALSE)
    file_name <- paste0('ICOS_', site_ID, '_FLUXNET_FLUXMET_DD_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    write.csv(a_DD, file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_connect2025", file_name), row.names = FALSE)
    file_name <- paste0('ICOS_', site_ID, '_FLUXNET_FLUXMET_YY_', 
                        FLUXNET_all_10yr$year_start[i], '-', FLUXNET_all_10yr$year_end[i], '.csv')
    write.csv(a_YY, file.path(dir_rawdata, 'Fluxnet_shuttle', "unzip_connect2025", file_name), row.names = FALSE)
  }
}
